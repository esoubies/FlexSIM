%--------------------------------------------------------------------------
% FlexSIM simulation script on simulated data
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
params.pathToFlexSIM = '../';                      % Path to the root of GitHub FlexSIM repo

% -- Display, saving and GPU acceleration
params.GPU = 0;                         % Boolean on whether to use GPU or not
params.displ = 0;
rng(0);
%% Data related parameters
% -- Properties of the SIM data stack
params.StackOrder= 'pa';               % Phase (p) and angle (a) convention. Choose one of ('paz', 'pza' or 'zap')
params.SzRoiBack=[];                   % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)
params.nbOr = 3;                       % Number of orientations
params.nbPh = 3;                       % Number of phases 

% -- OTF Approximation, Acquisition parameters & Precomputations
params.lamb = 530;                     % Emission wavelength
params.res = 64;                       % Pixel size (nm)
params.Na = 1.2;                       % Objective numerica aperture
params.damp = 1;
params.nl=1.51;                        % Refractive index of the objective medium (glass/oil)
params.ns=1.333;                       % Refractive index of the sample medium (water)
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency

%% FlexSIM parameters
% -- Patch-based processing
params.szPatch=0;                 % If >0, FlexSIM will perform pattern estimation and reconstruction by patches of size 'szPatch'
params.overlapPatch=0;            % Overlap between patches if szPatch>0
params.parallelProcess=0;         % If 1, paralellizes the loop over patches

% -- Parameters for patterns estimation
params.SzRoiPatt = [];            % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.limits = [0.8, 1.05];      % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0, 1.1];    % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.nMinima = 1;               % Number of starting points for the refinement steps
params.nPoints = 0;             % Number of points in the J evaluation grid. If set to 0, initialization is done via peak detection
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)
params.method = 2;                % Method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns

% -- Path and files
imgType = 0;                       % Select 0 for synthethic, 1 for cameraman, 2 for        
for MEP = [10000, 10000, 10000]
    for a = [1 1 1]
        params.MEP = MEP; 
        params.a = a;
        
        params.DataPath = fullfile(pwd,sprintf('SIM_Simu_%d_%d.tif', imgType, params.MEP));    % Path to the SIM stack        
        params.ph = [0 pi/3 2*pi/3] + rand*pi/3; % Phases used for acquisition simulation
%         params.ph = [0 pi/3 2*pi/3];
        params.or = [0 pi/3 2*pi/3] + rand*pi/3; % Orientationsused for acquisition simulation
        for or = 1:3
            params.k(or, :) = 2*pi*params.ns/params.lamb*[cos(params.or(or)), sin(params.or(or))]*params.Na/params.nl;
        end
        params.k
        params.ph
        %% Generate data & load data
        run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
        GenerateSIMData(imgType, params);
        
        y = double(loadtiff(params.DataPath));
        if params.GPU
            y = gpuArray(x);
        end
        params.sz=size(y);
        %% Pre-processing
        % - Remove background
        if isfield(params,'SzRoiBack') && ~isempty(params.SzRoiBack)
            PosRoiBack=DetectPatch(sum(gather(y),3),params.SzRoiBack,-1);  % Detect the ROI for background
            y = RemoveBackground(y,PosRoiBack,params.SzRoiBack);           % Remove constant background and normalize in [0,1]
        else
            PosRoiBack=[1,1];
        end
        % - Detect ROI for pattern estimation
        if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)
            PosRoiPatt=DetectPatch(sum(gather(y),3),params.SzRoiPatt,1);
        else
            PosRoiPatt=[1,1];
        end
        
        % - Reorder stack with FlexSIM conventions
        [y, wf] = OrderY(y, params);                                    % Reorder and extract data if necessary
        wfUp=imresize(mean(wf,3),[size(wf,1),size(wf,2)]*2);            % For displays
        
        
        wf_stack = wf; 
        %%  Pattern Estimation
        if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)   % Detect ROI and Crop to ROI
            y = y(PosRoiPatt(1):PosRoiPatt(1)+params.SzRoiPatt-1,PosRoiPatt(2):PosRoiPatt(2)+params.SzRoiPatt-1,:);
            wf_stack = wf_stack(PosRoiPatt(1):PosRoiPatt(1)+params.SzRoiPatt-1,PosRoiPatt(2):PosRoiPatt(2)+params.SzRoiPatt-1,:);
        end
        sz = size(y);                                      % Calculate size of ROI
        y = (y - min(y(:))) / (max(y(:)) - min(y(:)));     % Normalize stack images
        wf_stack= (wf_stack-min(wf_stack(:))) / (max(wf_stack(:)) - min(wf_stack(:)));
         
        % Select the images of the same orientation 
        imgIdxs = 1:params.nbPh*params.nbOr;    % Select the indexes to use (through imgIdxs)...
        if params.method                        % according to the selected phase,z angle convention
            imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]);  
        end 
        
        % Common quantities of interest
        [grids.I, grids.J] = meshgrid(0:sz(2)-1,0:sz(1)-1);           % Numerical mesh - multipurpose
        grids.X = grids.I*params.res; grids.Y=grids.J*params.res;     % Scaled (by resolution) mesh 
        OTF = GenerateOTF(params.Na, params.lamb, sz, params.res, 1); % Computation of the OTF
        if params.GPU
            OTF = gpuArray(OTF); 
            grids.I = gpuArray(grids.I);
            grids.J = gpuArray(grids.J);
            grids.X = gpuArray(grids.X);
            grids.Y = gpuArray(grids.Y); 
        end
        
        %% Loop over batch of images (1 batch = 1 orr + x phases)
        OrientCount = 1; 
        for idx = imgIdxs
            % Preprocessing Remove WF and Mask
%             wf = wf_stack(:,:,min(size(wf_stack,3),3));
            wf = wf_stack(:,:,OrientCount);
            [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);
            
            % CC eq-ph - Extract wavevector and phase with standard method - done in this script to avoid repetition 
            params.method = 2;
            fac=4; fftwf=fft2(padarray(wf,sz(1:2)*fac,'post')); fftG=fft2(padarray(G,sz(1:2)*fac,'post'));
            corrtmp=fftshift(ifft2((fft2(ifftshift(fftwf))).*conj(fft2(ifftshift(fftG)))));
            wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
            tt=mean(corrtmp.*wght,3);tt2=conj(tt)./abs(tt);
            map=MaskFT(real(mean(corrtmp.*wght,3).*tt2),FCut,params.limits);
            [~,id]=max(abs(map(:)));[i,j]=ind2sub(size(map),id);
            kCCEqPh(OrientCount, :) = -([j,i]-floor(size(map)/2)-1)*pi/params.res./size(map);
            kCCEqPh(OrientCount, :) = sign(kCCEqPh(OrientCount, 1))*kCCEqPh(OrientCount, :)
            phaseCCEqPh(OrientCount, :) = mod(-angle(tt2(id)),2*pi)/2

            % Using the function CrossCorr (no phase extraction)
%             [map,K1,K2] = CrossCorr(G,wf, params);
%             map=-map; % As map corresponds here to cross-correl that we want to maximize (hence minimize the opposite)                    
%             kCCEqPh(OrientCount, :) = ExtractLocMin(params,map,K1,K2);            
%             fftGtmp = fft2(G(:,:,1)); phaseCCEqPh(OrientCount, :) = mod(angle(fftGtmp(id)),2*pi)/2;

%             if all(sign([j,i]-size(fftpattMask)/2)==sign(k)), sg=1; else sg=-1; end
%             phaseCCEqPh(OrientCount, :) = mod(sg*(angle(fftpattMask(id))),2*pi)/2;
%             phaseCCEqPh(OrientCount, :)=GetPhaseAndAmp(kInit(OrientCount, :)',wf,G,grids,OTF,sz,params);
        
            % Refine wo filters
            kCCEqPhRef(OrientCount,:) = IterRefNoFilt(kCCEqPh(OrientCount, :),wf,G,grids,OTF,sz,params);
            phaseCCEqPhRef(OrientCount, :) =GetPhaseNoFilt(kCCEqPhRef(OrientCount, :),wf,G,grids,OTF,sz,params);
            
            % Refine (from init) with filters (Use standard function)
            kCCEqPhFilt(OrientCount,:) = IterRefinementWavevec(kCCEqPh(OrientCount, :)',wf,G,grids,OTF,sz,params);
            phaseCCEqPhFilt(OrientCount, :)=GetPhaseAndAmp(kCCEqPhFilt(OrientCount, :)',wf,G,grids,OTF,sz,params);

%             if kCCEqPhFilt(OrientCount, 1) >  1
%                 kCCEqPhRef(OrientCount,:) = IterRefNoFilt(kCCEqPh(OrientCount, :),wf,G,grids,OTF,sz,params);
%             end


            % CC (generalization) - Extract wavevector and phase with standard method
            params.method = 1;
            % Use initializations from CC eq-ph
            tt=conj(corrtmp)./abs(corrtmp); 
            map=MaskFT(real(mean(corrtmp.*tt,3)),FCut,params.limits);
%             [map,K1,K2] = CrossCorr(G,wf, params);
%             map=-map; % As map corresponds here to cross-correl that we want to maximize (hence minimize the opposite)                    
%             kCC(OrientCount, :)= ExtractLocMin(params,map,K1,K2);
            [~,id]=max(abs(map(:)));[i,j]=ind2sub(size(map),id);
            kCC(OrientCount, :) = ([j,i]-floor(size(map)/2)-1)*pi/params.res./size(map);
            kCC(OrientCount, :) = sign(kCC(OrientCount, 1))*kCC(OrientCount, :);
            if all(sign([j,i]-size(map)/2)==sign(kCC)), sg=1; else sg=-1; end  % to know if we detected the one with same sign as simulated k (if not need to change the sign of the arg in the next line)
            phaseCC(OrientCount, :) = mod(-sg*angle(tt(i,j,:)),2*pi)/2;            
%             fftGtmp = fft2(G(:,:,1)); phaseCC(OrientCount, :) = mod(angle(fftGtmp(id)),2*pi)/2;
%             if all(sign([j,i]-size(fftpattMask)/2)==sign(k)), sg=1; else sg=-1; end
%             phaseCC(OrientCount, :) = mod(sg*(angle(fftpattMask(id))),2*pi)/2;
%             phaseCC(OrientCount, :)=GetPhaseAndAmp(kCC(OrientCount, :)',wf,G,grids,OTF,sz,params);
        
            % Refine wo filters
            kCCRef(OrientCount,:) = IterRefNoFilt(kCC(OrientCount, :),wf,G,grids,OTF,sz,params);
            phaseCCRef(OrientCount, :)=GetPhaseNoFilt(kCCRef(OrientCount, :)',wf,G,grids,OTF,sz,params);           
            
            % Refine (from init) with filters (Use standard function)
            kCCFilt(OrientCount,:) = IterRefinementWavevec(kCC(OrientCount, :)',wf,G,grids,OTF,sz,params);
            phaseCCFilt(OrientCount, :)=GetPhaseAndAmp(kCCFilt(OrientCount, :)',wf,G,grids,OTF,sz,params);                          
        
            % Estiamte phases from the refined (with filters) wavevector estimation
    %         phArg(OrientCount) = mod((angle(...)),2*pi)/2;   % Simple arg{wavevecto} issue: FT doesn't have enough resolution
%             phJ(OrientCount, :)=GetPhaseNoFilt(kFilt(OrientCount,:),wf,G,grids,OTF,sz,params);
%             phFilt(OrientCount, :)=GetPhaseAndAmp(kFilt(OrientCount,:),wf,G,grids,OTF,sz,params);
            OrientCount=OrientCount+1;
        end
        % Include estimated params in WV
        params.kCCEqPh = kCCEqPh ; params.phCCEqPh = phaseCCEqPh;
        params.kCCEqPhRef = kCCEqPhRef; params.phCCEqPhRef = phaseCCEqPhRef; 
        params.kCCEqPhFilt = kCCEqPhFilt; params.phCCEqPhFilt = phaseCCEqPhFilt; 
        params.kCC = kCC; params.phCC = phaseCC; 
        params.kCCRef = kCCRef; params.phCCRef = phaseCCRef; 
        params.kCCFilt = kCCFilt; params.phCCFilt = phaseCCFilt; 

        clear kCCEqPh kCCEqPhRef kCCEqPhFilt kCC kCCRef kCCFilt 

        params = EvalRun(params);                            % Process the results
    end
end

% Compare with CrossCorr for comp wo refinement
% Discussion points - 
% CC (method 1) vs CCEqPh as a generalization
% Improvement with the refinement
% Improvement with the filter
% Download 2 data sets from https://www.nature.com/articles/s41377-021-00513-w#Sec15
% (Same link) Test the code from the same paper
% (Examples/SIM-...) Add a ReadMe explaining the conditions of the image, parameters, comparisons
% For plot, use boxplots summarizing all the MEP and contrast combinations (at the minimum, ensure the same combinations across different methods)


% 1. k
% 1. PH
% 1. FrobeNorm

function k = IterRefNoFilt(k_init,wf,G,grids,OTF,sz,params)
% -- Line search by backtracking parameters
tau_fact=1.5; tau_min=1e-10; nit_tot=100; cf_tolerance = 1e-6; tau_init=1e-2;  count_it = 1;    
k=k_init;    

% Calculate current cost
cf(count_it) = EvalJ(k, wf, G, params, grids, 0, 0, 0); %#ok<AGROW>    
% Gradient descent w.r.t. wavevector
for jj=1:nit_tot        
    tau = tau_init;
    [~, g] = EvalJ(k, wf, G, params, grids, 0, 0, 1);    % Get current gradient
    ktmp = k - tau * g;                                           % Update k
    cf_new = EvalJ(ktmp, wf, G, params, grids, 0, 0, 0); % Calculate new cost
    % Line search by backtracking
    if cf(count_it) > cf_new
        while cf(count_it) > cf_new
            k = ktmp;
            tau = tau * tau_fact;
            ktmp = k - tau * g;
            cf_new = EvalJ(ktmp, wf, G, params, grids, 0, 0, 0);
        end
        ktmp = k; tau = tau/tau_fact;
    else
        while cf(count_it)<=cf_new && (tau>tau_min)
            tau = tau/tau_fact;
            ktmp = k - tau * g;
            cf_new = EvalJ(ktmp, wf, G, params, grids, 0, 0, 0);
        end
        if (tau>tau_min)
            k=ktmp;
        end
    end
    % Stop grad descent if step size becomes negligible
    if (tau<tau_min) 
        break;
    end
    
    % Calculate and store new cost
    count_it=count_it+1;                  
    cf(count_it) = EvalJ(ktmp, wf, G, params, grids, 0, 0, 0);    
    % Stop grad descent if cost improvement becomes negligible
    if abs(cf(count_it)-cf(count_it-1))/abs(cf(count_it-1))<cf_tolerance
        break;
    end

end
end

function ph=GetPhaseNoFilt(k,wf,G,grids,OTF,sz,params)
% -- Build matrix A
A=BuildA(k, wf, 0, params, grids); AA = A'*A;
% -- Solve linear system
if params.method == 2
    s = AA\A'*G(:);
elseif params.method == 1
    s = zeros(2, params.nbPh);
    for phNb = 1:params.nbPh
        G_tmp= G(:,:,phNb);
        s(:,phNb) = AA\A'*G_tmp(:);
    end
else
    G_tmp= G_filt;
    s = AA\A'*G_tmp(:);
end
% -- Extract cos and sin components
ac=s(1,:); as=s(2,:);

% -- Convert (cos,sin) --> angle
tmp = atan(as./ac);           % Calculate phase, ac, and as
tmp = tmp + pi * (ac < 0);
ph=mod(tmp,2*pi)/2;   
end