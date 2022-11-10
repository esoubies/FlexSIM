%--------------------------------------------------------------------------
% FlexSIM simulation script on simulated data
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;
if exist('Results.mat', 'file'), delete('Results.mat'); end % Delete Results.mat if exists
%% Parameters
% -- General parameters
params.pathToFlexSIM = '../';     % Path to the root of GitHub FlexSIM repo
params.GPU = 0;                   % Boolean on whether to use GPU or not
params.displ = 0;                 % Display option
imgType = 0;                      % Select 0 for synthethic, 1 for cameraman, 2 for   

% -- Properties of the SIM data stack
params.StackOrder= 'pa';          % Phase (p) and angle (a) convention. Choose one of ('paz', 'pza' or 'zap')
params.nbOr = 3;                  % Number of orientations
params.nbPh = 3;                  % Number of phases 

% -- OTF Approximation
params.lamb = 530;                % Emission wavelength
params.res = 64;                  % Pixel size (nm)
params.Na = 1.2;                  % Objective numerica aperture
params.damp = 1;                  % Damping factor for the OTF
params.nl=1.51;                   % Refractive index of the objective medium (glass/oil)
params.ns=1.333;                  % Refractive index of the sample medium (water)

% -- Parameters for patterns estimation
params.limits = [0.8, 1.05];      % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0, 1.1];    % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)

%% Main Loop
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
for MEP = [1000, 1000, 1000]
    params.MEP = MEP; 
    for a = [1 1 1]      
        params.a = a;
        disp(['<strong>## MEP = ',num2str(MEP),' a = ',num2str(a),'</strong>']);
        % -- Generate random orientations and phases
        params.DataPath = fullfile(pwd,sprintf('SIM_Simu_%d_%d.tif', imgType, params.MEP));    % Path to the SIM stack        
        params.ph = [0 pi/3 2*pi/3] + rand*pi/4 + 0.1 ; % Phases used for acquisition simulation
        params.or = [0 pi/3 2*pi/3] + rand*pi/3 + 0.1 ; % Orientationsused for acquisition simulation
        for or = 1:3
            params.k(or, :) = 2*pi*params.ns/params.lamb*[cos(params.or(or)), sin(params.or(or))]*params.Na/params.nl;
        end

        % -- Generate data & load data
        GenerateSIMData(imgType, params);        
        y = double(loadtiff(params.DataPath));
        if params.GPU
            y = gpuArray(x);
        end
        params.sz=size(y);
        [y, wf_stack] = OrderY(y, params); % Compute wiefield images                                   

        % -- Common quantities of interest
        [grids.I, grids.J] = meshgrid(0:params.sz(2)-1,0:params.sz(1)-1);              % Numerical mesh - multipurpose
        grids.X = grids.I*params.res; grids.Y=grids.J*params.res;                      % Scaled (by resolution) mesh 
        OTF = GenerateOTF(params.Na, params.lamb, params.sz, params.res, params.damp); % Computation of the OTF
        
        % Loop over batch of images (1 batch = 1 orr + x phases)
        imgIdxs = 1:params.nbPh*params.nbOr;    % Select the indexes to use (through imgIdxs)...
        imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]);
        OrientCount = 1; 
        FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
        for idx = imgIdxs
            disp(['- Orrientation #',num2str(OrientCount)]);
            % - Preprocessing Remove WF, Mask, and compute cross-corr
            wf = wf_stack(:,:,min(size(wf_stack,3),3));
            [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);
            fac=4; fftwf=fft2(padarray(wf,params.sz(1:2)*fac,'post')); fftG=fft2(padarray(G,params.sz(1:2)*fac,'post'));
            corrtmp=fftshift(ifft2((fft2(ifftshift(fftwf))).*conj(fft2(ifftshift(fftG)))));
            
            % - Eq-ph assumption
            params.method = 2;
            % Extract initial wavevector and phase (on grid using cross-coor) - done in this script to avoid repetition            
            wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
            tt=mean(corrtmp.*wght,3);tt2=conj(tt)./abs(tt);
            map=MaskFT(real(mean(corrtmp.*wght,3).*tt2),FCut,params.limits);
            [~,id]=max(map(:));[i,j]=ind2sub(size(map),id);
            kCCEqPh(OrientCount, :) = ([j,i]-floor(size(map)/2)-1)*pi/params.res./size(map);
            phaseCCEqPh(OrientCount, :) = mod(-angle(tt2(id)),2*pi)/2;       
            disp(['   [Eq-ph Assump] Error Init wavevect                  : ',num2str(norm(kCCEqPh(OrientCount, :)*sign(kCCEqPh(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
            % Refine wo filters
            kCCEqPhRef(OrientCount,:) = IterRefNoFilt(kCCEqPh(OrientCount, :),wf,G,grids,OTF,params.sz,params);
            phaseCCEqPhRef(OrientCount, :) =GetPhaseNoFilt(kCCEqPhRef(OrientCount, :)',wf,G,grids,OTF,params.sz,params);        
            disp(['   [Eq-ph Assump] Error Refined wavevect (w/o Filt)    : ',num2str(norm(kCCEqPhRef(OrientCount, :)*sign(kCCEqPhRef(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
            % Refine with filters
            kCCEqPhFilt(OrientCount,:) = IterRefinementWavevec(kCCEqPh(OrientCount, :)',wf,G,grids,OTF,params.sz,params);
            phaseCCEqPhFilt(OrientCount, :)=GetPhaseAndAmp(kCCEqPhFilt(OrientCount, :)',wf,G,grids,OTF,params.sz,params);
            disp(['   [Eq-ph Assump] Error Refined wavevect (Filt)        : ',num2str(norm(kCCEqPhFilt(OrientCount, :)*sign(kCCEqPhFilt(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
            
            % - No Eq-ph assumption 
            params.method = 1;
            % Extract initial wavevector and phase (on grid using cross-coor)
            tt=conj(corrtmp)./abs(corrtmp); 
            map=MaskFT(real(mean(corrtmp.*tt,3)),FCut,params.limits);
            [~,id]=max(map(:));[i,j]=ind2sub(size(map),id);
            kCC(OrientCount, :) = ([j,i]-floor(size(map)/2)-1)*pi/params.res./size(map);
            kCC(OrientCount, :) = sign(kCC(OrientCount, 1))*kCC(OrientCount, :);
            if all(sign([j,i]-size(map)/2)==sign(kCC)), sg=1; else sg=-1; end  % to know if we detected the one with same sign as simulated k (if not need to change the sign of the arg in the next line)
            phaseCC(OrientCount, :) = mod(-sg*angle(tt(i,j,:)),2*pi)/2; 
            disp(['   [No Eq-ph Assump] Error Init wavevect               : ',num2str(norm(kCC(OrientCount, :)*sign(kCC(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
            % Refine wo filters
            kCCRef(OrientCount,:) = IterRefNoFilt(kCC(OrientCount, :),wf,G,grids,OTF,params.sz,params);
            phaseCCRef(OrientCount, :)=GetPhaseNoFilt(kCCRef(OrientCount, :)',wf,G,grids,OTF,params.sz,params);   
            disp(['   [No Eq-ph Assump] Error Refined wavevect (w/o Filt) : ',num2str(norm(kCCRef(OrientCount, :)*sign(kCCRef(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
            % Refine with filters
            kCCFilt(OrientCount,:) = IterRefinementWavevec(kCC(OrientCount, :)',wf,G,grids,OTF,params.sz,params);
            phaseCCFilt(OrientCount, :)=GetPhaseAndAmp(kCCFilt(OrientCount, :)',wf,G,grids,OTF,params.sz,params);                           
            disp(['   [No Eq-ph Assump] Error Refined wavevect (Filt)     : ',num2str(norm(kCCFilt(OrientCount, :)*sign(kCCFilt(OrientCount,1))-params.k(OrientCount, :)*sign(params.k(OrientCount,1))))]);
    
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