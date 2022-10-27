function [k_final, phase, a] = EstimatePatterns(params,PosRoiPatt,y,k_init, wf_stack)
%--------------------------------------------------------------------------
% Function [k_final, phase, a] = EstimatePatterns(params,PosRoiPatt,y,k_init, wf_stack)
%
% Estimate the parameters of 2D sinusoidal patterns of the form 
%    w(x) = 1 + a cos(<k,x> + phase),
% as described in [1].
%
% Inputs  : params  -> Structures with fields:
%                         - StackOrder: order of the SIM stack. Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
%                         - roi: region on interest on which the parameters will be estimated
%                         - lamb: Emission wavelength
%                         - Na: Objective numerical aperture
%                         - res: resolution of the SIM data stac
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%                         - SzRoiPatt: Size of the ROI used for patterns estimation
%                         - method: method used to estimate the parameters, 3 choices
%                                    0 - treat all images independently
%                                    1 - use all images with same orientation to estimate a unique wavevector
%                                    2 - 1 + assume equally spaced phases
%                         - ringMaskLim: 1x2 array (eg. [0.3, 1.1]) defining the ring used to mask the WF component and the high-freq of the data, givien as factor of fc=1
%                         - limits: 1x2 array (eg. [0.9, 1.1]) defining  the ring over which the J function is evaluated for initializing, givien as factor of fc=1
%                         - nMinima: Number of starting points for the refinement steps
%                         - nPoints: Number of points in the J evaluation grid
%                         - FilterRefinement: Number of times that the filter is upgraded (gradient descent cycles)
%                         - displ: if >1, displays intermediate results
%           PosRoiPatt -> Top-left corner of the ROI used for patterns estimation
%           y          -> SIM data stack
%           k_init     -> Initial wavevector. If provided, skips the grid search and runs only the
%                         iterative optimization. Set to 0 to perform grid search 
%           wf         -> Widefield image.                   
%
%
% Outputs : k_final -> array with the estimated wavevectors (if method=0, of size nbOr*nbPh. Otherwise, of size nbOr)
%           phase   -> array with the estimated phases (if method=0,1, of size nbOr*nbPh. Otherwise, of size nbOr)
%           a       -> Same size as phase, contains the amplitude of the SIM patterns
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also FlexSIM.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Preprocessing of stack & initialization of useful variables
% Crop to the ROI, keeping the original size before
if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)   % Detect ROI and Crop to ROI
    y = y(PosRoiPatt(1):PosRoiPatt(1)+params.SzRoiPatt-1,PosRoiPatt(2):PosRoiPatt(2)+params.SzRoiPatt-1,:);
    wf_stack = wf_stack(PosRoiPatt(1):PosRoiPatt(1)+params.SzRoiPatt-1,PosRoiPatt(2):PosRoiPatt(2)+params.SzRoiPatt-1,:);
end
sz = size(y);                                      % Calculate size of ROI
y = (y - min(y(:))) / (max(y(:)) - min(y(:)));     % Normalize stack images
wf_stack= (wf_stack-min(wf_stack(:))) / (max(wf_stack(:)) - min(wf_stack(:)));

if ~gather(k_init)
    compute_k_init=true;
else
    compute_k_init=false;
end
 
% Check the acquisition convention of the user and convert to ap(w)
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
    DispMsg(1,[' Batch of images: ', num2str(idx')]);     % Display info to the user
    DispMsg(1,'   - Remove WF and mask...');
    wf = wf_stack(:,:,min(size(wf_stack,3),3));
    [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);

    if params.displ > 1                                 % Debug mode - display conditioned images
        f = figure; 
        if params.method
            sgtitle(sprintf("Conditioned data for orientation #%d", OrientCount))
        else
            sgtitle(sprintf("Conditioned data for image #%", OrientCount))
        end
        WFPan = axes('position',[0.55,0.55,0.4,0.4], 'title','Widefield', 'Parent', f);
        WFFTPan = axes('position',[0.55,0.05,0.4,0.4], 'title','Widefield FT', 'Parent', f);
        if numel(size(G)) > 2            
            ImgPan = uipanel('position',[0.05,0.55,0.4,0.4], 'title','SIM Images', 'Parent', f);
            ImgFTPan = uipanel('position',[0.05,0.05,0.4,0.4], 'title','SIM Images FT', 'Parent', f);
            sliceViewer(gather(G), "Colormap",viridis, "Parent",ImgPan);
            sliceViewer(gather(log10(abs(fftshift(fft2(G))+1))), "Colormap",viridis, "Parent",ImgFTPan);                                                                  
            imshow(gather(wf), [], 'Parent', WFPan); colormap(viridis); title('WF');
            imshow(gather(log10(abs(fftshift(fft2(wf))+1))), [], 'Parent', WFFTPan); colormap(viridis); title('WF FT')
        else
            subplot(2, 2, 1); imshow(G, [], "Parent",ImgPan); colormap(viridis);
            subplot(2, 2, 3); imshow(log10(abs(fftshift(fft2(G))+1)), [], "Parent",ImgFTPan); colormap(viridis);
            subplot(2, 2, 2); imshow(wf, [], 'Parent', WFPan); colormap(viridis);
            subplot(2, 2, 4); imshow(log10(abs(fftshift(fft2(G))+1)), [], 'Parent', WFFTPan); colormap(viridis);
        end               
    end
    
    if compute_k_init
        if params.nPoints
            DispMsg(1,'   - Grid-based evaluation of J landscape...');
            [map,K1,K2] = GridEvalJ(params,wf,G,grids);
            
            if params.displ > 1                           % If requested, display J grid
                mJp=max(map(:));
                fg=figure; subplot(1,2,1); axis xy;hold on; view(45,45);
                surf(K1,K2,map,'FaceColor','interp','EdgeColor','interp');
                colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14); %axis([0 maxp -maxp maxp]);grid;
                sgtitle(sprintf('J landscape for orientation #%d', OrientCount))
                title('Landscape and initial local minima');
                subplot(1,2,2); axis xy;hold on; title('Zoom');
                surf(K1,K2,map,'FaceColor','interp','EdgeColor','interp');
                colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14)
                drawnow;
            end
        else
            DispMsg(1,'   - Cross-correl btw WF and data in Fourier...');
            [map,K1,K2] = CrossCorr(G,wf, params);
            map=-map; % As map corresponds here to cross-correl that we want to maximize (hence minimize the opposite)            
        end
        
        DispMsg(1,['   - Extracting ',num2str(params.nMinima),' candidate wave-vectors...']);
        k_init= ExtractLocMin(params,map,K1,K2);
                    
        DispMsg(1,'   - Refine position of candidate wave-vectors...');
        fprintf('%s','     - candidate #');
        for ithk = 1:params.nMinima
            if mod(ithk,ceil(params.nMinima/10))==0
                if ithk==params.nMinima, fprintf('%i\n',ithk);  else, fprintf('%i, ',ithk); end
            end
            k_init(ithk,:) = IterRefinementWavevec(k_init(ithk, :)',wf,G,grids,OTF,sz,params);
        end
        
        DispMsg(1,'   - Choosing the best wavevector...');    % Choose the best wavevector in terms of value of J
        if params.GPU
            Jp=zeros(1,params.nMinima,'double','gpuArray');
        else
            Jp=zeros(1,params.nMinima);
        end
        for iii = 1:params.nMinima
            Jp(iii) = EvalJ(k_init(iii,:), wf, G, params, grids, 0, 0, 0);
        end
        [~,optIdx] = min(Jp);
        k_final(OrientCount, :) = k_init(optIdx, :);
    else
        DispMsg(1,'   - Refine position of wavevector...');
        k_final(OrientCount, :) = IterRefinementWavevec(k_init(OrientCount, :)',wf,G,grids,OTF,sz,params);
    end
    
    DispMsg(1,'   - Computing phases and amplitutes...'); 
    [phase(OrientCount, :),a(OrientCount, :)]=GetPhaseAndAmp(k_final(OrientCount, :)',wf,G,grids,OTF,sz,params);
    
    OrientCount=OrientCount+1;
end
end