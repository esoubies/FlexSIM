function [k_final, ph_final, k_init, ph_init] = EstimatePatterns(params,PosRoiPatt,y,k_init, wf_stack)
%--------------------------------------------------------------------------
% Function [k_final, ph_final, k_init, ph_init] = EstimatePatterns(params,PosRoiPatt,y,k_init, wf_stack)
%
% Estimate the parameters of 2D sinusoidal patterns of the form 
%    w(x) = 1 + a cos(<k,x> + phase),
% as described in [1].
%
% Inputs  : params  -> Structures with fields:
%                         - lamb: Emission wavelength
%                         - Na: Objective numerical aperture
%                         - res: resolution of the SIM data stac
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%                         - SzRoiPatt: Size of the ROI used for patterns estimation
%                         - eqPh: Boolean, if true equally-spaced phases are assumed
%                         - ringMaskLim: 1x2 array (eg. [0.3, 1.1]) defining the ring used to mask the WF component and the high-freq of the data, givien as factor of fc=1
%                         - nMinima: Number of starting points for the refinement steps
%                         - displ: if >1, displays intermediate results
%           PosRoiPatt -> Top-left corner of the ROI used for patterns estimation
%           y          -> SIM data stack
%           k_init     -> Initial wavevector. If provided, skips the grid search and runs only the
%                         iterative optimization. Set to 0 to perform grid search 
%           wf         -> Widefield image.                   
%
%
% Outputs : k_final  -> array with the estimated wavevectors 
%           ph_final -> array with the estimated phases 
%
% [1] Handling Challenging Structured Illumination Microscopy Data with FlexSIM
%     E. Soubies et al, Preprint, 2023
%
% See also FlexSIM.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Preprocessing of stack & initialization of useful variables
% Crop to the ROI, keeping the original size before
if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)   % Detect ROI and Crop to ROI
    tmp_y=zeros([params.SzRoiPatt,params.SzRoiPatt,size(y,[3,4])]);
    tmp_wf=zeros([params.SzRoiPatt,params.SzRoiPatt,size(wf_stack,[3,4])]);
    for ii=1:params.nframes
        tmp_y(:,:,:,ii) = y(PosRoiPatt(ii,1):PosRoiPatt(ii,1)+params.SzRoiPatt-1,PosRoiPatt(ii,2):PosRoiPatt(ii,2)+params.SzRoiPatt-1,:,ii);
        tmp_wf(:,:,:,ii) = wf_stack(PosRoiPatt(ii,1):PosRoiPatt(ii,1)+params.SzRoiPatt-1,PosRoiPatt(ii,2):PosRoiPatt(ii,2)+params.SzRoiPatt-1,:,ii);
    end
    y=tmp_y;wf_stack=tmp_wf;
end
sz = size(y);                                      % Calculate size of ROI
nt=size(y,4);                                      % Number of time steps
y = (y - min(y,[],1:3)) ./ (max(y,[],1:3) - min(y,[],1:3));     % Normalize stack images
wf_stack= (wf_stack-min(wf_stack,[],1:3)) ./ (max(wf_stack,[],1:3) - min(wf_stack,[],1:3));

if ~gather(k_init)
    compute_k_init=true;
    k_init=zeros(1,2);
else
    compute_k_init=false;
end
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
 
% Check the acquisition convention of the user and convert to ap(w)
imgIdxs = 1:params.nbPh*params.nbOr;    % Select the indexes to use (through imgIdxs)...
imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]);  
 
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

%% Preprocess
DispMsg(params.verbose,'      Preprocess: Remove WF and mask...');
G=zeros(size(y));wf=zeros(size(wf_stack));
OrientCount = 1;
for idx = imgIdxs
    idwf=min(size(wf_stack,3),OrientCount);
    [G(:,:,idx,:),wf(:,:,idwf,:)] = RemoveWFandMask(y(:,:,idx,:),wf_stack(:,:,idwf,:),params);
    OrientCount=OrientCount+1;
end

%% Initialization via cross-corr
if compute_k_init
    DispMsg(params.verbose,'      Initialization: Cross-corr btw WF and data in Fourier...');
    OrientCount = 1;
    for idx = imgIdxs        
        idwf=min(size(wf_stack,3),OrientCount);

        DispMsg(params.verbose,['        - Orientation #', num2str(OrientCount),': Cross-correl btw WF and data in Fourier...']);
        [map,K1,K2] = CrossCorr(G(:,:,idx,:),wf(:,:,idwf,:), params);
        Jmap=-MaskFT(mean(abs(map).^2,3:length(sz)),FCut,params.ringRegionSearch);

        DispMsg(params.verbose,['                          Extract initial wavevector...']);
        k_init(OrientCount, :)= ExtractLocMin(1,Jmap,K1,K2);

        if nargout==4  && compute_k_init % To output the initial phase estimate if required
            id=intersect(find(k_init(OrientCount,1)==K1(:)),find(k_init(OrientCount,2)==K2(:)));
            [ii,jj]=ind2sub(size(K1),id);
            ph_init(OrientCount, :,:)=mod(angle(map(ii,jj,:,:)),2*pi)/2;
        end

        OrientCount=OrientCount+1;
    end
end

% Displays
if params.displ > 0
    % - Displays related to estimated parameters
    fig_patt=DisplayPattParams(y(:,:,:,1),params,k_init,[],-1,0,'Initial estimated wavevectors');
end

%% Refinement
if nt==1, DispMsg(params.verbose,'      Refine initial wavevector and compute phases...');
else,  DispMsg(params.verbose,'      Refine initial wavevector and compute phases for each time step.'); end
OrientCount = 1; 
for idx = imgIdxs
    if nt==1,  DispMsg(params.verbose,['        - Orientation #', num2str(OrientCount),' ...']);
    else DispMsg(params.verbose,['        - Orientation #', num2str(OrientCount),', time step:'],0); end
    idwf=min(size(wf_stack,3),OrientCount);
    parfor (idt = 1:nt,params.nbcores)
        if nt >1,  DispMsg(params.verbose,[' t',num2str(idt)],0); end
        k_final(OrientCount, :,idt) = IterRefinementWavevec(k_init(OrientCount, :)',wf(:,:,idwf,idt),G(:,:,idx,idt),grids,OTF,sz(1:3),params);
        ph_final(OrientCount, :,idt)=GetPhaseAndAmp(k_final(OrientCount, :,idt)',wf(:,:,idwf,idt),G(:,:,idx,idt),grids,OTF,sz(1:3),params);
    end
    if nt >1 && params.verbose, fprintf('\n'); end

    OrientCount=OrientCount+1;
end

close(fig_patt);
end