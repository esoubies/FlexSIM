function res = FlexSIM(params)
%--------------------------------------------------------------------------
% Function res = FlexSIM(params)
% 
% FlexSIM [1] should be called from a script defining the structure `params`,
% which  includes all the paths and optical parameters necessary for a 
% reconstruction. (See folder Examples)
%
% Inputs : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Outputs: res    -> Structure with the final image (in the field `res`)
%                    and other intermediate results, like patterns.  
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Data loading + routinary checks
y = double(loadtiff(params.DataPath));    % Read data 
CheckParams(params);                      % Check conformity of parameters
wf=mean(y,3);wf=imresize(wf,size(wf)*2);  % Widefield image;

% -- Displays
if params.displ > 0
    fig_y=-1;fig_patt=-1;fig_patt_par=-1;fig_rec=-1;  % Initialize figures
    fig_y=DisplayStack(y,'SIM Raw data',fig_y);
end

%% FlexSIM pipeline
% -- Extract patches (if required)
norm_patches=[];
if params.szPatch==0
    patches={y};
elseif params.szPatch>0 && params.enhanceContrast
    patches=Image2Patches(y,params.szPatch,params.overlapPatch); 
else
    [patches,norm_patches]=Image2Patches(y,params.szPatch,params.overlapPatch);
end

% -- Initializations
nbPatches=length(patches(:));
patterns=CellZeros(patches,[2,2,1],[1,2,3]);
rec=CellZeros(patches,[2,2],[1,2]);
% if params.method                       % To think about - I would skip the...
%     k=zeros(params.nbOr,2,nbPatches);  % initialization altogether, it's not looping that much
% else
%     k=zeros(params.nbOr*params.nbPh,2,nbPatches);

% -- Loop over patches
prefix_disp='';
for id_patch = 1:nbPatches
    if params.szPatch>0, prefix_disp=['[Patch #',num2str(id_patch),'/',num2str(nbPatches),']']; end
    sz_p=size(patches{id_patch});                             % Patch size
    
    % -- Pattern Estimation
    disp(['<strong>=== ',prefix_disp,' Patterns parameter estimation START</strong> ...']);
    
    [k(:,:,id_patch), phase, a] = EstimatePatterns(params, patches{id_patch});
    if params.estiPattLowFreq
        Lf = EstimateLowFreqPatterns(patches{id_patch},5);
    end
    
    % -- Generate Patterns for reconstruction
%     a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
    patterns{id_patch} = GenerateReconstructionPatterns(params,k(:,:,id_patch),phase,a,sz_p);
    
    % -- Displays
    if params.displ > 0
        % - Display of estimated patterns
        fig_patt=DisplayStack(Patches2Image(patterns,params.overlapPatch*2),'Estimated Patterns',fig_patt);
        % - Displays related to estimated parameters
        if params.szPatch==0
            fig_patt_par=DisplayPattParams(patches{id_patch},params,k(:,:,id_patch),phase,a,fig_patt_par,0);
        else
            fig_patt_par=DisplayPattParams(patches{id_patch},params,k(:,:,id_patch),phase,a,fig_patt_par,id_patch);
        end
    end
    
    % -- Reconstruction
    disp(['<strong>=== ',prefix_disp,' Reconstruction START</strong> ...']);
    rec{id_patch} = Reconstruct(patches{id_patch},patterns{id_patch},params);
    
    % -- Displays
    if params.displ > 0
        fig_rec=DisplayReconstruction(Patches2Image(rec,params.overlapPatch*2,norm_patches),wf,fig_rec);
    end
end

% -- Re-run reconstruction for patches with wrong patterns
if params.szPatch>0
    disp('<strong>=== Detect and correct patches with bad estimated patterns</strong> ...');
    kmed=median(k,3); % median wavevector
    relErr=sum(sum((k-kmed).^2,1),2)./sum(sum((k).^2,1),2);
    idx_err=find(relErr>1e-3);
    for ii=idx_err'
        prefix_disp=['[Correction Patch #',num2str(ii),']'];
        sz_p=size(patches{ii});    
        disp(['<strong>--- ',prefix_disp,' Patterns parameter correction  START</strong> ...']);
        [k(:,:,ii), phase, a] = EstimatePatterns(params, patches{ii},kmed);
        a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
        patterns{ii} = GenerateReconstructionPatterns(params,k(:,:,ii),phase,a,sz_p);
        disp(['<strong>--- ',prefix_disp,' New reconstruction START</strong> ...']);
        rec{ii} = Reconstruct(patches{ii},patterns{ii},params);
        % -- Displays
        if params.displ >0
            fig_patt_par=DisplayPattParams(patches{ii},params,k(:,:,ii),phase,a,fig_patt_par,ii);
            fig_patt=DisplayStack(Patches2Image(patterns,params.overlapPatch*2),'Estimated Patterns',fig_patt);
            fig_rec=DisplayReconstruction(Patches2Image(rec,params.overlapPatch*2,norm_patches),wf,fig_rec);
        end
    end
end

%% Save
% - Save reconstruction / patterns / reconst. parameters / pattern parameters
if params.sav
    prefix=params.DataPath(1:end-4);
    saveastiff(single(Patches2Image(rec,params.overlapPatch*2,norm_patches)),strcat(prefix,'_Rec.tif'));
    saveastiff(single(Patches2Image(patterns,params.overlapPatch*2)),strcat(prefix,'_Patt.tif'));
    save(strcat(prefix,'_Params'),'params');
end

% - Fill output variable
res.rec=Patches2Image(rec,params.overlapPatch*2,norm_patches);
res.patt=Patches2Image(patterns,params.overlapPatch*2);

end

