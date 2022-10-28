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

warning('off','all')
time0=tic;
%% Routinary checks + Data loading
CheckParams(params);                       % Check conformity of parameters
y = double(loadtiff(params.DataPath));     % Read data
if params.GPU
    y = gpuArray(y);
end
sz=size(y);

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

% -- Displays
if params.displ > 0
    fig_rec=-1;  % Initialize figures
    DisplayStack(gather(y),'SIM Raw data',-1);
    leg={};
    if ~isempty(params.SzRoiBack)
        rectangle('Position',[PosRoiBack(2) PosRoiBack(1) params.SzRoiBack params.SzRoiBack],'EdgeColor','r');
        line(NaN,NaN,'Color','r'); leg{length(leg)+1}='ROI Background';% Hack for legend display of the rectangle
    end
    if ~isempty(params.SzRoiPatt)
        rectangle('Position',[PosRoiPatt(2) PosRoiPatt(1) params.SzRoiPatt params.SzRoiPatt],'EdgeColor','b');
        line(NaN,NaN,'Color','b'); leg{length(leg)+1}='ROI Patterns';% Hack for legend display of the rectangle
    end
    if ~isempty(leg), legend(leg); end
    drawnow;
end

%% FlexSIM pipeline
% -- Pattern Estimation
disp('<strong>=== Patterns estimation</strong> ...');
[k, phase] = EstimatePatterns(params, PosRoiPatt, y, 0, wf);
if params.estiPattLowFreq
    Lf = EstimateLowFreqPatterns(params,y,wf,5);
else
    Lf=1;
end

% Generate Patterns for reconstruction
patterns = GenerateReconstructionPatterns(params,PosRoiPatt,k,phase,params.pattAmp,sz,Lf);

% Displays
if params.displ > 0
    % - Display of estimated patterns
    DisplayStack(patterns,'Estimated Patterns',-1);
    % - Displays related to estimated parameters
    DisplayPattParams(y,params,k,phase,-1,0);
end
        
% -- Reconstruction
% Extract patches (if required)
if params.szPatch==0
    patches={y};
    patches_patt={patterns};
elseif params.szPatch>0
    patches=Image2Patches(y,params.szPatch,params.overlapPatch);
    patches_patt=Image2Patches(patterns,params.szPatch*2,params.overlapPatch*2);
end

% Initializations
nbPatches=length(patches(:));
rec=CellZeros(patches,[2,2],[1,2]);

% Loop over patches
if nbPatches==1 || ~params.parallelProcess
    % No parallelization over patches
    disp(['<strong>===  Reconstruction</strong> ...']);
    for id_patch = 1:nbPatches
        if params.szPatch>0, disp(['<strong>- [Patch #',num2str(id_patch),'/',num2str(nbPatches),'] </strong>']); end      
        rec{id_patch} = Reconstruct(gather(patches{id_patch}),gather(patches_patt{id_patch}),params);
        
        % Displays
        if params.displ > 0
            fig_rec=DisplayReconstruction(Patches2Image(rec,params.overlapPatch*2),wfUp,fig_rec);
        end
        
        % Save reconstruction 
        if params.sav
            prefix=params.DataPath(1:end-4);
            saveastiff(single(Patches2Image(rec,params.overlapPatch*2)),strcat(prefix,'_Rec.tif'));
        end
    end
else
    % Parallelization over patches
    disp('<strong>=== Reconstruction in parallel patch-based mode ... </strong>');
    if ~isempty(gcp('nocreate')), delete(gcp('nocreate')); end
    nbcores=feature('numcores');
    parpool('local',nbcores);
    parfevalOnAll(@warning,0,'off','all');
    parfor (id_patch = 1:nbPatches,nbcores)
        t = getCurrentTask(); 
        disp(['-- [Worker #',num2str(t.ID),'] Process patch #',num2str(id_patch),'/',num2str(nbPatches)]);
        rec{id_patch} = Reconstruct(gather(patches{id_patch}),gather(patches_patt{id_patch}),params);
    end
    delete(gcp('nocreate'));
end

%% Save
% - Save reconstruction / patterns / reconst. parameters / pattern parameters
if params.sav
    prefix=params.DataPath(1:end-4);
    saveastiff(single(Patches2Image(rec,params.overlapPatch*2)),strcat(prefix,'_Rec.tif'));
    saveastiff(single(patterns),strcat(prefix,'_Patt.tif'));
    save(strcat(prefix,'_Params'),'params');
end

% - Fill output variable
res.k = k; res.phase = phase; 
res.rec=Patches2Image(rec,params.overlapPatch*2);
res.patt=patterns;

disp(['<strong>=== FlexSIM END. Elapsed time (s): ',num2str(toc(time0)),' </strong>']);
end

