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

warning off backtrace
time0=tic;
%% Routinary checks + Data loading
CheckParams(params);                       % Check conformity of parameters
y = double(loadtiff(params.DataPath));     % Read data
if params.GPU
    y = gpuArray(y);
end

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
    fig_y=-1;fig_patt=-1;fig_patt_par=-1;fig_rec=-1;  % Initialize figures
    fig_y=DisplayStack(gather(y),'SIM Raw data',fig_y);
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
% -- Extract patches (if required)
if params.szPatch==0
    patches={y};
    patches_wf={wf};
elseif params.szPatch>0
    patches=Image2Patches(y,params.szPatch,params.overlapPatch);
    patches_wf=Image2Patches(wf,params.szPatch,params.overlapPatch);
end

% -- Initializations
nbPatches=length(patches(:));
patterns=CellZeros(patches,[2,2,1],[1,2,3]);
rec=CellZeros(patches,[2,2],[1,2]);
Lf=CellZeros(patches,[2,2,1],[1,2,3]);

% -- Loop over patches
if nbPatches==1 || ~params.parallelProcess
    % No parallelization over patches
    prefix_disp='';
    for id_patch = 1:nbPatches
        if params.szPatch>0, prefix_disp=['[Patch #',num2str(id_patch),'/',num2str(nbPatches),']']; end
        sz_p=size(patches{id_patch});                             % Patch size
        
        % -- Pattern Estimation
        disp(['<strong>=== ',prefix_disp,' Patterns estimation</strong> ...']);
        
        [k(:,:,id_patch), phase(:,id_patch), a] = EstimatePatterns(params, PosRoiPatt, patches{id_patch}, 0, patches_wf{id_patch});
        if params.estiPattLowFreq
            Lf{id_patch} = EstimateLowFreqPatterns(params,patches{id_patch},patches_wf{id_patch},5);
        else
            Lf{id_patch}=1;
        end
        
        % -- Generate Patterns for reconstruction
        a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
        patterns{id_patch} = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,id_patch),phase(:,id_patch),a,sz_p,Lf{id_patch});
        
        % -- Displays
        if params.displ > 0
            % - Display of estimated patterns
            fig_patt=DisplayStack(Patches2Image(patterns,params.overlapPatch*2),'Estimated Patterns',fig_patt);
            % - Displays related to estimated parameters
            if params.szPatch==0
                fig_patt_par=DisplayPattParams(patches{id_patch},params,k(:,:,id_patch),phase(:,id_patch),a,fig_patt_par,0);
            else
                fig_patt_par=DisplayPattParams(patches{id_patch},params,k(:,:,id_patch),phase(:,id_patch),a,fig_patt_par,id_patch);
            end
        end
        
        % -- Reconstruction
        disp(['<strong>=== ',prefix_disp,' Reconstruction</strong> ...']);
        rec{id_patch} = Reconstruct(gather(patches{id_patch}),gather(patterns{id_patch}),params);
        
        % -- Displays
        if params.displ > 0
            fig_rec=DisplayReconstruction(Patches2Image(rec,params.overlapPatch*2),wfUp,fig_rec);
        end
        
        % - Save reconstruction / patterns
        if params.sav
            prefix=params.DataPath(1:end-4);
            saveastiff(single(Patches2Image(rec,params.overlapPatch*2)),strcat(prefix,'_Rec.tif'));
            saveastiff(single(Patches2Image(patterns,params.overlapPatch*2)),strcat(prefix,'_Patt.tif'));
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
            disp(['<strong>--- ',prefix_disp,' Patterns correction</strong> ...']);
            [k(:,:,ii), phase(:,id_patch), a] = EstimatePatterns(params,PosRoiPatt, patches{ii},kmed, patches_wf{ii});
            a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
            patterns{ii} = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,ii),phase(:,id_patch),a,sz_p,Lf{ii});
            disp(['<strong>--- ',prefix_disp,' New reconstruction</strong> ...']);
            rec{ii} = Reconstruct(gather(patches{ii}),gather(patterns{ii}),params);
            % -- Displays
            if params.displ >0
                fig_patt_par=DisplayPattParams(patches{ii},params,k(:,:,ii),phase(:,id_patch),a,fig_patt_par,ii);
                fig_patt=DisplayStack(Patches2Image(patterns,params.overlapPatch*2),'Estimated Patterns',fig_patt);
                fig_rec=DisplayReconstruction(Patches2Image(rec,params.overlapPatch*2),wfUp,fig_rec);
            end
        end
    end
else
    % Parallelization over patches
    disp('<strong>START FlexSIM in parallel patch-based mode ... </strong>');
    parfor id_patch = 1:nbPatches
        t = getCurrentTask(); 
        disp(['-- [Worker #',num2str(t.ID),'] Process patch #',num2str(id_patch),'/',num2str(nbPatches)]);
        sz_p=size(patches{id_patch});                             % Patch size
        
        % -- Pattern Estimation        
        [k(:,:,id_patch), phase(:,id_patch), a] = EstimatePatterns(params, PosRoiPatt, patches{id_patch}, 0, patches_wf{id_patch});
        if params.estiPattLowFreq
            Lf{id_patch} = EstimateLowFreqPatterns(params,patches{id_patch},patches_wf{id_patch},5);
        else
            Lf{id_patch}=1;
        end
        
        % -- Generate Patterns for reconstruction
        a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
        patterns{id_patch} = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,id_patch),phase(:,id_patch),a,sz_p,Lf{id_patch});
   
        % -- Reconstruction
        rec{id_patch} = Reconstruct(gather(patches{id_patch}),gather(patterns{id_patch}),params);
    end
    
    % -- Re-run reconstruction for patches with wrong patterns
    kmed=median(k,3); % median wavevector
    relErr=sum(sum((k-kmed).^2,1),2)./sum(sum((k).^2,1),2);
    idx_err=find(relErr>1e-3);
    parfor idx=1:nbPatches
        if sum(idx_err==idx)==1
            t = getCurrentTask();
            disp(['-- [Worker #',num2str(t.ID),'] Correct patch #',num2str(idx)]);
            sz_p=size(patches{idx});
            [k(:,:,idx), phase(:,idx), a] = EstimatePatterns(params,PosRoiPatt, patches{idx},kmed, patches_wf{idx});
            a=a./a; % TODO: Hardcode to 1 for now (to be as in previous version)
            patterns{idx} = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,idx),phase(:,idx),a,sz_p,Lf{idx});
            rec{idx} = Reconstruct(gather(patches{idx}),gather(patterns{idx}),params);
        end
    end  
end

%% Save
% - Save reconstruction / patterns / reconst. parameters / pattern parameters
if params.sav
    prefix=params.DataPath(1:end-4);
    saveastiff(single(Patches2Image(rec,params.overlapPatch*2)),strcat(prefix,'_Rec.tif'));
    saveastiff(single(Patches2Image(patterns,params.overlapPatch*2)),strcat(prefix,'_Patt.tif'));
    save(strcat(prefix,'_Params'),'params');
end

% - Fill output variable
res.k = k; res.phase = phase; 
res.rec=Patches2Image(rec,params.overlapPatch*2);
res.patt=Patches2Image(patterns,params.overlapPatch*2);

disp(['<strong>=== FlexSIM END. Elapsed time (s): ',num2str(toc(time0)),' </strong>']);
end

