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
if params.GPU, y = gpuArray(y); end
sz=size(y);
tmp = fft_best_dim(sz(1)*2+params.padSz) - sz(1)*2; 
if tmp~=params.padSz
    DispMsg(params.verbose,['Note: padSz has been increased from ',num2str(params.padSz),' to ',num2str(tmp),' for faster FFTs (multiple of powers of 2, 3 and/or 5)']);
    params.padSz=tmp;
end


%% Pre-processing
% - Remove background
if isfield(params,'SzRoiBack') && ~isempty(params.SzRoiBack)
    PosRoiBack=DetectPatch(mean(gather(y),[3,4]),params.SzRoiBack,-1);  % Detect the ROI for background
    y = RemoveBackground(y,PosRoiBack,params.SzRoiBack);           % Remove constant background and normalize in [0,1]
else
    PosRoiBack=[1,1];
end
% - Detect ROI for pattern estimation
if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)
    if ~isfield(params,'posRoiPatt') || isempty(params.posRoiPatt)
        PosRoiPatt=DetectPatch(mean(gather(y),[3,4]),params.SzRoiPatt,1);
    else
        PosRoiPatt=params.posRoiPatt;
    end
else
    PosRoiPatt=[1,1];
end


% - Reorder stack with FlexSIM conventions
[y, wf] = OrderYandExtWF(y, params);                            % Reorder and extract data if necessary
if isfield(params,'frameRange') && ~isempty(params.frameRange)  % Keep only frames specified by user
    y=y(:,:,:,params.frameRange);             
    wf=wf(:,:,:,params.frameRange);
end
params.nframes=size(y,4);                                       % Number of time frames
wfUp=imresize(mean(wf,3),[size(wf,1),size(wf,2)]*2);            % For displays

% -- Displays
if params.displ > 0
    fig_rec=-1;  % Initialize figures
    DisplayStack(gather(y(:,:,:,1)),'SIM Raw data',-1);
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
% -- Estimate sinusoidal patterns components
DispMsg(params.verbose,'<strong>=== Patterns estimation</strong> ...');
DispMsg(params.verbose,'---> Estimate sinusoidal patterns components ...');
[k, phase] = EstimatePatterns(params, PosRoiPatt, y, 0, wf);
if  isfield(params,'timeAvgPattParams') && params.timeAvgPattParams
    k=mean(k,3);phase=mean(phase,3);
end

% -- Loop over frames
rec=zeros([sz(1:2)*2,params.nframes]);
figStackPatt=-1;figDispPatt=-1;fig_rec=-1;
for it=1:params.nframes   
    if params.nframes>1, DispMsg(params.verbose,['<strong>=== Process temporal frame #',num2str(it),'/',num2str(params.nframes),'</strong> ...']); end

    % Estimate low frequency patterns components
    if params.estiPattLowFreq
        DispMsg(params.verbose,'---> Estimate low frequency patterns components ...');
        Lf = EstimateLowFreqPatterns(params,y(:,:,:,it),wf(:,:,:,it),5);
    else
        Lf=1;
    end
    % Generate Patterns for reconstruction
    if params.estiPattLowFreq && params.padSz>0
        Lf= padarray(Lf,[1,1]*params.padSz,'replicate','post');
    end
    DispMsg(params.verbose,'---> Generate reconstruction patterns ...');
    patterns = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,min(it,size(k,3))),phase(:,:,min(it,size(phase,3))),params.pattAmp,sz+params.padSz/2,Lf);

    % Displays
    if params.displ > 0 && it==1
        % - Display of estimated patterns
        figStackPatt=DisplayStack(patterns,'Estimated Patterns',figStackPatt);
        % - Displays related to estimated parameters
        figDispPatt=DisplayPattParams(y(:,:,:,it),params,k(:,:,min(it,size(k,3))),phase(:,:,min(it,size(phase,3))),figDispPatt,0);
    end

    % -- Reconstruction
    if params.nframes==1,  DispMsg(params.verbose,'<strong>===  Reconstruction</strong> ...'); 
    else, if params.sepOrr==0, DispMsg(params.verbose,'---> Reconstruct ...'); end; end
    
    if params.szPatch==0
        rec(:,:,it) = Reconstruct(gather(y(:,:,:,it)),gather(patterns),params);
    elseif params.szPatch>0
        % Extract patches (if required)
        patches=Image2Patches(y(:,:,:,it),params.szPatch,params.overlapPatch);
        patches_patt=Image2Patches(patterns,params.szPatch*2+params.padSz,params.overlapPatch*2+params.padSz);

        % Initializations
        nbPatches=length(patches(:));
        rec_tmp=CellZeros(patches,[2,2],[1,2]);
        % Loop over patches
        if params.parallelProcess && (params.szPatch>0)
            nbcores=feature('numcores');
            parfevalOnAll(@warning,0,'off','all');
        else
            nbcores=0;
        end
        parfor (id_patch = 1:nbPatches,nbcores)
            if nbcores>0
                t = getCurrentTask();
                DispMsg(params.verbose,['-- [Worker #',num2str(t.ID),'] Process patch #',num2str(id_patch),'/',num2str(nbPatches)]);
            else
                if params.szPatch>0, fprintf('<strong>- [Process patch #%i/%i]</strong> ...',id_patch,nbPatches); end
                if params.verbose==2 && nbPatches>1, fprintf('\n'); end
            end

            rec_tmp{id_patch} = Reconstruct(gather(patches{id_patch}),gather(patches_patt{id_patch}),params);
        end

        rec(:,:,it)=Patches2Image(rec_tmp,params.overlapPatch*2);
    end

    % Displays
    if params.displ > 0
        fig_rec=DisplayReconstruction(rec(:,:,it),wfUp(:,:,:,it),fig_rec);
    end

    %% Save
    % - Save reconstruction / patterns / reconst. parameters / pattern parameters
    if params.sav
        prefix=params.DataPath(1:end-4);
        saveastiff(rec,strcat(prefix,'_Rec.tif'));
        saveastiff(patterns,strcat(prefix,'_Patt.tif'));
        save(strcat(prefix,'_Params'),'params');
    end
    DispMsg(params.verbose,['== FlexSIM Elapsed time (s): ',num2str(toc(time0))]);
end

% - Fill output variable
res.k = k; res.phase = phase; 
res.rec=rec;
res.patt=patterns;
DispMsg(params.verbose,['<strong>=== FlexSIM END. Elapsed time (s): ',num2str(toc(time0)),' </strong>']);


