function FlexSIM(params)
%--------------------------------------------------------------------------
% Function  FlexSIM(params)
% 
% FlexSIM [1] should be called from a script defining the structure `params`,
% which  includes all the paths and optical parameters necessary for a 
% reconstruction. (See folder Examples)
%
% Inputs : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
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

% -- Set stuff for parallel processing
if params.parallelProcess
    params.nbcores=feature('numcores');
    parfevalOnAll(@warning,0,'off','all');
    params.paraLoopOrr=((params.sepOrr) && (params.nframes==1));
    params.paraLoopFrames=(params.nframes>1);
else
    params.nbcores=0;params.paraLoopOrr=0;params.paraLoopFrames=0;
end


% -- Displays
if params.displ > 0
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
    k=repmat(mean(k,3),[1 1 params.nframes]);phase=repmat(mean(phase,3),[1 1 params.nframes]);
end

% Displays
if params.displ > 0
    % - Displays related to estimated parameters
    DisplayPattParams(y(:,:,:,1),params,k(:,:,1),phase(:,:,1),-1,0);
end

% -- Loop over frames
rec=zeros([sz(1:2)*2,params.nframes]);
if params.paraLoopFrames, DispMsg(params.verbose,'<strong>===  Reconstruction</strong> ...'); end
parfor (it = 1:params.nframes,params.nbcores*params.paraLoopFrames)
    if params.nframes>1
        if ~params.paraLoopFrames
            DispMsg(params.verbose,['<strong>=== Process temporal frame #',num2str(it),'/',num2str(params.nframes),'</strong> ...']);
        else
            t = getCurrentTask();
            DispMsg(params.verbose,['-- [Worker #',num2str(t.ID),'] Process temporal frame #',num2str(it),'/',num2str(params.nframes),' ...']);
        end
    end
    % Estimate low frequency patterns components
    if params.estiPattLowFreq
        if ~params.paraLoopFrames, DispMsg(params.verbose,'---> Estimate low frequency patterns components ...'); end
        Lf = EstimateLowFreqPatterns(params,y(:,:,:,it),wf(:,:,:,it),5);
    else
        Lf=1;
    end
    % Generate Patterns for reconstruction
    if params.estiPattLowFreq && params.padSz>0
        Lf= padarray(Lf,[1,1]*params.padSz,'replicate','post');
    end
    if ~params.paraLoopFrames, DispMsg(params.verbose,'---> Generate reconstruction patterns ...'); end
    patterns = GenerateReconstructionPatterns(params,PosRoiPatt,k(:,:,it),phase(:,:,it),params.pattAmp,sz+params.padSz/2,Lf);    

    % -- Reconstruction
    if params.nframes==1,  DispMsg(params.verbose,'<strong>===  Reconstruction</strong> ...'); 
    else, if params.sepOrr==0 && ~params.paraLoopFrames, DispMsg(params.verbose,'---> Reconstruct ...'); end; end
    rec(:,:,it) = Reconstruct(gather(y(:,:,:,it)),gather(patterns),params);

    % - Save 
    if params.sav && (params.nframes>1)
        if params.nframes>1
            saveastiff(rec(:,:,it),[params.DataPath(1:end-4),'_Rec_Frame_',num2str(it),'.tif']);
            saveastiff(patterns,[params.DataPath(1:end-4),'_Patt_Frame_',num2str(it),'.tif']);
        else
            saveastiff(patterns,[params.DataPath(1:end-4),'_Patt.tif']);
        end
    end
    if params.nframes>1
        if ~params.paraLoopFrames
            DispMsg(params.verbose,['== FlexSIM Elapsed time (s): ',num2str(toc(time0))]);
        else
            DispMsg(params.verbose,['-- [Worker #',num2str(t.ID),'] Process temporal frame #',num2str(it),'/',num2str(params.nframes),' done. FlexSIM Elapsed time (s): ',num2str(toc(time0))]);
        end
    end
end

% Displays
if params.displ > 0
    DisplayReconstruction(rec(:,:,1),wfUp(:,:,:,1),-1);
end

% Save 
if params.sav
    save(strcat(params.DataPath(1:end-4),'_Params'),'params'); 
    saveastiff(rec,strcat(params.DataPath(1:end-4),'_Rec.tif'));
    if params.nframes>1
        patt=zeros([sz(1:2)*2+params.padSz,params.nbOr*params.nbPh,params.nframes]);
        for it=1:params.nframes
            delete([params.DataPath(1:end-4),'_Rec_Frame_',num2str(it),'.tif']);
            patt(:,:,:,it)=loadtiff([params.DataPath(1:end-4),'_Patt_Frame_',num2str(it),'.tif']);
            delete([params.DataPath(1:end-4),'_Patt_Frame_',num2str(it),'.tif']);
        end
        saveastiff(reshape(patt,[size(patt,[1,2]),prod(size(patt,[3,4]))]),[params.DataPath(1:end-4),'_Patt.tif']);
    end
end
DispMsg(params.verbose,['<strong>=== FlexSIM END. Elapsed time (s): ',num2str(toc(time0)),' </strong>']);


