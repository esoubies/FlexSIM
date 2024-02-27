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
% [1] Handling Challenging Structured Illumination Microscopy Data with FlexSIM
%     E. Soubies et al, Preprint, 2023
%
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

warning('off','all')
time0=tic;
if params.GPU
    useGPU(1);
else
    useGPU(0);
end
%% Routinary checks + Data loading
CheckParams(params);                       % Check conformity of parameters
if strcmp(params.DataPath(end-3:end),'.tif')
    y = double(loadtiff(params.DataPath));     % Read data
elseif strcmp(params.DataPath(end-3:end),'.nd2')
    if isfield(params,'channel') && ~isempty(params.channel)
        y = nd2read(params.DataPath,params.channel);
    else
        y = nd2read(params.DataPath);
    end
end
% Set default parameters related to the treatment of temporal stacks
if ~isfield(params,'frameRange'), params.frameRange=[]; end
if ~isfield(params,'cstTimePatt'), params.cstTimePatt=0; end
if ~isfield(params,'framePattEsti'), params.framePattEsti=[]; end
if ~isfield(params,'doRefinement'), params.doRefinement=1; end
if ~isempty(params.framePattEsti), assert(params.cstTimePatt==1,'If framePattEsti can be non empty only when cstTimePatt is true.'); end
if params.cstTimePatt, assert(params.eqPh==1,'Currently cstTimePatt=1 is possible only with eqPh=1.'); end

%% Pre-processing
% - Reorder stack with FlexSIM conventions
[y, wf] = OrderYandExtWF(y, params);                            % Reorder and extract data if necessary
if ~isempty(params.frameRange)  % Keep only frames specified by user
    assert(all((params.frameRange>=1).*(params.frameRange<=size(y,4))),'Parameter frameRange is not compatible with the size of the data.')
    y=y(:,:,:,params.frameRange);             
    if ~isempty(wf), wf=wf(:,:,:,params.frameRange); end
end
params.nframes=size(y,4);                                       % Number of time frames
% - Crop if required
if isfield(params,'SzRoi') && ~isempty(params.SzRoi)
    if ~isfield(params,'posRoi') || isempty(params.posRoi)
        PosRoi=DetectPatch(mean(y,[3,4]),params.SzRoi,1);
    else
        PosRoi=params.posRoi;
    end
    y=y(PosRoi(1):PosRoi(1)+params.SzRoi-1,PosRoi(2):PosRoi(2)+params.SzRoi-1,:,:);
end
% - Get size 
sz=size(y);
tmp = fft_best_dim(sz(1)*2+params.padSz) - sz(1)*2; 
if tmp~=params.padSz
    DispMsg(params.verbose,['Note: padSz has been increased from ',num2str(params.padSz),' to ',num2str(tmp),' for faster FFTs (multiple of powers of 2, 3 and/or 5)']);
    params.padSz=tmp;
end

% - Remove background
if isfield(params,'SzRoiBack') && ~isempty(params.SzRoiBack)
    for ii=params.nframes:-1:1 % Inverse order so that the last PosRoiBack corresponds to the first frame for display below
        PosRoiBack=DetectPatch(mean(y(:,:,:,ii),3),params.SzRoiBack,-1);  % Detect the ROI for background
        y(:,:,:,ii) = RemoveBackground(y(:,:,:,ii),PosRoiBack,params.SzRoiBack);           % Remove constant background and normalize in [0,1]
        if ~isempty(wf),  wf(:,:,:,ii) = RemoveBackground(wf(:,:,:,ii),PosRoiBack,params.SzRoiBack); end
    end
else
    PosRoiBack=[1,1];
end
if isempty(wf)
    for ii=1:params.nbOr
        wf(:,:,ii,:)=mean(y(:,:,(ii-1)*params.nbPh+1:ii*params.nbPh,:),3);
    end
end
wfUp=imresize(mean(wf,3),[size(wf,1),size(wf,2)]*2);            % For displays
% - Detect ROI for pattern estimation
if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)
    if ~isfield(params,'posRoiPatt') || isempty(params.posRoiPatt)
        PosRoiPatt=zeros(params.nframes,2);
        for ii=1:params.nframes
            PosRoiPatt(ii,:)=DetectPatch(mean(y(:,:,:,ii),3),params.SzRoiPatt,1);
        end
    else
        PosRoiPatt=repmat(params.posRoiPatt,[params.nframes,1]);
    end
else
    PosRoiPatt=ones([params.nframes,1]);
end

% -- Set stuff for parallel processing
if params.parallelProcess
    params.nbcores=feature('numcores');
    parfevalOnAll(@warning,0,'off','all');
    params.paraLoopOrr=((params.sepOrr) && (params.nframes==1));
    params.paraLoopFrames=(params.nframes>1);
else
    params.nbcores=0;params.paraLoopOrr=0;params.paraLoopFrames=0;
    if ~isempty(gcp('nocreate'))
       delete(gcp('nocreate'));
    end
end

% -- Displays
if params.displ > 0
    id=1;
    DisplayStack(y(:,:,:,id),'SIM Raw data',-1);
    leg={};
    if ~isempty(params.SzRoiBack)
        rectangle('Position',[PosRoiBack(2) PosRoiBack(1) params.SzRoiBack params.SzRoiBack],'EdgeColor','r');
        line(NaN,NaN,'Color','r'); leg{length(leg)+1}='ROI Background';% Hack for legend display of the rectangle
    end
    if ~isempty(params.SzRoiPatt)
        rectangle('Position',[PosRoiPatt(id,2) PosRoiPatt(id,1) params.SzRoiPatt params.SzRoiPatt],'EdgeColor','b');
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
if params.cstTimePatt
    med=median(phase,3);tt=cat(2,abs(phase-med),abs(phase+pi -med),abs(phase-pi-med));
    [~,id]=min(tt,[],2);
    phase_min=phase.*(id==1)+(phase+pi).*(id==2)+(phase-pi).*(id==3);
    phase=repmat(mean(phase_min,3),[1 1 params.nframes]);
    if params.displ > 0
        cmap=[0 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
        if isempty(params.framePattEsti), tt=1:params.nframes; else tt=params.framePattEsti; end 
        figure;
        subplot(1,2,1); set(gca,'fontsize',14);hold all;
        leg={};
        for ii=1:params.nbOr
            plot(tt,reshape(k(ii,1,:),1,length(tt))','*','linewidth',2,'markersize',10,'color',cmap(ii,:));
            plot(tt,reshape(k(ii,2,:),1,length(tt))','o','linewidth',2,'markersize',10,'color',cmap(ii,:));
            plot(tt,repmat(mean(k(ii,1,:),3),[length(tt), 1]),'--','linewidth',2,'color',cmap(ii,:));
            plot(tt,repmat(mean(k(ii,2,:),3),[length(tt), 1]),'-.','linewidth',2,'color',cmap(ii,:));
            leg{end+1}=['Or #',num2str(ii),' kx'];
            leg{end+1}=['Or #',num2str(ii),' ky'];
            leg{end+1}=['Or #',num2str(ii),' mean kx'];
            leg{end+1}=['Or #',num2str(ii),' mean ky'];
        end
        xlabel('Frames');title('Wavevector');legend(leg(:));grid;
        subplot(1,2,2); set(gca,'fontsize',14);hold all; title('Phases');
        leg={};
        for ii=1:params.nbOr
            plot(tt,reshape(phase_min(ii,1,:),1,length(tt))','*','linewidth',2,'markersize',10,'color',cmap(ii,:));hold all;
            plot(tt,repmat(mean(phase_min(ii,1,:),3),[length(tt), 1]),'--','linewidth',2,'color',cmap(ii,:));
            leg{end+1}=['Or #',num2str(ii),' ph-offset'];
            leg{end+1}=['Or #',num2str(ii),' mean'];
        end
        xlabel('Frames'); ylabel('Phases');legend(leg(:));grid;
    end    
    k=repmat(mean(k,3),[1 1 params.nframes]);
end

% Displays
if params.displ > 0
    % - Displays related to estimated parameters
    DisplayPattParams(y(:,:,:,1),params,k(:,:,1),phase(:,:,1),-1,0,'Refined patterns parameters');
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
    patterns = GenerateReconstructionPatterns(params,PosRoiPatt(it,:),k(:,:,it),phase(:,:,it),params.pattAmp,sz+params.padSz/2,Lf);    

    % -- Reconstruction
    if params.nframes==1,  DispMsg(params.verbose,'<strong>===  Reconstruction</strong> ...'); 
    else, if params.sepOrr==0 && ~params.paraLoopFrames, DispMsg(params.verbose,'---> Reconstruct ...'); end; end
    rec(:,:,it) = Reconstruct(y(:,:,:,it),patterns,params);

    % - Save 
    if params.sav 
        if params.nframes>1
            saveastiff(rec(:,:,it),[params.DataPath(1:end-4),'_Rec_Frame_',num2str(it),'.tif']);
            saveastiff(gather(patterns),[params.DataPath(1:end-4),'_Patt_Frame_',num2str(it),'.tif']);
        else
            saveastiff(gather(patterns),[params.DataPath(1:end-4),'_Patt.tif']);
        end
    end
    if params.nframes>1
        if ~params.paraLoopFrames
            DispMsg(params.verbose,['== FlexSIM Elapsed time (s): ',num2str(toc(time0))]);
        else
            t = getCurrentTask();
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
        saveastiff(reshape(patt,[size(patt,1),size(patt,2),size(patt,3)*size(patt,4)]),[params.DataPath(1:end-4),'_Patt.tif']);
    end
end
DispMsg(params.verbose,['<strong>=== FlexSIM END. Elapsed time (s): ',num2str(toc(time0)),' </strong>']);


