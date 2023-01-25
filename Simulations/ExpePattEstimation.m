%--------------------------------------------------------------------------
%
%
%--------------------------------------------------------------------------
clear; close all; clc; 
run('../InstallFlexSIM.m');
rng(1);

%% Parameters
% - General
displ=0;                     % Boolean, if true will display intermediate images
sav=0;                       % Boolean, if true save the results
% - Physical parameters
Na=1.2;                      % Numerical aperture
res=64;                      % Resolution
lamb=530;                    % Wavelength
nl=1.51;                     % RI oil
ns=1.333;                    % RI immersion medium (water)
nbPh=3;                      % Number of phases shifts considered in the simulation
pattAmp=0.9;                 % Pattern modulation amplitude
% - To define loops
MEP_test=[3 5 10 100 1000];  % Maximum expected number of photon per pixel
Nrep=50;                    % Number of repetitions to be averaged for each MEP 

%% Pre-computed quantities
x = double(loadtiff('phantom.tif')); x=x/max(x(:));  % Load phamtom sample
% MAKE AND AUTOMATIC DOWNLOAD FROM THE NET
sz = size(x);                                        % Size of the sample
OTF=GenerateOTF(Na,lamb,sz,res,1);                   % Generate the OTF
% Function to generate patterns
[X,Y]=meshgrid(0:sz(2)-1,0:sz(1)-1); X=X*res; Y=Y*res;
gen_patt = @(k,p,a) 1 + a*(cos(2*p)*cos(2*(k(1)*X+k(2)*Y)) - sin(2*p)*sin(2*(k(1)*X+k(2)*Y)));
% Parameters for pattern estimation
paramsEsti.lamb=lamb;paramsEsti.Na=Na;paramsEsti.res=res;paramsEsti.nbOr=1;paramsEsti.GPU=0;
paramsEsti.nbPh=nbPh;paramsEsti.nMinima=1;paramsEsti.displ=0;paramsEsti.verbose=0;paramsEsti.damp=1;
paramsEsti.ringMaskLim=[0.1,1.05];paramsEsti.overFac=2;
% Function used in displays
disp_err = @(name,tab,idMEP,rep) fprintf('   - Error %s: %0.2e  [Mean (Std) over rep: %0.2e (%0.2e)]\n',name,tab(rep,idMEP),mean(tab(1:rep,idMEP)),std(tab(1:rep,idMEP)));
% Colors for boxplots
firstcol=[0 0.4470 0.7410];
seconcol =[0.8500 0.3250 0.0980];
newcolors = min([firstcol/1.3;firstcol/1.3;firstcol*1.2;seconcol/1.3;seconcol*1.2],1);
% Displays
if displ
    figure(1);
    subplot(1,2,1);imdisp(x,'Input Image',0);
    subplot(1,2,2);imdisp(log(1+abs(fftshift(fft2(x)))),'FFT Input Image',0);
    figure(2);imdisp(fftshift(OTF),'OTF',0);colormap parula;
end
if sav, delete *.png ; end% To remove previous runs

%% Main Loop (over a_test and MEP_test)
% -- Initializations
err_k_init_EqPh=zeros([Nrep,length(MEP_test)]);
err_k_ref_EqPh=zeros([Nrep,length(MEP_test)]);
err_ph_init_EqPh=zeros([Nrep,length(MEP_test)]);
err_ph_ref_EqPh=zeros([Nrep,length(MEP_test)]);
err_frob_init_EqPh=zeros([Nrep,length(MEP_test)]);
err_frob_ref_EqPh=zeros([Nrep,length(MEP_test)]);
err_k_init_NoEqPh=zeros([Nrep,length(MEP_test)]);
err_k_ref_NoEqPh=zeros([Nrep,length(MEP_test)]);
err_ph_init_NoEqPh=zeros([Nrep,length(MEP_test)]);
err_ph_ref_NoEqPh=zeros([Nrep,length(MEP_test)]);
err_frob_init_NoEqPh=zeros([Nrep,length(MEP_test)]);
err_frob_ref_NoEqPh=zeros([Nrep,length(MEP_test)]);
patt = ones([sz 3]);
patt_init = ones([sz 3]);
patt_ref = ones([sz 3]);

% -- Random generation of orientation and phases
ph_array  = (0:nbPh-1)*pi/3 + rand(Nrep,1)*ones(1,nbPh)*pi/3;   % Phases
orr_array = rand(Nrep,1)*pi;                                    % orientations
% -- Loops
for idMEP=1:length(MEP_test)
    for rep=1:Nrep
        % === SIM data simulation
        ph=ph_array(rep,:);
        orr=orr_array(rep);
        k  = 2*pi*ns/lamb*[cos(orr), sin(orr)]*Na/nl;  % Modulation wavevector

        % - Patterns Generation
        for id=1:nbPh
            patt(:,:,id) = gen_patt(k,ph(id),pattAmp);
        end

        % - SIM data generation + widefield
        y_noNoise = real(ifft2(fft2(x.*patt).*OTF));y_noNoise=y_noNoise/max(y_noNoise(:));
        scaling = 1e12/MEP_test(idMEP);
        y = imnoise(y_noNoise/scaling,'poisson')*scaling; % Add poisson noise
        wf = mean(y,3);

        % - Display
        if displ
            figure(3);
            subplot(2,2,1); imdisp(patt(:,:,1),'Patt #1',0);
            subplot(2,2,2); imdisp(log(1+abs(fftshift(fft2(patt(:,:,1))))),'FFT Patt #1',0);
            subplot(2,2,3); imdisp(y(:,:,1),'SIM data #1',0);
            subplot(2,2,4); imdisp(log(1+abs(fftshift(fft2(y(:,:,1))))),'FFT SIM data #1',0);
        end

        % === Pattern Estimation
        % - Equally-spaced phases
        paramsEsti.eqPh=1;
        tstart=tic;
        [k_ref, ph_ref, k_init, ph_init] = EstimatePatterns(paramsEsti, [1,1], y, 0, wf); 
        tteq=toc(tstart);
        for id=1:nbPh
            patt_init(:,:,id) = gen_patt(k_init,ph_init+(id-1)*pi/3,pattAmp); 
            patt_ref(:,:,id) = gen_patt(k_ref,ph_ref+(id-1)*pi/3,pattAmp); 
        end
        err_k_init_EqPh(rep,idMEP) = norm((k_init-k).*sz(2:-1:1)*res/pi);
        err_ph_init_EqPh(rep,idMEP) = min([abs(ph(1)-ph_init),abs(ph(1)-(ph_init-pi)),abs(ph(1)-(ph_init+pi))])/pi*180;
        err_frob_init_EqPh(rep,idMEP) = norm(patt(:)-patt_init(:));
        err_k_ref_EqPh(rep,idMEP) = norm((k_ref-k).*sz(2:-1:1)*res/pi);
        err_ph_ref_EqPh(rep,idMEP) = min([abs(ph(1)-ph_ref),abs(ph(1)-(ph_ref-pi)),abs(ph(1)-(ph_ref+pi))])/pi*180;
        err_frob_ref_EqPh(rep,idMEP) = norm(patt(:)-patt_ref(:));

        % - No equally-spaced phases
        paramsEsti.eqPh=0;paramsEsti.useFilter=1;
        tstart=tic;
        [k_ref, ph_ref, k_init, ph_init] = EstimatePatterns(paramsEsti, [1,1], y, 0, wf);
        ttnoeq=toc(tstart);
        if sign(k_ref(1))~=sign(k(1))
            k_ref = -k_ref;k_init = -k_init;
            ph_ref=mod(-ph_ref,pi);ph_init=mod(-ph_init,pi);
        end
        for id=1:nbPh
            patt_init(:,:,id) = gen_patt(k_init,ph_init(id),pattAmp); 
            patt_ref(:,:,id) = gen_patt(k_ref,ph_ref(id),pattAmp); 
        end
        err_k_init_NoEqPh(rep,idMEP) = norm((k_init-k).*sz(2:-1:1)*res/pi);
        err_ph_init_NoEqPh(rep,idMEP) = mean(min([abs(ph-ph_init);abs(ph-(ph_init-pi));abs(ph-(ph_init+pi))],[],1))/pi*180;
        err_frob_init_NoEqPh(rep,idMEP) = norm(patt(:)-patt_init(:));
        err_k_ref_NoEqPh(rep,idMEP) = norm((k_ref-k).*sz(2:-1:1)*res/pi);
        err_ph_ref_NoEqPh(rep,idMEP) = mean(min([abs(ph-ph_ref);abs(ph-(ph_ref-pi));abs(ph-(ph_ref+pi))],[],1))/pi*180;
        err_frob_ref_NoEqPh(rep,idMEP) = norm(patt(:)-patt_ref(:));

        % === Displays
        clc;
        disp(['<strong>ExpePattEstimation: MEP = ',num2str(MEP_test(idMEP)),' | Rep #',num2str(rep),'/',num2str(Nrep),'</strong>']);
        disp(' Equally-spaced phases assumption')
        disp_err('k initial (px)',err_k_init_EqPh,idMEP,rep);
        disp_err('k refined (px)',err_k_ref_EqPh,idMEP,rep);
        disp_err('ph initial (°)',err_ph_init_EqPh,idMEP,rep);
        disp_err('ph refined (°)',err_ph_ref_EqPh,idMEP,rep);
        disp_err('frob initial  ',err_frob_init_EqPh,idMEP,rep);
        disp_err('frob refined  ',err_frob_ref_EqPh,idMEP,rep);
        disp(['Time: ',num2str(tteq)])
        disp(' No equally-spaced phases assumption')
        disp_err('k initial (px)',err_k_init_NoEqPh,idMEP,rep);
        disp_err('k refined (px)',err_k_ref_NoEqPh,idMEP,rep);
        disp_err('ph initial (°)',err_ph_init_NoEqPh,idMEP,rep);
        disp_err('ph refined (°)',err_ph_ref_NoEqPh,idMEP,rep);
        disp_err('frob initial  ',err_frob_init_NoEqPh,idMEP,rep);
        disp_err('frob refined  ',err_frob_ref_NoEqPh,idMEP,rep);
        disp(['Time: ',num2str(ttnoeq)])
    end

    % === Displays 
    % - Examples of noise levels
    figure(4);
    y_crop=y(160:260,360:460); 
    subplot(2,length(MEP_test),idMEP);imdisp(y_crop,['MEP = ',num2str(MEP_test(idMEP))],0);
    ffty=log(1+abs(fftshift(fft2(y)))); kpx=k.*sz(2:-1:1)*res/pi+floor(sz/2)+1;
    ffty_crop=ffty(kpx(2)-50:kpx(2)+49,kpx(1)-50:kpx(1)+49);
    subplot(2,length(MEP_test),length(MEP_test)+idMEP);imdisp(ffty_crop,'',0);

    % - Box plots    
    cat = categorical(kron(ones(length(MEP_test),1),kron((1:6)',ones(Nrep,1))),1:6,{'tt1','Initial (Eq-Ph)','Refined (Eq-Ph)','Initial (No Eq-Ph)','Refined (No Eq-Ph)','tt'});
    mep_ = categorical(kron((1:length(MEP_test))',ones(6*Nrep,1)),1:length(MEP_test),cellfun(@(x) ['MEP = ',x],strsplit(num2str(MEP_test)),'UniformOutput',false));
    figure(5);
    k_err = [ones(Nrep,length(MEP_test))*NaN;err_k_init_EqPh;err_k_ref_EqPh;err_k_init_NoEqPh;err_k_ref_NoEqPh;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,k_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$\|k - \hat{k}\|$ [pixels]','Interpreter','latex'); colororder(newcolors);ll=legend;ll.NumColumns=4;ylim([0 max(k_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    figure(6);
    ph_err = [ones(Nrep,length(MEP_test))*NaN;err_ph_init_EqPh;err_ph_ref_EqPh;err_ph_init_NoEqPh;err_ph_ref_NoEqPh;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,ph_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$|\phi - \hat{\phi}|$ (degree)','Interpreter','latex'); colororder(newcolors);ll=legend;ll.NumColumns=4;ylim([0 max(ph_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    figure(7);
    frob_err = [ones(Nrep,length(MEP_test))*NaN;err_frob_init_EqPh;err_frob_ref_EqPh;err_frob_init_NoEqPh;err_frob_ref_NoEqPh;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,frob_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$\|w - \hat{w}\|$','Interpreter','latex');colororder(newcolors);ll=legend;ll.NumColumns=4;ylim([0 max(frob_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    drawnow;

    % - Save examples of data
    if sav
    imwrite(uint8(y_noNoise(:,:,1)/max(y_noNoise(:))*255),'Noiseless_SIM_data.png');
    ffty_noNoise=log(1+abs(fftshift(fft2(y_noNoise(:,:,1)))));
    imwrite(uint8(ffty_noNoise/max(ffty_noNoise(:))*255),'FF_Noiseless_SIM_data.png');
    imwrite(uint8(y_crop/max(y_crop(:))*255),['y_crop_MEP_',num2str(MEP_test(idMEP)),'.png']);
    imwrite(uint8(ffty_crop/max(ffty_crop(:))*255),['ffty_crop_MEP_',num2str(MEP_test(idMEP)),'.png']);
    end
end
figure(5);set(gca,'yscale','log')
figure(6);set(gca,'yscale','log')
figure(7);set(gca,'yscale','log')

