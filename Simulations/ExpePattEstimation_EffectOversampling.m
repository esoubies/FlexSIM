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
% - Physical parameters
Na=1.2;                      % Numerical aperture
res=64;                      % Resolution
lamb=530;                    % Wavelength
nl=1.51;                     % RI oil
ns=1.333;                    % RI immersion medium (water)
nbPh=3;                      % Number of phases shifts considered in the simulation
pattAmp=0.9;                 % Pattern modulation amplitude
% - To define loops
MEP_test=[3 10 1000];  % Maximum expected number of photon per pixel
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
paramsEsti.ringMaskLim=[0.1,1.05];paramsEsti.eqPh=1;
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

%% Main Loop (over a_test and MEP_test)
% -- Initializations
err_k_fac2=zeros([Nrep,length(MEP_test)]);
err_k_fac4=zeros([Nrep,length(MEP_test)]);
err_k_fac8=zeros([Nrep,length(MEP_test)]);
err_k_fac16=zeros([Nrep,length(MEP_test)]);
err_ph_fac2=zeros([Nrep,length(MEP_test)]);
err_ph_fac4=zeros([Nrep,length(MEP_test)]);
err_ph_fac8=zeros([Nrep,length(MEP_test)]);
err_ph_fac16=zeros([Nrep,length(MEP_test)]);
err_frob_fac2=zeros([Nrep,length(MEP_test)]);
err_frob_fac4=zeros([Nrep,length(MEP_test)]);
err_frob_fac8=zeros([Nrep,length(MEP_test)]);
err_frob_fac16=zeros([Nrep,length(MEP_test)]);
patt = ones([sz 3]);
patt_est = ones([sz 3]);

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
        % - Common computations
        y = (y - min(y(:))) / (max(y(:)) - min(y(:)));   
        wf= (wf-min(wf(:))) / (max(wf(:)) - min(wf(:)));
        [G,wf] = RemoveWFandMask(y,wf,paramsEsti);
        FCut = 2*Na/lamb*res;        

        % - Oversampling factor 2
        tstart=tic;
        paramsEsti.overFac=2;
        [map,K1,K2] = CrossCorr(G,wf, paramsEsti);
        Jmap=-MaskFT(mean(abs(map).^2,3),FCut,paramsEsti.ringMaskLim);
        k_est= ExtractLocMin(paramsEsti,Jmap,K1,K2);
        id=intersect(find(k_est(1)==K1(:)),find(k_est(2)==K2(:)));
        [ii,jj]=ind2sub(size(K1),id);
        ph_est=mod(angle(map(ii,jj,:)),2*pi)/2;
        tt2=toc(tstart);
        for id=1:nbPh, patt_est(:,:,id) = gen_patt(k_est,ph_est+(id-1)*pi/3,pattAmp); end
        err_k_fac2(rep,idMEP) = norm((k_est-k).*sz(2:-1:1)*res/pi);
        err_ph_fac2(rep,idMEP) = min([abs(ph(1)-ph_est),abs(ph(1)-(ph_est-pi)),abs(ph(1)-(ph_est+pi))])/pi*180;
        err_frob_fac2(rep,idMEP) = norm(patt(:)-patt_est(:));

        % - Oversampling factor 4
        tstart=tic;
        paramsEsti.overFac=4;
        [map,K1,K2] = CrossCorr(G,wf, paramsEsti);
        Jmap=-MaskFT(mean(abs(map).^2,3),FCut,paramsEsti.ringMaskLim);
        k_est= ExtractLocMin(paramsEsti,Jmap,K1,K2);
        id=intersect(find(k_est(1)==K1(:)),find(k_est(2)==K2(:)));
        [ii,jj]=ind2sub(size(K1),id);
        ph_est=mod(angle(map(ii,jj,:)),2*pi)/2;
        tt4=toc(tstart);
        for id=1:nbPh, patt_est(:,:,id) = gen_patt(k_est,ph_est+(id-1)*pi/3,pattAmp); end
        err_k_fac4(rep,idMEP) = norm((k_est-k).*sz(2:-1:1)*res/pi);
        err_ph_fac4(rep,idMEP) = min([abs(ph(1)-ph_est),abs(ph(1)-(ph_est-pi)),abs(ph(1)-(ph_est+pi))])/pi*180;
        err_frob_fac4(rep,idMEP) = norm(patt(:)-patt_est(:));

        % - Oversampling factor 8
        tstart=tic;
        paramsEsti.overFac=8;
        [map,K1,K2] = CrossCorr(G,wf, paramsEsti);
        Jmap=-MaskFT(mean(abs(map).^2,3),FCut,paramsEsti.ringMaskLim);
        k_est= ExtractLocMin(paramsEsti,Jmap,K1,K2);
        id=intersect(find(k_est(1)==K1(:)),find(k_est(2)==K2(:)));
        [ii,jj]=ind2sub(size(K1),id);
        ph_est=mod(angle(map(ii,jj,:)),2*pi)/2;
        tt8=toc(tstart);
        for id=1:nbPh, patt_est(:,:,id) = gen_patt(k_est,ph_est+(id-1)*pi/3,pattAmp); end
        err_k_fac8(rep,idMEP) = norm((k_est-k).*sz(2:-1:1)*res/pi);
        err_ph_fac8(rep,idMEP) = min([abs(ph(1)-ph_est),abs(ph(1)-(ph_est-pi)),abs(ph(1)-(ph_est+pi))])/pi*180;
        err_frob_fac8(rep,idMEP) = norm(patt(:)-patt_est(:));

        % - Oversampling factor 16
        tstart=tic;
        paramsEsti.overFac=16;
        [map,K1,K2] = CrossCorr(G,wf, paramsEsti);
        Jmap=-MaskFT(mean(abs(map).^2,3),FCut,paramsEsti.ringMaskLim);
        k_est= ExtractLocMin(paramsEsti,Jmap,K1,K2);
        id=intersect(find(k_est(1)==K1(:)),find(k_est(2)==K2(:)));
        [ii,jj]=ind2sub(size(K1),id);
        ph_est=mod(angle(map(ii,jj,:)),2*pi)/2;
        tt16=toc(tstart);
        for id=1:nbPh, patt_est(:,:,id) = gen_patt(k_est,ph_est+(id-1)*pi/3,pattAmp); end
        err_k_fac16(rep,idMEP) = norm((k_est-k).*sz(2:-1:1)*res/pi);
        err_ph_fac16(rep,idMEP) = min([abs(ph(1)-ph_est),abs(ph(1)-(ph_est-pi)),abs(ph(1)-(ph_est+pi))])/pi*180;
        err_frob_fac16(rep,idMEP) = norm(patt(:)-patt_est(:));

        % === Displays
        clc;
        disp(['<strong>ExpePattEstimation: MEP = ',num2str(MEP_test(idMEP)),' | Rep #',num2str(rep),'/',num2str(Nrep),'</strong>']);
        disp(' Oversampling x2')
        disp_err('k  (px)',err_k_fac2,idMEP,rep);
        disp_err('ph  (°)',err_ph_fac2,idMEP,rep);
        disp_err('frob   ',err_frob_fac2,idMEP,rep);
        disp(['Time: ',num2str(tt2)])
        disp(' Oversampling x4')
        disp_err('k  (px)',err_k_fac4,idMEP,rep);
        disp_err('ph  (°)',err_ph_fac4,idMEP,rep);
        disp_err('frob   ',err_frob_fac4,idMEP,rep);
        disp(['Time: ',num2str(tt4)])
        disp(' Oversampling x8')
        disp_err('k  (px)',err_k_fac8,idMEP,rep);
        disp_err('ph  (°)',err_ph_fac8,idMEP,rep);
        disp_err('frob   ',err_frob_fac8,idMEP,rep);
        disp(['Time: ',num2str(tt8)])
        disp(' Oversampling x16')
        disp_err('k  (px)',err_k_fac16,idMEP,rep);
        disp_err('ph  (°)',err_ph_fac16,idMEP,rep);
        disp_err('frob   ',err_frob_fac16,idMEP,rep);
        disp(['Time: ',num2str(tt16)])
    end

    % === Displays 
    % - Box plots    
    cat = categorical(kron(ones(length(MEP_test),1),kron((1:6)',ones(Nrep,1))),1:6,{'tt1','Over sampling x2','Over sampling x4','Over sampling x8','Over sampling x16','tt'});
    mep_ = categorical(kron((1:length(MEP_test))',ones(6*Nrep,1)),1:length(MEP_test),cellfun(@(x) ['MEP = ',x],strsplit(num2str(MEP_test)),'UniformOutput',false));
    figure(5);
    k_err = [ones(Nrep,length(MEP_test))*NaN;err_k_fac2;err_k_fac4;err_k_fac8;err_k_fac16;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,k_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$\|k - \hat{k}\|$ [pixels]','Interpreter','latex'); ll=legend;ll.NumColumns=4;ylim([0 max(k_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    figure(6);
    ph_err = [ones(Nrep,length(MEP_test))*NaN;err_ph_fac2;err_ph_fac4;err_ph_fac8;err_ph_fac16;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,ph_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$|\phi - \hat{\phi}|$ (degree)','Interpreter','latex'); ll=legend;ll.NumColumns=4;ylim([0 max(ph_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    figure(7);
    frob_err = [ones(Nrep,length(MEP_test))*NaN;err_frob_fac2;err_frob_fac4;err_frob_fac8;err_frob_fac16;ones(Nrep,length(MEP_test))*NaN];
    bb=boxchart(mep_,frob_err(:),'GroupByColor',cat,'BoxWidth',0.75,'BoxFaceAlpha',1,'BoxEdgeColor','k'); legend;set(gca,'fontsize',14);grid;
    ylabel('$\|w - \hat{w}\|$','Interpreter','latex');ll=legend;ll.NumColumns=4;ylim([0 max(frob_err(:))]*1.2);
    bb(1).HandleVisibility='off';bb(end).HandleVisibility='off';
    drawnow;
end
figure(5);set(gca,'yscale','log')
figure(6);set(gca,'yscale','log')
figure(7);set(gca,'yscale','log')