function [corr,K1,K2] = CrossCorr(G,wf,params)
%--------------------------------------------------------------------------
% Function [corr,K1,K2] = CrossCorr(G,wf,params)
%
% Compute initial wave-vector through cross-correlation between the
% widefield image and the raw SIM images in Fourier domain.
% 
% Inputs : G       -> SIM images with the same orientation pattern  
%          wf      -> widefield image
%          params  -> Structures with fields:
%                         - lamb: Emission wavelength
%                         - Na: Objective numerica aperture
%                         - res: resolution of the SIM data stac
%         
% Outputs : corr   -> Correlation map
%           K1, K2 -> grid on which the correl map is evaluated
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
sz=size(wf);fac=4;
fftwf=fft2(padarray(wf,sz(1:2)*fac,'post'));
fftG=fft2(padarray(G,sz(1:2)*fac,'post'));
sz=size(fftwf);



% -- Perform cross-correlation btw the fft of wf and the fft of processed data G
% With padding for cross-corel
% fftwfpad=padarray(padarray(fftwf,ceil(sz(1:2)/2),'pre'),floor(sz(1:2)/2),'post');
% fftGpad=padarray(padarray(fftG,ceil(sz(1:2)/2),'pre'),floor(sz(1:2)/2),'post');
% corrtmp=fftshift(ifft2((fft2(ifftshift(fftwfpad))).*conj(fft2(ifftshift(fftGpad)))));
% corrtmp=corrtmp(ceil(sz(1)/2)+1:end-floor(sz(1)/2),ceil(sz(2)/2)+1:end-floor(sz(2)/2),:);

% Without padding for cross-corel
corrtmp=fftshift(ifft2((fft2(ifftshift(fftwf))).*conj(fft2(ifftshift(fftG)))));

if params.method==2
    wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
    tt=mean(corrtmp.*wght,3);tt2=conj(tt)./abs(tt);
    corr=MaskFT(real(mean(corrtmp.*wght,3).*tt2),FCut,params.limits);
    %corr=real(mean(corrtmp.*wght,3).*tt2);
    
    % Extract wave-vec and phase (for tests)
%     [~,id]=max(corr(:));[i,j]=ind2sub(size(corr),id);
%     kest=([j,i]-floor(size(corr)/2)-1)./size(corr).*size(wf);
%     disp(['Estimated wavec-vec [px w.r.t the ROI !]: [',num2str(kest(1)),',',num2str(kest(2)),']']);
%     ph=mod(-angle(tt2(id)),2*pi)/2;
%     disp(['Estimated phase offset (i.e. ph(0)): ',num2str(ph)])
else
    tt=conj(corrtmp)./abs(corrtmp); 
    corr=MaskFT(real(mean(corrtmp.*tt,3)),FCut,params.limits);
    %corr=real(mean(corrtmp.*tt,3));
    
    % Extract wave-vec and phase (for tests)
%     [~,id]=max(corr(:));  [i,j]=ind2sub(size(corr),id);
%     kest=([j,i]-floor(size(corr)/2)-1)./size(corr).*size(wf);
%     disp(['Estimated wavec-vec [px w.r.t the ROI !]: [',num2str(kest(1)),',',num2str(kest(2)),']']);
%     ph=mod(-angle(tt(i,j,:)),2*pi)/2;ph=ph(:)';
%     disp(['Estimated phase #1: ',num2str(ph(1))])
%     disp(['Estimated phase #2: ',num2str(ph(2))])
%     disp(['Estimated phase #3: ',num2str(ph(3))])
end


% Generate grid and get radius
if mod(sz(1),2)==0
    ll_v=linspace(-0.5,0,sz(1)/2+1);
    lr_v=linspace(0,0.5,sz(1)/2);
    lv=[ll_v,lr_v(2:end)];
else
    lv=linspace(-0.5,0.5,sz(1));
end
if mod(sz(2),2)==0
    ll_h=linspace(-0.5,0,sz(2)/2+1);
    lr_h=linspace(0,0.5,sz(2)/2);
    lh=[ll_h,lr_h(2:end)];
else
    lh=linspace(-0.5,0.5,sz(2));
end
[K1,K2]=meshgrid(lh*pi/params.res,lv*pi/params.res);

end