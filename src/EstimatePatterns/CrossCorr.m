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

% Pre-computations
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
fftwf=fftshift(fft2(wf));
fftG=fftshift(fft2(G));
sz=size(fftwf);

% Perform cross-correlation btw the fft of wf and the fft of processed data G
fftwfpad=padarray(padarray(fftwf,ceil(sz(1:2)/2),'pre'),floor(sz(1:2)/2),'post');
fftGpad=padarray(padarray(fftG,ceil(sz(1:2)/2),'pre'),floor(sz(1:2)/2),'post');
corrtmp=fftshift(ifft2(fft2(ifftshift(fftwfpad)).*fft2(ifftshift(fftGpad))));
corrtmp=corrtmp(ceil(sz(1)/2)+1:end-floor(sz(1)/2),ceil(sz(2)/2)+1:end-floor(sz(2)/2),:);
if params.method==2
    wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
    corr=MaskFT((abs(sum(corrtmp.*wght,3))),FCut,params.limits);
else
    corr=MaskFT((sum(abs(corrtmp),3)),FCut,params.limits);
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