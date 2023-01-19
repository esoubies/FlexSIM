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

% Without padding for cross-corel
corrtmp=fftshift(ifft2((fft2(ifftshift(fftwf))).*conj(fft2(ifftshift(fftG)))));

if params.eqPh
    wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
    tt=mean(corrtmp.*wght,3);tt2=conj(tt)./abs(tt);
    corr=MaskFT(real(mean(corrtmp.*wght,3).*tt2),FCut,params.ringMaskLim);
else
    tt=conj(corrtmp)./abs(corrtmp); 
    corr=MaskFT(real(mean(corrtmp.*tt,3)),FCut,params.ringMaskLim);
end


% Generate grid and get radius
lv = ((1:sz(1)) - floor(sz(1)/2)-1)*pi/params.res/sz(1);
lh = ((1:sz(2)) - floor(sz(2)/2)-1)*pi/params.res/sz(2);
[K1, K2] = meshgrid(lh, lv);
end