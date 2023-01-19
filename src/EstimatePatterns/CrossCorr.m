function [Jmap,K1,K2] = CrossCorr(G,wf,params)
%--------------------------------------------------------------------------
% Function [Jmap,K1,K2] = CrossCorr(G,wf,params)
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
% Outputs : Jmap   -> Grid-based evaluation of J in [1]
%           K1, K2 -> grid on which the correl map is evaluated
%
% [1] FlexSIM: ADD REF TO PAPER
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com), 
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
fftshift_ = @(x) fftshift(fftshift(x,1),2);
ifftshift_ = @(x) ifftshift(ifftshift(x,1),2);
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
sz=size(wf);fac=4;
fftwf=fft2(padarray(wf,sz(1:2)*fac,'post'));
fftG=fft2(padarray(G,sz(1:2)*fac,'post'));
sz=size(fftwf);

% -- Perform cross-correlation btw the fft of wf and the fft of processed data G
% tmp=fft2(ifftshift(fftwf)).*conj(fft2(ifftshift_(fftG)));
% corr=fftshift_(ifft2(tmp));
% corrSymConj=fftshift_(ifft2(conj(tmp)));
% if params.eqPh
%     wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
%     tmp=mean((corr + corrSymConj).*wght,3);
%     Jmap=MaskFT(abs(tmp).^2,FCut,params.ringMaskLim);
% else
%     tmp=corr + corrSymConj;
%     Jmap=MaskFT(mean(abs(tmp).^2,3),FCut,params.ringMaskLim);
% end

corrtmp=fftshift(ifft2((fft2(ifftshift(fftwf))).*conj(fft2(ifftshift(fftG)))));
if params.eqPh
    wght=reshape(exp(-2*1i*[0:params.nbPh-1]*pi/params.nbPh),[1,1,params.nbPh]);
    tt=mean(corrtmp.*wght,3);tt2=conj(tt)./abs(tt);
    Jmap=MaskFT(real(mean(corrtmp.*wght,3).*tt2),FCut,params.ringMaskLim);
else
    tt=conj(corrtmp)./abs(corrtmp); 
    Jmap=MaskFT(real(mean(corrtmp.*tt,3)),FCut,params.ringMaskLim);
end

% Generate grid and get radius
lv = ((1:sz(1)) - floor(sz(1)/2)-1)*pi/params.res/sz(1);
lh = ((1:sz(2)) - floor(sz(2)/2)-1)*pi/params.res/sz(2);
[K1, K2] = meshgrid(lh, lv);
end