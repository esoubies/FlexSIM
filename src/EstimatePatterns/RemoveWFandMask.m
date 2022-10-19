function [G,wf] = RemoveWFandMask(y,wf,params)
%--------------------------------------------------------------------------
% function [G,wf] = RemoveWFandMask(y,wf,params)
%
% From a stack of SIM images emove the scaled WF + filter freq within a ring of interest
%
% Inputs  : y       -> SIM data stack
%           wf      -> widefield image
%           params  -> Structures with fields:
%                         - lamb: Emission wavelength
%                         - Na: Objective numerica aperture
%                         - res: resolution of the SIM data stac
%                         - ringMaskLim: 1x2 array (eg. [0.3, 1.1]) defining the ring used to mask the WF component and the high-freq of the data, givien as factor of fc=1
%           
% Outputs:  G       -> stack with WF removed and FT masked
%           wf      -> corresponding WF masked
%
% See also EstimatePatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Precomputations
sz=size(y);                                   % Size of the input stack
if length(sz) == 2; sz(3) = 1; end
[I,J]=meshgrid(1:sz(2),1:sz(1));              % Create a LP filter in Fourier domain to remove WF
FCut = 2*params.Na/params.lamb*params.res;    % Cut-off frequency
p = ([sz(2), sz(1)]+1)/2; r = sz(1)/8;               % Build mask
mask=double(((I-p(1)).^2+(J-p(2)).^2) < r^2); 

% -- Initializations
G=zeros(sz);
if params.GPU
    G = gpuArray(G);
end
fft_wf = fftshift(fft2(wf));

% -- Loop over images of the stack
for img = 1:sz(3)                       
    ffty=fftshift(fft2(y(:,:,img)));      % Calculate FFT
    a=real(OptWght(ffty,fft_wf, mask));   % Calculate argmin_a |a*ffty - fft_wf|^2
    % Remove the scaled WF + filter freq within a ring of interest
    G(:,:,img) = real(ifft2(ifftshift(MaskFT((a*ffty-fft_wf), FCut, params.ringMaskLim))));
end
% -- Mask widefield with the same ring as the images
wf = real(ifft2(ifftshift(MaskFT(fft_wf, FCut, params.ringMaskLim))));

end