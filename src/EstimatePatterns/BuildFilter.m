function [attFilt, Filt] = BuildFilter(k, sz, OTF, params, grids)
%--------------------------------------------------------------------------
% [att_filt, filt] = BuildFilter(k, sz, params, grids)
%
% Builds filters that will be applied to A (`filt`) and to the WF (`attFilt`)
%
% Inputs :  k       -> The wavevector to construct the filter around
%           sz      -> The size of the filters
%           OTF     -> Optical transfer function of the system
%           params  -> Structure containing the input parameters of the system
%           grids   -> Structure that stores intermediate and final grids
%
% Outputs : attFilt -> Filter for the WF image
%           Filt    -> Filter for A 
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also FlexSIM.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

kPix = k .* sz(2:-1:1)' * params.res / pi;      % Get the wavevector in pixel units
if isfield(params,'ringMaskLim')             % Get radius of the circles that will make up the filters
    rr = min(0.5,1-params.ringMaskLim(1))* norm(kPix);
else
    rr = 0.5 * norm(kPix);
end
I = grids.I; J = grids.J;
OTF0=double(OTF.*ifftshift((sqrt((I-kPix(1)-floor(sz(2)/2)-1).^2+(J-kPix(2)-floor(sz(1)/2)-1).^2)<rr)+(sqrt((I+kPix(1)-floor(sz(2)/2)-1).^2+(J+kPix(2)-floor(sz(1)/2)-1).^2)<rr))>0);
OTFshift=ifftshift(imtranslate(fftshift(OTF),kPix(:)')+imtranslate(fftshift(OTF),-kPix(:)'));
attFilt = OTFshift.*OTF0;
Filt = OTF0.*OTF; 
end