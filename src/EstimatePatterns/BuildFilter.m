function [OTFshiftCrop, OTFCrop] = BuildFilter(k, sz, OTF, params, grids)
%--------------------------------------------------------------------------
% [OTFshiftCrop, OTFCrop] = BuildFilter(k, sz, OTF, params, grids)
%
% Builds filters that will be applied to A (`filt`) and to the WF (`attFilt`)
%
% Inputs :  k       -> The wavevector to construct the filter around
%           sz      -> The size of the filters
%           OTF     -> Optical transfer function of the system
%           params  -> Structure containing the input parameters of the system
%           grids   -> Structure that stores intermediate and final grids
%
% Outputs : OTFshiftCrop -> Filter for b (OTF shifted and then cropped around peaks)
%           OTFCrop      -> Filter for A (OTF croped around peaks)
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also FlexSIM.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

FCut = 2*params.Na/params.lamb*params.res;    % Cut-off frequency
kPix = k .* sz(2:-1:1)' * params.res / pi;      % Get the wavevector in pixel units
if isfield(params,'ringMaskLim')             % Get radius of the circles that will make up the filters
    rr = min(0.5,1-params.ringMaskLim(1))* norm(kPix);
else
    rr = 0.5 * norm(kPix);
end
I = grids.I+1; J = grids.J+1;
% Masked OTF (as done for the WF in RemoveWFandMask)
MaskedOTF=MaskFT(fftshift(OTF), FCut, params.ringMaskLim);
% Build the filter corresponding to the shifted and cropped OTF (for b)
OTF0=double(OTF.*ifftshift((sqrt((I-kPix(1)-floor(sz(2)/2)-1).^2+(J-kPix(2)-floor(sz(1)/2)-1).^2)<rr)+(sqrt((I+kPix(1)-floor(sz(2)/2)-1).^2+(J+kPix(2)-floor(sz(1)/2)-1).^2)<rr))>0);
OTFshiftCrop=ifftshift(imtranslate(gather(MaskedOTF),kPix(:)')+imtranslate(gather(MaskedOTF),-kPix(:)')).*OTF0;
% Build the filter corresponding to the cropped OTF
OTFCrop = OTF.*double(OTFshiftCrop>0); 
end