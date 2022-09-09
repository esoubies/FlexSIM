function [attFilt, filt] = BuildFilter(k, sz, params)
%--------------------------------------------------------------------------
% [att_filt, filt] = BuildFilter(k, sz, params, displ)
%
% Builds filters that will be applied to A (`filt`) and to the WF (`attFilt`)
%
% Inputs :  k       → The wavevector to construct the filter around
%           sz      → The size of the filters
%           params  → Structure containing the parameters of the system
%
% Outputs : attFilt → Filter for the WF image
%           Filt    → Filter for A 
%--------------------------------------------------------------------------
kPix = k .* sz(1:2)' * params.res / pi; 
if isfield(params,'ringMaskLim')
    rr = min(0.5,1-params.ringMaskLim(1))* norm(kPix);
else
    rr = 0.5 * norm(kPix);
end
OTF = params.OTF; I = params.I; J = params.J;
OTF0=double(OTF.*ifftshift((sqrt((I-kPix(1)-floor(sz(1)/2)-1).^2+(J-kPix(2)-floor(sz(2)/2)-1).^2)<rr)+(sqrt((I+kPix(1)-floor(sz(1)/2)-1).^2+(J+kPix(2)-floor(sz(2)/2)-1).^2)<rr))>0);
OTFshift=ifftshift(imtranslate(fftshift(OTF),kPix(:)')+imtranslate(fftshift(OTF),-kPix(:)'));
attFilt = OTFshift.*OTF0;
filt = OTF0.*OTF; 
end