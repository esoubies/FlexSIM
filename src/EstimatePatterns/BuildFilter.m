function [OTFshiftCrop, OTFCrop] = BuildFilter(k, sz, OTF, params, grids)
%--------------------------------------------------------------------------
% [OTFshiftCrop, OTFCrop] = BuildFilter(k, sz, OTF, params, grids)
%
% Builds filters that will be applied to A (OTFCrop) and to b (OTFshiftCrop)
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
kPix = k .* sz(2:-1:1)' * params.res / pi;    % Get the wavevector in pixel units
rr = max(0.5,min(params.maskWF+0.3,1))*norm(kPix);  % Get radius of the circles that will make up the filters (at least a band of 0.3 should remains)
I = grids.I+1; J = grids.J+1;
% Build the filter corresponding to the shifted and cropped OTF (for b)
OTF0=double(OTF.*ifftshift((sqrt((I-kPix(1)-floor(sz(2)/2)-1).^2+(J-kPix(2)-floor(sz(1)/2)-1).^2)<rr)+(sqrt((I+kPix(1)-floor(sz(2)/2)-1).^2+(J+kPix(2)-floor(sz(1)/2)-1).^2)<rr))>0);
OTFshift1= GenerateOTF(params.Na,params.lamb,sz,params.res,params.damp,kPix./ sz(2:-1:1)');
OTFshift2= GenerateOTF(params.Na,params.lamb,sz,params.res,params.damp,-kPix./ sz(2:-1:1)');
OTFshiftCrop=MaskFT(fftshift(OTFshift1), FCut, [params.maskWF,1.1],kPix)+MaskFT(fftshift(OTFshift2), FCut, [params.maskWF,1.1],-kPix);
OTFshiftCrop=ifftshift(OTFshiftCrop).*OTF0;
% Build the filter corresponding to the cropped OTF
OTFCrop = OTF.*double(OTFshiftCrop>0); 
end