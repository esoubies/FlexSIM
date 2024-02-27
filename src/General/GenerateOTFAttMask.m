function OTFatt = GenerateOTFAttMask(Na,lamb,sz,res,str,width)
%--------------------------------------------------------------------------
% function OTFatt = GenerateOTFAttMask(Na,lamb,sz,res,str,width)
%
% Generates the OTF attenuation mask proposed in [1].
%
% Inputs : Na    -> Numerical Aperture
%          lamb  -> Illumination wavelength
%          sz    -> Size of the returned image
%          res   -> Pixel size
%          str   -> strength of the attenuation (in [0,1])
%          width -> width of the attenuation (>0)
%
% Outputs : OTFatt -> Array containing the OTF attenuation mask
%
% [1] Phase optimisation for structured illumination microscopy. 
%     Optics express, 21(2), 2032-2049. (2013)
%     Wicker, K., Mandula, O., Best, G., Fiolka, R., & Heintzmann, R.
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% Generate grids
lv = ((1:sz(1)) - floor(sz(1)/2)-1)/sz(1);
lh = ((1:sz(2)) - floor(sz(2)/2)-1)/sz(2);
[X,Y]=meshgrid(lh,lv);
[~,rho]=cart2pol(X,Y);
X=gpuCpuConverter(X);
Y=gpuCpuConverter(Y);
rho=gpuCpuConverter(rho);

% Generate a normalized OTF
fc=2*Na/lamb*res;        % cut-off frequency
OTFatt=ifftshift((1-str*exp(-rho.^2/(2*(width*fc)^2))));
end
