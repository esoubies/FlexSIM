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
[X,Y]=meshgrid(lh,lv);
[~,rho]=cart2pol(X,Y);

% Generate a normalized OTF
fc=2*Na/lamb*res;        % cut-off frequency
OTFatt=ifftshift((1-str*exp(-rho.^2/(2*(width*fc)^2))).*(double(rho<fc)));
end