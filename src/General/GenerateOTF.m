function OTF = GenerateOTF(Na,lamb,sz,res,damp)
%--------------------------------------------------------------------------
% function OTF = GenerateOTF(Na,lamb,sz,res,damp)
%
% Generates the Optical Transfer Function as a function of the acquisition 
% optical parameters. 
%
% Inputs : Na   -> Numerical Aperture
%          lamb -> Illumination wavelength
%          sz   -> Size of the returned image
%          res  -> Pixel size
%          damp -> Damping parameter (in [0,1] where 1 = no damping) to 
%                  attenuate the medium frequency response of the OTF
%
% Outputs : OTF -> Array containing the Optical Transfer Function
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Generate grid and get radius
if mod(sz(1),2)==0
    ll_v=linspace(-0.5,0,sz(1)/2+1);
    lr_v=linspace(0,0.5,sz(1)/2);
    ll_h=linspace(-0.5,0,sz(2)/2+1);
    lr_h=linspace(0,0.5,sz(2)/2);
    [X,Y]=meshgrid([ll_h,lr_h(2:end)],[ll_v,lr_v(2:end)]);
else
    ll_v=linspace(-0.5,0.5,sz(1));
    ll_h=linspace(-0.5,0.5,sz(2));
    [X,Y]=meshgrid(ll_h,ll_v);
end
[~,rho]=cart2pol(X,Y);

% Generate a normalized OTF
fc=2*Na/lamb*res;        % cut-off frequency
OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc).*damp.^(abs(rho)/fc));
OTF=double(OTF)/max(max(OTF));    
end