function OTF = GenerateOTF(Na,lamb,sz,res,damp)
%--------------------------------------------------------------------------
% function OTF = GenerateOTF(Na,lamb,sz,res,damp,displ)
%
% Generates the Optical Transfer Function as a function of the acquisition 
% optical parameters. 
%
% Inputs : Na    : Numerical Aperture
%          lamb  : Illumination wavelength
%          sz    : Size of the returned image
%          res   : Pixel size
%          damp  : Damping parameter (in [0,1] where 1 = no damping) to 
%                  attenuate the medium frequency response of the OTF
%
% Outputs : OTF → Array containing the Optical Transfer Function
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Generate grid and get radius
if mod(sz(1),2)==0
    ll=linspace(-0.5,0,sz(1)/2+1);
    lr=linspace(0,0.5,sz(1)/2);
    [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
else
    ll=linspace(-0.5,0.5,sz(1));
    [X,Y]=meshgrid(ll,ll);
end
[~,rho]=cart2pol(X,Y);

% Generate a normalized OTF
fc=2*Na/lamb*res;        % cut-off frequency
OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc).*damp.^(abs(rho)/fc));
OTF=double(OTF)/max(max(OTF));    
end