function OTF = GenerateOTF(Na,lamb,sz,res,damp,shift)
%--------------------------------------------------------------------------
% function OTF = GenerateOTF(Na,lamb,sz,res,damp,shift)
%
% Generates the Optical Transfer Function as a function of the acquisition 
% optical parameters. 
%
% Inputs : Na    -> Numerical Aperture
%          lamb  -> Illumination wavelength
%          sz    -> Size of the returned image
%          res   -> Pixel size
%          damp  -> Damping parameter (in [0,1] where 1 = no damping) to 
%                   attenuate the medium frequency response of the OTF
%          shift -> optional parameter to shift the center (0,0) of the OTF
%
% Outputs : OTF -> Array containing the Optical Transfer Function
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

if nargin <6
    shift=[0,0];
end

lv = ((1:sz(1)) - floor(sz(1)/2)-1)/sz(1);
lh = ((1:sz(1)) - floor(sz(1)/2)-1)/sz(1);
[X,Y]=meshgrid(lh,lv);
X=X-shift(1);Y=Y-shift(2);
[~,rho]=cart2pol(X,Y);

% Generate a normalized OTF
fc=2*Na/lamb*res;        % cut-off frequency
rhotrunc=rho.*(rho<fc);       % to avoid temporary complex values outside the OTF support (with the acos)
OTF=ifftshift(1/pi*(2*acos(abs(rhotrunc)/fc)-sin(2*acos(abs(rhotrunc)/fc))).*(rho<fc).*damp.^(abs(rho)/fc));
OTF=double(OTF)/max(max(OTF));    
end