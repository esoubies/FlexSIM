function [PSF,OTF] = GeneratePSF(Na,lamb,sz,res,type,damp,displ)
%--------------------------------------------------------------------------
% function [PSF,OTF] = GeneratePSF(Na,lamb,sz,res,type,damp,displ)
%
% Generates the PSF
%
% Inputs : Na    : Numerical Aperture
%          lamb  : Illumination wavelength
%          sz    : size of the returned image
%          res   : Pixel size
%          type  : 0 -> ideal sinc PSF (unit disk OTF)
%                  1 -> more realistic model
%          damp  : dampening parameter (in [0,1] with 1 = no dampening) to attenuate the medium frequency response of the OTF
%          displ : display (default 0)
%
% Outputs : PSF : image of the PSF
%           OTF : image of the OTF (Fourier of PSF)
%--------------------------------------------------------------------------

if nargin <7
    displ=0;
end
if nargin <6
    damp=1;
    displ=0;
end

fc=2*Na/lamb*res;        % cut-off frequency
if mod(sz(1),2)==0
    ll=linspace(-0.5,0,sz(1)/2+1);
    lr=linspace(0,0.5,sz(1)/2);
    [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
else
    ll=linspace(-0.5,0.5,sz(1));
    [X,Y]=meshgrid(ll,ll);
end
[~,rho]=cart2pol(X,Y);
if  type == 0
    OTF=fftshift((rho<fc));
else
    OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc).*damp.^(abs(rho)/fc));
end
OTF=double(OTF)/max(max(OTF));    
PSF=real(ifft2(OTF));

if displ
    figure;subplot(1,2,1);imagesc(log(1+abs(fftshift(OTF)))); axis image; title('OTF'); 
%     viscircles(floor(sz(1:2)/2)+1,fc*sz(1)); 
    colormap(viridis); colorbar
    subplot(1,2,2);imagesc(fftshift(PSF)); axis image; title('PSF'); caxis([0 0.05*max(PSF(:))]); colorbar
end
end

