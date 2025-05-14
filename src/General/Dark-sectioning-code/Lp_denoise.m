function [Lo] = Lp_denoise(image,params)

%% 基本参数
Nx = params.Nx;
Ny = params.Ny;
NA = params.NA;
emwavelength = params.emwavelength;
pixel_size = params.pixelsize;

%% 其他参数
res = 0.5 * emwavelength / NA/ params.factor;     % resolution
sigmaLP = Ny / (res / pixel_size(1));    % objective cut-off frequency ???

%% 滤波
lp = lpgauss(Nx,Ny,sigmaLP);

fft_image = fftshift(fft2(image));
Lo = real(ifft2(ifftshift(fft_image.*fftshift(lp))));

end

function [ out ] = hpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.
out=1-lpgauss(H,W,SIGMA);
end

function [ out ] = lpgauss(H,W,SIGMA)
%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.
H = double(H);
W = double(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
temp0 = -floor(W/2);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
temp = -(x.^2/(kcx^2)+y.^2/(kcy^2));
out = ifftshift(exp(temp));
% out = ifftshift(exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));
end