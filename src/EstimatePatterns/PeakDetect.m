function [k] = PeakDetect(G,params)
%--------------------------------------------------------------------------
% Function k = PeakDetect(G,params)
%
% Takes as input a preconditioned set of SIM images, adds them and finds a
% detected wavevector from the location of the pixel with the highest 
% intensity within the right half of the image and within a radius of 
% [0.5 - 1.1] of the cutoff frequency. 
% 
% Inputs : G        -> SIM images with the same orientation pattern  
%           params  -> Structures with fields:
%                         - lamb: Emission wavelength
%                         - Na: Objective numerica aperture
%                         - res: resolution of the SIM data stac
% Outputs : k       -> Detected wavevector 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

FCut = 2*params.Na/params.lamb*params.res;    % Cut-off frequency
fftg = abs(fftshift(fft2(G))); 
fftg = sum(fftg, 3);
fftg = MaskFT(fftg, FCut, [0.5, 1.1]); 
sz = size(fftg); 
if mod(sz(2), 2) == 1                                    % Odd in x scenario
    fftg_r = fftg(:, ceil(sz(2)/2):end);                % from kx == 0 to max
    wv_x = linspace(0,floor(sz(2)/2),ceil(sz(2)/2));
else
    fftg_r = fftg(:, sz(2)/2:end);                      % from kx == 
    wv_x = linspace(-0.5,(sz(2)-1)/2, size(fftg_r, 2));    
end
wv_y = linspace(-(sz(1)-1)/2,(sz(1)-1)/2,sz(1));

[~,idxMax]=max(fftg_r(:));                     % Extract and store maxima
[ii,jj]=ind2sub(size(fftg_r),idxMax);
k = [wv_x(jj) wv_y(ii)]; 
k = k*pi./sz/params.res;
k = k*sign(k(1));
if params.displ == 2
    figure; sgtitle('Wavevector Initialization through Peak Detection') 
    imshow(abs(fftg), []); colormap(viridis); title('$\Sigma FT\{G\}$', Interpreter='latex'); impixelinfo        
    hold on; plot(jj+size(G,2)-length(wv_x) , ii, 'ro', 'MarkerSize', 8, 'LineWidth', 3); legend('Detected $\mathbf{k}$', Interpreter='latex')
end
end