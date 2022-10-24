function [k] = CrossCorr(G,params)
%--------------------------------------------------------------------------
% Function k = CrossCorr(G,params)
%
% Takes as input a preconditioned set of SIM images with equally spaced 
% phases and estimates the wavector through the cross correlation method 
% wavevector.
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

% Solve linear system to extract components
[m, n] = meshgrid(0:2, 0:params.nbPh-1);  % Calculate M matrix
phi = 2*pi*n/params.nbPh;
M = exp(i*m.*phi); M = pinv(M);
fftG = fft2(G);                            % Get FFT
sz = size(G); n = sz(1) * sz(2);           % Calculate necessary dimensions
fftG = reshape(fftG, [n, sz(3)])';         % Solve matrix multiplication
fftC = M * fftG;
fftC = reshape(fftC', sz);
corr = conv2(fftshift(abs(fftC(:,:,1))), ...    % Cross correlate shifted components
    fftshift(abs(fftC(:,:,2))), 'same'); 
FCut = 2*params.Na/params.lamb*params.res; % Filter center out
corr = MaskFT(corr, FCut, [0.7, 1.1]); 

wv_x = linspace(-(sz(2)-1)/2,(sz(2)-1)/2,sz(2)); % X and Y indices of centered
wv_y = linspace(-(sz(1)-1)/2,(sz(1)-1)/2,sz(1)); % image

[~,idxMax]=max(corr(:));                     % Extract and store maxima in pxl coordinates
[ii,jj]=ind2sub(size(G(:,:,1)),idxMax); 
k = [wv_x(jj) wv_y(ii)]; 
sz = sz(1:2);
k = k*pi./sz/params.res;                     % Change units and apply positive x convention
k = k*sign(k(1));

if params.displ == 2
    figure; sgtitle('Wavevector Initialization through Cross-Correlation') 
    subplot(2, 2, 1); imshow(fftshift(abs(fftC(:,:,1))), []); colormap(viridis); title('$C_{0}(\mathbf{k})$', Interpreter='latex'); impixelinfo
    subplot(2, 2, 2); imshow(fftshift(abs(fftC(:,:,2))), []); colormap(viridis); title('$C_{-p_{\theta}}(\mathbf{k})$', Interpreter='latex'); impixelinfo
    subplot(2, 2, 3); imshow(fftshift(abs(fftC(:,:,3))), []); colormap(viridis); title('$C_{p_{\theta}}(\mathbf{k})$', Interpreter='latex'); impixelinfo
    subplot(2, 2, 4); imshow(corr, []); colormap(viridis); title('$C_{0}(\mathbf{k})\ast C_{p_{\theta}}(\mathbf{k})$', Interpreter='latex'); impixelinfo
    hold on; plot(jj, ii, 'ro', 'MarkerSize', 8, 'LineWidth', 3); legend('Detected $\mathbf{k}$', Interpreter='latex');
end
end