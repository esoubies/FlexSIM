function res = FlexSIM(params)
%--------------------------------------------------------------------------
% Function res = FlexSIM(params)
%
% FlexSIM encloses the pipeline of the end-to-end (raw images to final) SIM
% reconstruction. It should be called from a script defining the structure 
% `params`, which includes all the paths and optical parameters necessary
% for a reconstruction. 
%
% Inputs  : params → Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Outputs : res    → Structure with the final image (in the field `res`)
%                    and other intermediate results, like patterns.  
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Data loading + routinary checks
addpath(genpath(pwd))                         % Add subfolders to path
y = double(loadtiff(params.DataPath));        % Read data and add to params
params.y = y; 
params = Tests(params); 

%% Pattern Estimation
disp('=== Patterns parameter estimation START ...');
% Estimate pattern parameters and save in fileTxt
tic
[k, phase, a, results] = EstimatePatterns(params, params.displ);
toc
paramsarray = horzcat([k, phase, a]);
disp('=== Pattern parameter succesfully estimated and saved.');

[patterns, params] = GenerateReconstructionPatterns(params, results); 
saveastiff(patterns,pattname);
disp('=== Patterns generated and saved as TIFF.');   

%% Display of results
% Sample Patterns 
figure; imagesc(results.patterns); axis image; title('OTF'); 
%     viscircles(floor(sz(1:2)/2)+1,fc*sz(1)); 
colormap(viridis); colorbar
subplot(1,2,2);imagesc(fftshift(PSF)); axis image; title('PSF'); caxis([0 0.05*max(PSF(:))]); colorbar

% OTF
if displ
    figure;subplot(1,2,1);imagesc(log(1+abs(fftshift(OTF)))); axis image; title('OTF'); 
%     viscircles(floor(sz(1:2)/2)+1,fc*sz(1)); 
    colormap(viridis); colorbar
    subplot(1,2,2);imagesc(fftshift(PSF)); axis image; title('PSF'); caxis([0 0.05*max(PSF(:))]); colorbar
end

% All the wavevectors found
if displ    
    figure; imshow(log10(sum(abs(fftshift(fft2(y))),3)+1), []); colormap(viridis); 
    colorbar; axis on; hold on;
    % Plot cross at row 100, column 50
    for i = 1:params.nbOr
        tmp = k(i, :) * sz(1) * params.res / pi + sz(1)/2+1; 
        plot(tmp(1), tmp(2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    end
    title('Sum of data Fourier amplitudes + detected wavevectors'); drawnow;
end

end

%% Reconstruction