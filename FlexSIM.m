function res = FlexSIM(params)
%--------------------------------------------------------------------------
% Function res = FlexSIM(params)
% 
% FlexSIM [1] should be called from a script defining the structure `params`,
% which  includes all the paths and optical parameters necessary for a 
% reconstruction. (See folder Examples)
%
% Inputs : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Outputs: res    -> Structure with the final image (in the field `res`)
%                    and other intermediate results, like patterns.  
%
% [1] FlexSIM: ADD REF TO PAPER
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Data loading + routinary checks
y = double(loadtiff(params.DataPath));    % Read data and add to params
params.y=y;
Tests(params);                            % Check conformity of parameters

%% Pattern Estimation
disp('=== Patterns parameter estimation START ...');
% - Estimate pattern parameters
[k, phase, a, results] = EstimatePatterns(params, params.displ);
disp('=== Patterns parameter estimation END.');
% - Generate Patterns for reconstruction
results.a = results.a./results.a;a=a./a; % Hardcode to 1 for now (to be as in previous version)
[patterns, params] = GenerateReconstructionPatterns(params, results); 

% - Displays
imdisp(patterns(:,:,1),'Estimated Pattern #1',1);  % Default display
if params.displ ==1                                % Detailled display
    sz=size(y);
    fig=figure;imshow(log10(sum(abs(fftshift(fft2(y))),3)+1), []); colormap(viridis);
    colorbar; axis on; hold on; pause(0.1);
    for i = 1:params.nbOr
        tmp = k(i, :) * sz(1) * params.res / pi + sz(1)/2+1;
        figure(fig);plot(tmp(1), tmp(2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    end
    title('Sum of data Fourier amplitudes + detected wavevectors'); drawnow;
end

%% Reconstruction
rec = Reconstruct(y,patterns,params);

% - Displays
figure;subplot(1,2,1);imagesc(rec); axis image; title('Reconstructed image');colormap gray;
subplot(1,2,2);imagesc(log(1+abs(fftshift(fftn(rec))))); axis image; title('Reconstructed image FFT');

%% Save
% - Save reconstruction / patterns / parameters
if params.sav
    prefix=params.DataPath(1:end-4);
    saveastiff(single(rec),strcat(prefix,'_Rec.tif'));
    saveastiff(single(patterns),strcat(prefix,'_Patt.tif'));
    save(strcat(prefix,'_Params'),'params');
end

% - Fill output variable
res.rec=rec;
res.patt=patterns;

end

