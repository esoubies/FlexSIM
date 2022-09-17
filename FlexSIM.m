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
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Data loading + routinary checks
y = double(loadtiff(fullfile(params.WorkingDir, params.DataPath)));    % Read data 
sz_y=size(y);
CheckParams(params);                                                   % Check conformity of parameters

figure;sliceViewer(y);title('SIM Raw data');  % Display data
drawnow;set(gcf,'Visible','on');

%% Pattern Estimation
disp('=== Patterns parameter estimation START ...');
% - Estimate pattern parameters
[k, phase, a] = EstimatePatterns(params, y);
disp('=== Patterns parameter estimation END.');
% - Generate Patterns for reconstruction
a=a./a; % Hardcode to 1 for now (to be as in previous version) 
[patterns] = GenerateReconstructionPatterns(params,k,phase,a,sz_y); 

% - Displays
figure;sliceViewer(patterns, "Colormap", viridis);        % Default display of patterns
title('Estimated Patterns'); colorbar; 
if params.displ > 0                                       % Detailled displays s
    figure; subplot(2, 1, 1);                             % Subplot to display wavevectors...
    sz=size(y);                                           % as image, and...
    imdisp(log10(sum(abs(fftshift(fft2(y))),3)+1),'Detected wavevectors',0); colormap(viridis);
    colorbar; hold on; 
    for i = 1:params.nbOr
        tmp = k(i, :) * sz(1) * params.res / pi + sz(1)/2+1;
        plot(tmp(1), tmp(2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    end
    tFig = subplot(2, 1, 2);                    % As table in the same subplot as the wv 
    patternParams = horzcat(k, phase, a);       % Initialize the data
    T = array2table(patternParams, 'variableNames', {'Kx', 'Ky', 'Phase [rad]', 'Amplitude'});
%     hUI=uitable(tFig,'Data',T{:,:}); 
    hUI=uitable('Data',T{:,:}, 'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',tFig.Position); % TODO: 'ColumnWidth', 'fit' throwing error
end
drawnow;

%% Reconstruction
rec = Reconstruct(y,patterns,params);

% - Displays
figure;
s1=subplot(1,2,1);imdisp(rec,'SIM Reconstruction',0);
s2=subplot(1,2,2);imdisp(imresize(mean(y,3),size(rec)),'Widefield',0);
Link = linkprop([s1, s2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

%% Save
% - Save reconstruction / patterns / reconst. parameters / pattern parameters
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

