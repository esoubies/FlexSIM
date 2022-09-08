%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FlexSIM.m - to be called from the file params.m and not modified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd))                   % Add subfolders to path
% run("src\General\Tests.m") / Tests
% run("src\General\Preprocess.m") / Preprocess
% disp('=== Patterns parameter estimation START ...');
% [k, phase, params] = EstimatePatterns(y, wfs, params, displ); 
% run("src\EstimatePatterns\EstimatePatterns.m")

% Pattern Estimation

% Tests(params)
tic
[y, wfs, params] = Preprocess(dataname, psfname, params, displ);

if EstiPatt
    disp('=== Patterns parameter estimation START ...');
    % Estimate pattern parameters and save in fileTxt
    [k, phase, a, params] = EstimatePatterns(y, wfs, params, displ);
    toc
    paramsarray = horzcat([k*params.sz(1)*params.res/pi, phase, a]);
    writematrix(paramsarray,fileTxt,'Delimiter',' ')
    disp('=== Pattern parameter succesfully estimated and saved.');
else
    % Read from fileTxt and
    T = table2array(readtable(fileTxt,'Delimiter',' '));
    params.k = T(:,1:2)/params.sz(1)/params.res*pi;
    params.phase = T(:,3:3+params.nbPh-1);
    params.a = T(:,3+params.nbPh:end);
end

if GenPatt
    [patterns, params] = GenerateReconstructionPatterns(params, 0); 
    saveastiff(patterns,pattname);
    disp('=== Patterns generated and saved as TIFF.');
end

%% Reconstruction