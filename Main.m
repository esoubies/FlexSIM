%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params.m - Set the parameters of your acquisition and reconstruction %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic (necessary) parameters
% -- File and directory paths
basedir='../FlexTestMacro'; % basedir='../FlexTestSimu';
dataname = strcat(basedir,'/SIMStack.tif');    % File name with raw SIM data
fileTxt = strcat(basedir,'/PattParams.txt');   % Name to save/read estimated parameters
pattname = strcat(basedir,'/Patterns.tif');    % Name to save/read estimated patterns
psfname = strcat(basedir,'/PSF.tif');            % Name to save/read  psf
gtname = [];                                     % No ground truth
outFolder = strcat(basedir,'/');               % Folder to save results

% -- Select what will be computed
GenPatt = 1;            % If 1, generates the images of the estimated patterns, otherwise read file containing previous estimation 
EstiPatt = 1;           % If 1, estimate the patterns parameters, otherwise read file containing previous estimation
SIMrec = 0;             % If 1, performs the SIM reconstruction
displ = 1;               % Select the level of intermediate results; 0 - final result; 1 + patterns; 2 -> all relevant intermediate results

% -- Optical  and acquisition parameters
params.AcqConv= 'paz';  % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.lamb = 510;      % Emission wavelength
params.res = 64.5;      % Pixel size (nm)
params.Na = 1.49;       % Objective numerica aperture
params.damp = 0.3;      % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF
params.nbOr = 3;        % Number of orientations
params.nbPh = 3;        % and phases 
% params.pattAmp = 1;     % patterns amplitude (in [0,1]) ---- TO DEPRECATE, should be calculated

%% Advanced parameters
% -- Parameters for patterns estimation
params.roi = [170 150 116];       % Select ROI ([initial y-coord, initial x-coord, size])
params.limits = [0.95, 1.05];     % Ring over which the J function is evaluated for initializing (fc = 1)
params.store_all_ph = 1;
params.ringMaskLim = [0.4, 1.1];  % Mask to finish hiding WF component
params.n_minima = 3;              % Number of starting points for the refinement steps
params.n_pts = 200;               % Number of points in the initialization grid
params.method = 2;                % method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases

% -- Parameters for SIM Reconstruction 
sav = 1;                 % Boolean if true save the result
sepOrr = 1;              % Boolean if true treat each orientation separately
nbOrr = params.nbOr; nbPh = params.nbPh;
% Data
valback = 0;              % Background value
downFact = [2 2];         % Downsampling factors
% Operator + Objective Function
lamb = 8e-4;             % Hyperparameter (can be an array to loop)
padSz = 50;              % Padding size for the optimization variable (to accound for boundaries effects)
Reg=1;                 % Choice regul: 1 for Tikhonov, 2 for TV, 3 for Good roughness
% SIM reconstruction
maxIt = 100;            % Max iterations
StepTol = 5e-4;          % Relative error tolerance on iterates

run('FlexSIM.m')