%--------------------------------------------------------------------------
% FlexSIM simulation script on simulated data
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
imgType = 1;          
params.MEP = 2000; 
params.sz = [256 256];

params.DataPath = fullfile(pwd,sprintf('SIM_Simu_%d_%d.tif', imgType, params.MEP));    % Path to the SIM stack
params.pathToFlexSIM = '../';                      % Path to the root of GitHub FlexSIM repo

% -- Display, saving and GPU acceleration
params.displ = 2;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=1;                       % 0: minimal text displays / 1: detailled text displays
params.sav = 1;                         % Boolean if true save the result
params.GPU = 0;                         % Boolean on whether to use GPU or not

%% Data related parameters
% -- Properties of the SIM data stack
params.StackOrder= 'pa';               % Phase (p) and angle (a) convention. Choose one of ('paz', 'pza' or 'zap')
params.SzRoiBack=[];                   % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)
params.nbOr = 3;                       % Number of orientations
params.nbPh = 3;                       % Number of phases 

% -- OTF Approximation & Acquisition parameters
params.lamb = 530;                     % Emission wavelength
params.res = 64;                       % Pixel size (nm)
params.Na = 1.2;                       % Objective numerica aperture
params.damp = 0.8;                     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF
params.nl=1.51;                        % Refractive index of the objective medium (glass/oil)
params.ns=1.333;                       % Refractive index of the sample medium (water)
params.ph = [0 pi/3 2*pi/3] + rand*pi/3; % Phases used for acquisition simulation
params.or = [0 pi/3 2*pi/3] + rand*pi/3; % Orientationsused for acquisition simulation
params.a = 0.9;

%% FlexSIM parameters
% -- Patch-based processing
params.szPatch=0;                 % If >0, FlexSIM will perform pattern estimation and reconstruction by patches of size 'szPatch'
params.overlapPatch=0;            % Overlap between patches if szPatch>0
params.parallelProcess=0;         % If 1, paralellizes the loop over patches

% -- Parameters for patterns estimation
params.SzRoiPatt = [];            % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.limits = [0.8, 1.05];      % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0, 1.1];    % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.nMinima = 2;               % Number of starting points for the refinement steps
params.nPoints = 0;             % Number of points in the J evaluation grid. If set to 0, initialization is done via peak detection
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)
params.method = 2;                % Method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
                                  
% -- Parameters for image Reconstruction 
params.sepOrr = 1;                % Boolean if true treat each orientation separately
params.padSz=20;                  % Padding size for the optimization variable (to account for boundaries effects)
params.mu =  1e-3;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 1e-3;            % Relative error tolerance between two iterates (stopping criteria)

%% Generate data and run FlexSIM
GenerateSIMData(imgType, params);
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM
res = EvalRun(params, res);                            % Process the results