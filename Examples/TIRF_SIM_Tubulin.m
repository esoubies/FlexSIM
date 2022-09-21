%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the TIRF-SIM Tubulin data provided at
% https://github.com/fairSIM/test-datasets/blob/master/TIRF-SIM-Georgia.md
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'TIRF_SIM_Tubulin.tif');   % Path to the SIM stack 
params.pathToFlexSIM = '../';                             % Path to the root of GitHub FlexSIM repo

% -- Display and saving
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.sav = 1;                         % Boolean if true save the result

%% Data related parameters
% -- Properties of the SIM data stack
params.AcqConv= 'paz';                  % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.nbOr = 3;                        % Number of orientations
params.nbPh = 3;                        % Number of phases 

% -- OTF Approximation
params.lamb = 525;     % Emission wavelength
params.res = 63;       % Pixel size (nm)
params.Na = 1.49;      % Objective numerica aperture
params.damp = 0.3;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% FlexSIM parameters
% -- General
params.szPatch=0;                 % If >0, FlexSIM will perform pattern estimation and reconstruction by patches of size 'szPatch'
params.overlapPatch=64;           % Overlap between patches if szPatch>0
params.enhanceContrast=0;         % When 0 (and patch-based) rescale the reconstruction to widefield intensity variations. When 1, keep enhanced contrast due to patch-based rec

% -- Parameters for patterns estimation
params.roi = [];                  % Select ROI for pattern estimation ([initial y-coord, initial x-coord, size])
params.limits = [0.95, 1.05];     % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0.3, 1.1];  % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.nMinima = 2;               % Number of starting points for the refinement steps
params.nPoints = 150;             % Number of points in the J evaluation grid
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)
params.method = 2;                % Method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
                                  
% -- Parameters for image Reconstruction 
params.sepOrr = 0;                % Boolean if true treat each orientation separately
params.padSz=50;                  % Padding size for the optimization variable (to account for boundaries effects)
params.mu =  5e-5;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





