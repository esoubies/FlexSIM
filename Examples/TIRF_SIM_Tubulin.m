%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the TIRF-SIM Tubulin data provided at
% https://github.com/fairSIM/test-datasets/blob/master/TIRF-SIM-Georgia.md
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% Parameters
% -- General
params.DataPath = './TIRF_SIM_Tubulin.tif';   % Path to the SIM stack 
params.displ = 1;                             % Displaying choice, from 0 to 2 with increasing number of display
params.sav = 1;                               % Boolean if true save the result

% -- Optical  and acquisition parameters
params.AcqConv= 'paz'; % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.lamb = 525;     % Emission wavelength
params.res = 63;       % Pixel size (nm)
params.Na = 1.49;      % Objective numerica aperture
params.damp = 0.3;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF
params.nbOr = 3;       % Number of orientations
params.nbPh = 3;       % Number of phases 

% -- Parameters for patterns estimation
params.roi = [];                  % Select ROI for pattern estimation ([initial y-coord, initial x-coord, size])
params.limits = [0.95, 1.05];     % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0.3, 1.1];  % Mask to finish hiding WF component
params.nMinima = 1;               % Number of starting points for the refinement steps
params.nPoints = 150;             % Number of points in the J evaluation grid
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)
params.method = 2;                % Method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases
                                  
% -- Parameters for image Reconstruction 
params.sepOrr = 0;         % Boolean if true treat each orientation separately
params.padSz=50;           % Padding size for the optimization variable (to account for boundaries effects)
params.mu =  5e-5;         % Regularization parameter
params.regType=1;          % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness
params.maxIt = 100;        % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;     % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
res = FlexSIM(params); 





