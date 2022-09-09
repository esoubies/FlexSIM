%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main.m - Set the parameters of your acquisition and reconstruction %%%%
% 
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
%% Basic (necessary) parameters
% -- File and directory paths
params.DataPath = '../FlexTestMacro/SIMStack.tif'; % datapath = '../FlexTestSimu/SIMStack.tif'; 
params.displ = 1;             % Select the level of intermediate results, 0 → final result, 
                       % 1 → + patterns; 2 → all relevant intermediate results

% -- Optical  and acquisition parameters
params.AcqConv= 'paz'; % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.lamb = 510;     % Emission wavelength
params.res = 64.5;     % Pixel size (nm)
params.Na = 1.49;      % Objective numerica aperture
params.damp = 0.3;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF
params.nbOr = 3;       % Number of orientations
params.nbPh = 3;       % and phases 

%% Advanced parameters
% -- Parameters for patterns estimation
params.roi = [170 150 116];       % Select ROI ([initial y-coord, initial x-coord, size])
params.limits = [0.95, 1.05];     % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0.4, 1.1];  % Mask to finish hiding WF component
params.nMinima = 3;               % Number of starting points for the refinement steps
params.nPoints = 200;             % Number of points in the J evaluation grid
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

res = FlexSIM(params); 