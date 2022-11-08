%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the SIM Microtubules data provided at
% https://codeocean.com/capsule/8064732/tree/v1
%
% Reference paper for the dataset:
% High-speed image reconstruction for optically sectioned, super-resolution structured illumination microscopy. 
% Advanced Photonics, 4(2). Wang, Z., et al. (2022). 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'SIM_Microtubules.tif');    % Path to the SIM stack
params.pathToFlexSIM = '../../';                           % Path to the root of GitHub FlexSIM repo
if ~exist(params.DataPath, 'file')
    websave(params.DataPath,'https://files.codeocean.com/files/verified/052eb33f-ffcb-4817-9ac0-7b1562ea02e1_v1.0/code/0-testdata/COS7_Microtubulin_520nm_NA1.49_Mag90x_Frame1.tif');
end

% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=1;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean if true save the result
params.GPU = 0;                         % Boolean on whether to use GPU or not

%% Data related parameters
% -- Background noise
params.SzRoiBack=51;        % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder= 'pa';    % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 
params.pattAmp=1;           % Amplitude of the patterns in [0,1]

% -- OTF Approximation
params.lamb = 520;     % Emission wavelength
params.res = 72.2;     % Pixel size (nm)
params.Na = 1.49;      % Objective numerica aperture
params.damp = 1;       % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = [];            % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.limits = [0.5, 0.9];       % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0, 1.1];    % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.nMinima = 1;               % Number of starting points for the refinement steps
params.nPoints = 0;               % Number of points in the J evaluation grid. If set to 0, initialization is done via cross-correlation in Fourier domain.
params.FilterRefinement = 1;      % Number of times that the filter is upgraded (gradient descent cycles)
params.method = 2;                % Method : 0 - treat all images independently
                                  %          1 - use all images with same orientation to estimate a unique wavevector
                                  %          2 - 1 + assume equally spaced phases
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
                                  
%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=128;               % If >0, FlexSIM will perform the reconstruction by patches of size 'szPatch'
params.overlapPatch=32;           % Overlap between patches if szPatch>0
params.parallelProcess=0;         % If 1, paralellizes the loop over patches

% -- OTF Attenuation
params.OTFAttStr=0.9999;          % Strenght of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.5;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.sepOrr = 0;                % Boolean if true treat each orientation separately
params.padSz=10;                  % Padding size for the optimization variable (to account for boundaries effects)
params.mu =  3e-4;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





