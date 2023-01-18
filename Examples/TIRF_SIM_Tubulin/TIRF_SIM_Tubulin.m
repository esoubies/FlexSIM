%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the TIRF-SIM Tubulin data provided at
% https://github.com/fairSIM/test-datasets/blob/master/TIRF-SIM-Georgia.md
%
% Reference paper for the dataset:
% Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ. 
% Nature communications, 7(1), 1-6. Müller, M., Mönkemöller, V., Hennig, S., Hübner, W., & Huser, T. (2016).
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'TIRF_SIM_Tubulin.tif');   % Path to the SIM stack 
params.pathToFlexSIM = '../../';                             % Path to the root of GitHub FlexSIM repo
if ~exist(params.DataPath, 'file')
    websave(params.DataPath, 'https://github.com/fairSIM/test-datasets/releases/download/TIRF-SIM-Georgia/TIRF_Tubulin_525nm.tif');
end
% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=2;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean if true save the result
params.GPU = 0;                         % Boolean on whether to use GPU or not

%% Data related parameters
% -- Background estimation
params.SzRoiBack=31;        % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 
params.pattAmp=1;           % Amplitude of the patterns in [0,1]

% -- OTF Approximation
params.lamb = 525;        % Emission wavelength
params.res = 63;          % Pixel size (nm)
params.Na = 1.49;         % Objective numerica aperture
params.damp = 0.3;        % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = [];            % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.limits = [0.95, 1.05];     % Ring over which the J function is evaluated for initializing (fc = 1)
params.ringMaskLim = [0.4, 1.1];  % Lower and upper limit of mask to finish hiding WF component, givien as factor of fc
params.nMinima = 1;               % Number of starting points for the refinement steps
params.eqPh = 1;                  % If true, equally-spaced phases are assumed
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
                                  
%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=0;                 % If >0, FlexSIM will perform the reconstruction by patches of size 'szPatch'
params.overlapPatch=0;            % Overlap between patches if szPatch>0
params.parallelProcess=0;         % If 1, paralellizes the loop over patches

% -- OTF Attenuation
params.OTFAttStr=0.999;           % Strenght of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.3;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.sepOrr = 0;                % Boolean if true treat each orientation separately
params.padSz=50;                  % Padding size for the optimization variable (to account for boundaries effects)
params.mu =  5e-5;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





