%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the MIA SIM Membrane data provided at
% https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/AtheiSIM%2014-52_488.tif
%
% Reference paper for the dataset:
% ML-SIM: universal reconstruction of structured illumination microscopy images using transfer learning. 
% Biomedical optics express, 12(5), 2720-2733. Christensen, C. N., Ward, E. N., Lu, M., Lio, P., & Kaminski, C. F. (2021).
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'4_Stand_Membrane.tif');    % Path to the SIM stack
params.pathToFlexSIM = '../../';                                % Path to the root of GitHub FlexSIM repo
if ~exist(params.DataPath, 'file')
    websave(params.DataPath, 'https://github.com/charlesnchr/ML-SIM/raw/master/Test_data/AtheiSIM%2014-52_488.tif');
end

% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=2;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean on whether to save or not the reconstructed image and estimated patterns
params.GPU = 0;                         % Boolean on whether to use GPU or not
params.parallelProcess=0;               % Boolean on whether to parallelize over patches (if szPatch>0) or over orientation (if sepOrr=1). Requires the parallel computing toolbox. 

%% Data related parameters
% -- Background noise
params.SzRoiBack=51;        % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % SIM stack order ('ap','pa','apw','paw','wap','wpa') with a=angles, p=phases, w=widefield. For 'ap' and 'pa' a widefield per orientation is computed as the summing along phases. 
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 

% -- OTF Approximation
params.lamb = 530;     % Emission wavelength (nm)
params.res = 64;       % Pixel size (nm)
params.Na = 1.2;       % Objective numerica aperture
params.damp = 0.3;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = 257;           % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.ringMaskLim = [0.2, 1.1];  % Lower and upper limits of Fourier ring mask (given as factor of the cutoff freq)
params.eqPh = 1;                  % Boolean, if true equally-spaced phases are assumed
params.nMinima = 1;               % Number of starting points for the refinement steps
params.estiPattLowFreq=1;         % If true, estimate the low-freq component of the patterns params.estiPattLowFreq=1;         % If true, estimate the low-freq component of the patterns
params.pattAmp=0.8;               % Amplitude of the patterns in [0,1]
                   
%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=0;                 % Size (px) of patches (if 0, no patch-based processing)
params.overlapPatch=0;            % Overlap between patches

% -- OTF Attenuation
params.OTFAttStr=0.99;            % Strenght of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.3;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.apodize = 0;               % Boolean on whether to use apodization on boundaries
params.sepOrr = 0;                % Boolean on whether to treat each orientation separately
params.padSz=20;                  % Padding size (px) used in the forward operator
params.mu =  3e-4;                % Regularization parameter
params.regType=1;                 % Regularization function: 1 for 1st-order Tikhonov, 2 for Total Variation, 3 for Good roughness.

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 1e-3;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





