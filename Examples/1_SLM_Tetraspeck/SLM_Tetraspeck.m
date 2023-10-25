%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the SLM-SIM Tetraspeck data provided at
% https://github.com/fairSIM/test-datasets/blob/master/SLM-SIM-Bielefeld.md
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
params.DataPath = fullfile(pwd,'1_SLM_Tetraspeck.tif');    % Path to the SIM stack
params.pathToFlexSIM = '../../';                                % Path to the root of GitHub FlexSIM repo
try
    if ~exist(params.DataPath, 'file')
        websave(params.DataPath, 'https://github.com/fairSIM/test-datasets/releases/download/SLM-SIM-Bielefeld/SLM-SIM_Tetraspeck200_680nm.tif');
    end
catch
    error('Automatic downloading and extraction of raw data failed [Script need to be adapted to your OS and installed packages]')
end

% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=2;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean on whether to save or not the reconstructed image and estimated patterns
params.GPU = 0;                         % Boolean on whether to use GPU or not
params.parallelProcess=0;               % Boolean on whether to parallelize over patches (if szPatch>0) or over orientation (if sepOrr=1). Requires the parallel computing toolbox. 

%% Data related parameters
% -- Background estimation
params.SzRoiBack=13;        % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % SIM stack order ('ap','pa','apw','paw','wap','wpa') with a=angles, p=phases, w=widefield. For 'ap' and 'pa' a widefield per orientation is computed as the summing along phases. 
params.nbOr = 4;            % Number of orientations
params.nbPh = 3;            % Number of phases 

% -- OTF Approximation
params.lamb = 680;     % Emission wavelength (nm)
params.res = 97.6;     % Pixel size (nm)
params.Na = 1.2;       % Objective numerica aperture
params.damp = 0.4;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = 257;           % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.maskWF = 0;                % Radius (as a factor of the cutoff freq) of the disk used to mask central Fourier frequencies (attenuate residual WF contrib prior to cross-corr)
params.ringRegionSearch = [0 1];  % Lower and upper limits of Fourier ring region to search peaks (given as factor of the cutoff freq)
params.eqPh = 1;                  % Boolean, if true equally-spaced phases are assumed
params.nMinima = 1;               % Number of initial wavevectors extracted from cross-correlation maps
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
params.pattAmp=1;                 % Amplitude of the patterns in [0,1]
                       
%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=0;                 % Size (px) of patches (if 0, no patch-based processing)
params.overlapPatch=0;            % Overlap between patches

% -- OTF Attenuation
params.OTFAttStr=0;               % Strength of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0;              % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.apodize = 0;               % Boolean on whether to use apodization on boundaries
params.sepOrr = 0;                % Boolean on whether to treat each orientation separately
params.padSz=20;                  % Padding size (px) used in the forward operator
params.mu =  5e-5;                % Regularization parameter
params.regType=1;                 % Regularization function: 1 for 1st-order Tikhonov, 2 for Total Variation, 3 for Good roughness.

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 1e-3;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





