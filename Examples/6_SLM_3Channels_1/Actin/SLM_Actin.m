%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the SIM Actine data provided at
% https://static-content.springer.com/esm/art%3A10.1186%2Fs43593-022-00035-x/MediaObjects/43593_2022_35_MOESM4_ESM.7z
%
% Reference paper for the dataset:
% Structured illumination microscopy based on principal component analysis.. 
% eLight, 3(1). Qian, J., et al. (2023). 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'6_SLM_Actin_1.tif');    % Path to the SIM stack
params.pathToFlexSIM = '../../../';                   % Path to the root of GitHub FlexSIM repo
try
    if ~exist(params.DataPath, 'file')
        websave([pwd,'/tmp'],'https://static-content.springer.com/esm/art%3A10.1186%2Fs43593-022-00035-x/MediaObjects/43593_2022_35_MOESM4_ESM.7z');
        if isunix || ismac,  ! 7zr x tmp.7z
        elseif ispc
        end
        y = double(loadtiff([pwd,'/Supplementary code and data/Test/RawImage/Ours/561_1.tif'])); 
        saveastiff(single(y),params.DataPath);
        rmdir('Supplementary code and data/','s');
        delete tmp.7z
    end
catch
    error('Automatic downloading and extraction of raw data failed [Script need to be adapted to your OS and installed packages]')
end

% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=2;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean on whether to save or not the reconstructed image and estimated patterns
params.GPU = 0;                         % Boolean on whether to use GPU or not

%% Data related parameters
% -- Background noise
params.SzRoiBack=151;       % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % SIM stack order ('ap','pa','apw','paw','wap','wpa') with a=angles, p=phases, w=widefield. For 'ap' and 'pa' a widefield per orientation is computed as the summing along phases. 
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 

% -- OTF Approximation
params.lamb = 607;     % Emission wavelength (nm)
params.res = 65;       % Pixel size (nm)
params.Na = 1.4;       % Objective numerica aperture
params.damp = 0.3;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = 511;           % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.ringMaskLim = [0.2,1];     % Lower and upper limits of Fourier ring mask (given as factor of the cutoff freq)
params.eqPh = 1;                  % Boolean, if true equally-spaced phases are assumed
params.nMinima = 1;               % Number of starting points for the refinement steps
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
params.pattAmp=0.6;                 % Amplitude of the patterns in [0,1]
                                  
%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=0;               % Size (px) of patches (if 0, no patch-based processing)
params.overlapPatch=0;           % Overlap between patches
params.parallelProcess=0;         % Boolean on whether to parallelize the loop over patches

% -- OTF Attenuation
params.OTFAttStr=0.99;            % Strenght of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.1;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.apodize = 0;               % Boolean on whether to use apodization on boundaries
params.sepOrr = 0;                % Boolean on whether to treat each orientation separately
params.padSz=10;                  % Padding size (px) used in the forward operator
params.mu =  5e-4;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM





