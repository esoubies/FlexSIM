%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the SIM Microtubules data provided at
% https://github.com/LanMai/OpenSIM/raw/master/SIMexpt/sim01z4.tif
%
% Reference paper for the dataset:
% Structured illumination microscopy image reconstruction algorithm. 
% IEEE Journal of Selected Topics in Quantum Electronics, 22(4), 50-63. Lal, A., Shan, C., & Xi, P. (2016). 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'5_Stand_Microtubules.tif');    % Path to the SIM stack
params.pathToFlexSIM = '../../';                                % Path to the root of GitHub FlexSIM repo
try
    if ~exist(params.DataPath, 'file')
        websave(params.DataPath, 'https://github.com/LanMai/OpenSIM/raw/master/SIMexpt/sim01z4.tif');
        y = mean(double(loadtiff(params.DataPath)),3);
        y=cat(3,reshape(y(1:512,:),[512,512,3]),cat(3,reshape(y(513:1024,:),[512,512,3]),reshape(y(1025:end,:),[512,512,3])));
        saveastiff(single(y(:,:,1:9)),params.DataPath); % Remove the 3 last slides of the stack that are not SIM data
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
% -- Background noise
params.SzRoiBack=151;       % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % SIM stack order ('ap','pa','apw','paw','wap','wpa') with a=angles, p=phases, w=widefield. For 'ap' and 'pa' a widefield per orientation is computed as the summing along phases. 
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 

% -- OTF Approximation
params.lamb = 515;     % Emission wavelength (nm)
params.res = 60;       % Pixel size (nm)
params.Na = 1.49;       % Objective numerica aperture
params.damp = 0.8;     % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = [];           % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.maskWF = 0.2;              % Radius (as a factor of the cutoff freq) of the disk used to mask central Fourier frequencies (attenuate residual WF contrib prior to cross-corr)
params.ringRegionSearch = [0 1];  % Lower and upper limits of Fourier ring region to search peaks (given as factor of the cutoff freq)
params.eqPh = 1;                  % Boolean, if true equally-spaced phases are assumed
params.nMinima = 1;               % Number of starting points for the refinement steps
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
params.pattAmp=0.8;                 % Amplitude of the patterns in [0,1]

%% Parameters for image Reconstruction 
% -- Patch-based processing
params.szPatch=0;               % Size (px) of patches (if 0, no patch-based processing)
params.overlapPatch=0;           % Overlap between patches

% -- OTF Attenuation
params.OTFAttStr=0.99;            % Strenght of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.1;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Operators and Costs
params.apodize = 1;               % Boolean on whether to use apodization on boundaries
params.sepOrr = 0;                % Boolean on whether to treat each orientation separately
params.padSz=10;                  % Padding size (px) used in the forward operator
params.mu =  8e-4;                % Regularization parameter
params.regType=1;                 % Choice regul: 1 for Tikhonov (i.e., Wiener), 2 for Total Variation, 3 for Good roughness

% -- Optim
params.maxIt = 100;               % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
res = FlexSIM(params);                                 % Run FlexSIM


