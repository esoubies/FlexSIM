%--------------------------------------------------------------------------
% FlexSIM reconstruction script of the TIRF-SIM Tubulin data provided at
% https://github.com/fairSIM/test-datasets/blob/master/TIRF-SIM-Georgia.md
%
% Reference papers for the dataset:
% Super-resolution video microscopy of live cells by structured illumination.
% Nature Methods, 6(5). P. Kner, B.B. Chhun, E.R. Griffis, L. Winoto, and M.G.L. Gustafsson, (2009).
%
% Open-source image reconstruction of super-resolution structured illumination microscopy data in ImageJ. 
% Nature communications, 7(1), 1-6. Müller, M., Mönkemöller, V., Hennig, S., Hübner, W., & Huser, T. (2016).
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
clear; close all; clc;

%% General parameters
% -- Path and files
params.DataPath = fullfile(pwd,'1_TIRF_Tubulin.tif');   % Path to the SIM stack 
params.pathToFlexSIM = '../../';                             % Path to the root of GitHub FlexSIM repo
try
    if ~exist(params.DataPath, 'file')
        websave(params.DataPath, 'https://github.com/fairSIM/test-datasets/releases/download/TIRF-SIM-Georgia/TIRF_Tubulin_525nm.tif');
    end
catch
    error('Automatic downloading and extraction of raw data failed [Script need to be adapted to your OS and installed packages]')
end

% -- Display, saving and GPU acceleration
params.displ = 1;                       % Displaying choice, from 0 to 2 with increasing number of display
params.verbose=2;                       % 0: minimal text displays / 1: detailled text displays / 2: more details
params.sav = 1;                         % Boolean on whether to save or not the reconstructed image and estimated patterns
params.parallelProcess=0;               % Boolean on whether to parallelize over orientation (if sepOrr=1) or over frames (temporal stack). Requires the parallel computing toolbox. 
params.GPU=0;                           % Boolean on whether to use GPU. Requires the parallel computing toolbox. 

%% Physical parameters and pre-processing
% -- Background estimation
params.SzRoiBack=31;        % Size (odd number or empty) of the ROI for background estimation (position automatically detected so as to minimize the intensity within the ROI)

% -- Patterns
params.StackOrder = 'pa';   % SIM stack order ('ap','pa','apw','paw','wap','wpa') with a=angles, p=phases, w=widefield. For 'ap' and 'pa' a widefield per orientation is computed as the summing along phases. 
params.nbOr = 3;            % Number of orientations
params.nbPh = 3;            % Number of phases 

% -- OTF Approximation
params.lamb = 525;        % Emission wavelength (nm)
params.res = 63;          % Pixel size (nm)
params.Na = 1.49;         % Objective numerical aperture
params.damp = 0.3;        % damping parameter (in [0,1], 1= no damping) to attenuate middle freq in the approx OTF

%% Parameters for patterns estimation
params.SzRoiPatt = [];            % Size (odd number or empty) of the ROI for pattern estimation (position automatically detected so as to maximize the intensity within the ROI)
params.maskWF = 0.4;              % Radius (as a factor of the cutoff freq) of the disk used to mask central Fourier frequencies (attenuate residual WF contrib prior to cross-corr)
params.ringRegionSearch = [0 1];  % Lower and upper limits of Fourier ring region to search peaks (given as factor of the cutoff freq)
params.eqPh = 1;                  % Boolean, if true equally-spaced phases are assumed
params.estiPattLowFreq=0;         % If true, estimate the low-freq component of the patterns
params.pattAmp=1;                 % Amplitude of the patterns in [0,1]
        
%% Parameters for image Reconstruction 
% -- OTF Attenuation
params.OTFAttStr=0.999;           % Strength of the OTF attenuation (in [0,1]). If 0 no OTF attenuation.
params.OTFAttwdth=0.3;            % Width of the OTF attenuation (>0). If 0 no OTF attenuation.

% -- Cost function
params.apodize = 0;               % Boolean on whether to use apodization on boundaries
params.sepOrr = 0;                % Boolean on whether to treat each orientation separately
params.padSz=20;                  % Padding size (px) used in the forward operator
params.mu =  1e-6;                % Regularization parameter
params.regType=1;                 % Regularization function: 1 for 1st-order Tikhonov, 2 for Total Variation, 3 for Good roughness.

% -- Optim
params.maxIt = 40;                % Maximum number of iterations (stopping criteria)
params.stepTol = 5e-4;            % Relative error tolerance between two iterates (stopping criteria)

%% Run FlexSIM
run(strcat(params.pathToFlexSIM, '/InstallFlexSIM.m')) % Take care of paths & GlobalBioIm
FlexSIM(params);                                 % Run FlexSIM





