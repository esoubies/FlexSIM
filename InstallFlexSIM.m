function InstallFlexSIM()
%--------------------------------------------------------------------------
% Function InstallFlexSIM()
%
% Ensures that all the necessray functions of FlexSIM and GlobalBioIm are
% in the MATLab path.
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
addpath(genpath(pwd))               % Add source files to path
global syst_type                    % Keep track of system type
if ispc
    syst_type = 'W'; 
    MinGW = 'MATLAB Support for MinGW-w64 C/C++ Compiler'; 
    addOns = matlab.addons.installedAddons;
    if ~any(ismember(addOns.Name, MinGW))   % Check the existence of the necessary compiler
        error('MinGW-w64 C compiler (required for reconstruction) not found. Install MAtlab MinGW add-on by downloading the .mlpkginstall file from https://ch.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler and dragging it to the commnand line.')
%         disp(['MinGW-w64 C compiler (required for reconstruction) not found. Installing MAtlab add-on from https://www.mingw-w64.org/'])
%         MinGWAddon = 'MinGW.mlpkginstall'; 
%         websave(MinGWAddon, 'https://ch.mathworks.com/login?uri=https%3A%2F%2Fch.mathworks.com%2Fmatlabcentral%2Ffileexchange%2F52848-matlab-support-for-mingw-w64-c-c-compiler&form_type=community')
%         matlab.addons.install(MinGWAddon, true)
    end
elseif ismac
    syst_type = 'M'; 
else
    syst_type = 'L'; 
end

if exist('GlobalBioIm', 'dir')             % Check existance for GlobalBioIm in path
    GBIPath = what('GlobalBioIm').path;
    addpath(genpath(GBIPath))
else
    disp('GlobalBioIm not found in computer. Dowloading in current directory...');
    system('git clone https://github.com/Biomedical-Imaging-Group/GlobalBioIm') % Download GlobalBioIm in the same dir as  FlexSIM
    run('GlobalBioIm/setGlobalBioImPath.m')
    savepath
end
