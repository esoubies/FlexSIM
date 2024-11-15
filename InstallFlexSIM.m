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
addpath(pwd);
addpath(genpath([pwd,'/src']));      % Add source files to path
global syst_type                    % Keep track of system type
if ispc
    syst_type = 'W'; 
elseif ismac
    syst_type = 'M'; 
else
    syst_type = 'L'; 
end

if exist('GlobalBioIm', 'dir')             % Check existance for GlobalBioIm in path
    GBIPath = getfield(what('GlobalBioIm'),'path');
    addpath(genpath(GBIPath))
else
    disp('GlobalBioIm not found in computer. Dowloading in current directory...');
    websave([pwd,'/GlobalBioIm'], 'https://github.com/Biomedical-Imaging-Group/GlobalBioIm/archive/refs/heads/master.zip');
    unzip('GlobalBioIm.zip')
    movefile('GlobalBioIm-master', 'GlobalBioIm')
    delete('GlobalBioIm.zip')
    run('GlobalBioIm/setGlobalBioImPath.m')
    savepath
end
