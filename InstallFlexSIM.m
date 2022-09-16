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
elseif ismac
    syst_type = 'M'; 
else
    syst_type = 'L'; 
end

if exist('GlobalBioIm', 'dir')             % Check existance for GlobalBioIm in path
    GBIPath = what('GlobalBioIm').path;
    addpath(genpath(GBIPath))
%     if exist('LinOp', 'dir')               % Check existance of specific functions to be used
%         if ~exist('LinOpConv', 'file')     % Check existance            
%             addpath(genpath('LinOp'))
%             savepath
%         end
%     else
%         addpath(genpath('GlobalBioIm'))
%         savepath
%     end    
else
    % Search on the whole computer?
%     disp('GlobalBioIm not found in PATH. Searching on computer ...');
%     switch syst_type
%         case 'W'
%             system()
%         case 'M'
%         case 'L'     
%     path_list = dir(fullfile('C:','**','GlobalBioIm'));  
%     if numel(path_list) == 0                             % If it is not found on the whole computer
    disp('GlobalBioIm not found in computer. Dowloading in current directory...');
    system('git clone https://github.com/Biomedical-Imaging-Group/GlobalBioIm') % Download GlobalBioIm in the same dir as  FlexSIM
    run('GlobalBioIm/setGlobalBioImPath.m')
%     else
%         run('.../setGlobalBioImPath.m')
    savepath
end
