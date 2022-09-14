function params = Tests(params)
%--------------------------------------------------------------------------
% Function res = Tests(params)
%
% Tests is the first function called by FlexSIM, and does basic searches,
% variable intialization, etc. for the proper functioning of FlexSIM
%
% Inputs  : params → Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Outputs : params → Same structure as in the input, revised   
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% Initial tests to ensure that the necessary parameters for FlexSIM to run
% are present
assert(isfield(params, "nMinima"),'Missing parameter `nMinima` (positive integer)')
assert(isfield(params, "FilterRefinement"),'Missing parameter `FilterRefinement` (positive integer)')
assert(isfield(params, "ringMaskLim"),'Missing parameter `ringMaskLim` (array of two doubles )')
% Window vs Linux/Mac
% for n = ["nMinima", "FilterRefinement", "ringMaskLim" ]
%     assert(isfield(params, n),'Missing parameter '+n)
% end