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

% Declare a set of     
if ~isfield(params, "nMinima")
    params.nMinima = 1;
end
if ~isfield(params, "FilterRefinement")
    params.FilterRefinement = 1;
end
if ~isfield(params, "ringMaskLim")
    params.ringMaskLim = [0, 1];
end