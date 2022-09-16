function CheckParams(params)
%--------------------------------------------------------------------------
% Function res = Tests(params)
%
% Tests is the first function called by FlexSIM, and does basic searches,
% variable intialization, etc. for the proper functioning of FlexSIM
%
% Inputs  : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% Initial tests to ensure that the necessary parameters for FlexSIM to run are present

% Tests on general parameters
assert(contains(params.DataPath, 'tif'), 'The SIM data should be in tiff format!')                     

% Test optical and acquisition parameters
assert(ismember(params.AcqConv, ["paz", "pza" ,"zap"]), ...
    'Choose one of {`paz`, `pza` or `zap`} as acquisition convention (p(hase), a(ngle) and Z)')
assert(50 < params.lamb && params.lamb < 1000, ...
    'Acquisition wavelength [nm] expected for `params.lamb` (visible light range accepted)')
assert(1 < params.Na && params.Na < 3, ...
    'Numerical aperture accepted ranges from 1 to 3')
assert(0 < params.damp && params.damp < 1, ...
    'Damping parameter is expected in the range [0, 1]')
assert(isfield(params, "nbOr"), ...
    'Missing parameter `nbOr` (positive integer)')
assert(mod(params.nbOr, 1) == 0 && params.nbOr> 0, 'Parameter `nbOr` should be a positive integer')
assert(isfield(params, "nbPh"), ...
    'Missing parameter `nbPh` (positive integer)')
assert(mod(params.nbPh, 1) == 0 && params.nbPh > 0, 'Parameter `nbPh` should be a positive integer')


% Parameters for patterns estimation
assert(ismember(numel(params.roi), [0, 3]), ...
    'The region of interest (params.roi) should either be an empty array or have the form `[initial y-coord, initial x-coord, size]`')
assert(isfield(params, "nMinima"),'Missing parameter `nMinima` (positive integer)')
assert(numel(params.limits) == 2, ...
    '`params.limits` should have a lower and upper bound for the FT masking, thus should have 2 elements')
assert(all(params.limits < 2), ['The limits for the wave vector search (params.limits) should be a ' ...
    'multiple of the cutoff frequency. The range accepted is [0, 2]'])
assert(numel(params.ringMaskLim) == 2, ...
    '`params.ringMaskLim` should have a lower and upper bound for the FT masking, thus should have 2 elements')
assert(all(params.ringMaskLim < 2), ['The limits for the WF and image masking (params.ringMaskLim) should ' ...
    'be a multiple of the cutoff frequency. The range accepted is [0, 2]'])
assert(mod(params.nMinima, 1) == 0 && params.nMinima > 0, 'Parameter `nMinima` should be a positive integer')
if params.nPoints < 50
    warning('The grid for the J landscape evaluation is too rough. Results might be inaccurate')
elseif params.nPoints > 300
    warning('The grid for the J landscape evaluation is too fine. Computation might be slow')
end
assert(isfield(params, "FilterRefinement"), ...
    'Missing parameter `FilterRefinement` (positive integer)')
assert(mod(params.FilterRefinement, 1) == 0 && params.FilterRefinement> 0, 'Parameter `nMinima` should be a positive integer')
assert(mod(params.method, 1) == 0 && params.method> 0 && params.method < 3, ...
    'Parameter `method` should be in {0, 1, 2}')



