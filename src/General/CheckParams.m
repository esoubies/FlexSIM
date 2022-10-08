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
if params.nbPh == 1
    assert(params.method == 0, ...  % Ensure that the widefield is there
         "When providing only one phase, set `params.method == 1`.")
    assert(any(ismember(char(params.StackOrder), 'w')), ...  % Ensure that the widefield is there
         "Widefield image is necessary when providing one image pero orientation. " + ...
         "Ensure that the acquisition convention is one of `paw, apw, wap, wpa`")
end

% Test optical and acquisition parameters
assert(ismember(params.StackOrder, ["ap", "pa", "paw", "apw", "wap", "wpa"]), ...
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
assert(isfield(params, "SzRoiBack"),'Missing parameter `SzRoiBack` (should be either empty or an odd number)');
assert(isempty(params.SzRoiBack) || mod(params.SzRoiBack,2)==1,'Parameter SzRoiBack should be either empty or an odd number');
assert(isfield(params, "SzRoiPatt"),'Missing parameter `SzRoiPatt` (should be either empty or an odd number)');
assert(isempty(params.SzRoiPatt) || mod(params.SzRoiPatt,2)==1,'Parameter SzRoiPatt should be either empty or an odd number');

assert(isfield(params, "nMinima"),'Missing parameter `nMinima` (positive integer)')
assert(numel(params.limits) == 2, ...
    '`params.limits` should have a lower and upper bound for the FT masking, thus should have 2 elements')
assert(all(params.limits < 2), ['The limits for the wave vector search (params.limits) should be a ' ...
    'multiple of the cutoff frequency. The range accepted is [0, 2]'])
assert(numel(params.ringMaskLim) == 2, ...
    '`params.ringMaskLim` should have a lower and upper bound for the FT masking, thus should have 2 elements')
assert(all(params.ringMaskLim < 2), ['The limits for the WF and image masking (params.ringMaskLim) should ' ...
    'be a multiple of the cutoff frequency. The range accepted is [0, 2]'])
assert(max(params.limits(1)-0.5,0)>=params.ringMaskLim(1), ['The WF masking `ringMaskLim(1)'' is too large compared to the search region given in `limits''. `ringMaskLim(1)''' ...
    'should be set to ',num2str(max(params.limits(1)-0.5,0))])
assert(mod(params.nMinima, 1) == 0 && params.nMinima > 0, 'Parameter `nMinima` should be a positive integer')
if params.nPoints < 50
    warning('The grid for the J landscape evaluation is too rough. Results might be inaccurate')
elseif params.nPoints > 300
    warning('The grid for the J landscape evaluation is too fine. Computation might be slow')
end
assert(isfield(params, "FilterRefinement"), ...
    'Missing parameter `FilterRefinement` (positive integer)')
assert(mod(params.FilterRefinement, 1) == 0 && params.FilterRefinement> 0, 'Parameter `nMinima` should be a positive integer')
assert(mod(params.method, 1) == 0 && params.method > -1 && params.method < 3, ...
    'Parameter `method` should be in {0, 1, 2}')



