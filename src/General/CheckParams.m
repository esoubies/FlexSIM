function CheckParams(params)
%--------------------------------------------------------------------------
% Function res = CheckParams(params)
%
% Tests is the first function called by FlexSIM, and does basic searches,
% variable intialization, etc. for the proper functioning of FlexSIM
%
% Inputs  : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies  (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

missg="Missing parameter ";
invld="Invalid parameter ";
%% General parameters
% -- Path and files
prefix="[Paths and files] ";
% DataPath
msg="Should be in tif format.";
assert(isfield(params, "DataPath"),prefix + missg + "`DataPath`. " + msg);
assert(strcmp(params.DataPath(end-3:end),'.tif'), prefix + invld + "`DataPath'. " + msg);  
% pathToFlexSIM
assert(isfield(params, "DataPath"),prefix + missg + "`pathToFlexSIM`. ");

% -- Display and saving
prefix="[Display and saving] ";
% dipsl and sav
msg="Should be a boolean.";
assert(isfield(params, "displ"),prefix + missg + "`displ`. " + msg);
assert(params.displ==0 || params.displ==1, prefix + invld + "`displ'. " + msg);  
assert(isfield(params, "sav"),prefix + missg + "`sav`. " + msg);
assert(params.sav==0 || params.sav==1, prefix + invld + "`sav'. " + msg);  


%% Data related parameters
% -- Properties of the SIM data stack
prefix="[SIM data Stack] ";
% nbOr and nbPh
msg="Should be a positive integer.";
assert(isfield(params, "nbOr"),prefix + missg + "`nbOr`. " + msg);
assert(mod(params.nbOr, 1) == 0 && params.nbOr> 0,prefix + invld + "`nbOr`. " + msg);
assert(isfield(params, "nbPh"),prefix + missg + "`nbPh`. " + msg);
assert(mod(params.nbPh, 1) == 0 && params.nbPh > 0,prefix + invld + "`nbPh`. " + msg);
% StackOrder
msg="Choose one of {`ap`, `pa`, `paw`, `apw`, `wap`, `wpa`} as acquisition convention (p(hase), a(ngle) and w(idefield)).";
assert(isfield(params, "StackOrder"),prefix + missg + "`StackOrder`. " + msg);
assert(ismember(params.StackOrder, ["ap", "pa", "paw", "apw", "wap", "wpa"]),prefix + invld + "`StackOrder`. " + msg);
if params.nbPh == 1       % Ensure that the widefield is there
    assert(any(ismember(char(params.StackOrder), 'w')),prefix + invld + "`StackOrder`. " + ...
        "Widefield image is necessary when providing one image per orientation. " + ...
        "Ensure that the acquisition convention is one of `paw, apw, wap, wpa`")
end
% SzRoiBack
msg="Should be either empty or an odd number.";
assert(isfield(params, "SzRoiBack"),prefix + missg + "`SzRoiBack`. " + msg);
assert(isempty(params.SzRoiBack) || mod(params.SzRoiBack,2)==1,prefix + invld + "`SzRoiBack`. " + msg);

% -- OTF Approximation
prefix="[OTF Approximation]";
% lamb, res, Na
msg="Should be a positive real.";
assert(isfield(params, "lamb"),prefix + missg + "`lamb`. " + msg);
assert(params.lamb> 0,prefix + invld + "`lamb`. " + msg);
assert(isfield(params, "res"),prefix + missg + "`res`. " + msg);
assert( params.res > 0,prefix + invld + "`res`. " + msg);
assert(isfield(params, "Na"),prefix + missg + "`Na`. " + msg);
assert( params.Na > 0,prefix + invld + "`Na`. " + msg);
% damp
msg="Should be a real in (0,1].";
assert(isfield(params, "damp"),prefix + missg + "`damp`. " + msg);
assert(params.damp> 0 && params.damp <=1 ,prefix + invld + "`damp`. " + msg);

%% FlexSIM parameters
% -- Patch-based processing
prefix="[Patch-based processing] ";
% szPatch
msg="Should be a non-negative integer.";
assert(isfield(params, "szPatch"),prefix + missg + "`szPatch`. " + msg);
assert(params.szPatch>= 0,prefix + invld + "`szPatch`. " + msg);
if ~isempty(params.SzRoiBack) && params.szPatch>0, assert(params.szPatch>=params.SzRoiBack,"Incompatible `szPatch` and `SzRoiBack` paramater. Should be set as szPatch > SzRoiBack"); end
if params.szPatch>0 % Check the two other parameters only if patch-based process is active (szPatch>0)
% overlapPatch
assert(isfield(params, "overlapPatch"),prefix + missg + "`overlapPatch`. " + msg);
assert(params.overlapPatch>= 0,prefix + invld + "`overlapPatch`. " + msg);
end

% -- Parameters for patterns estimation
prefix="[Patterns Estimation] ";
% SzRoiPatt
msg="Should be either empty or an odd number.";
assert(isfield(params, "SzRoiPatt"),prefix + missg + "`SzRoiPatt`. " + msg);
assert(isempty(params.SzRoiPatt) || mod(params.SzRoiPatt,2)==1,prefix + invld + "`SzRoiPatt`. " + msg);
if params.szPatch >0 && ~isempty(params.SzRoiPatt), assert(params.szPatch>=params.SzRoiPatt,"Incompatible `szPatch` and `SzRoiPatt` paramater. Should be set as szPatch > SzRoiPatt"); end
% limits
msg="Should be a vector of length 2 with values within [0,2].";
assert(isfield(params, "limits"),prefix + missg + "`limits`. " + msg);
assert(numel(params.limits) == 2,prefix + invld + "`limits`. " + msg);
assert(all(params.limits<= 2) && all(params.limits>=0) ,prefix + invld + "`limits`. " + msg);
% ringMaskLim
assert(isfield(params, "ringMaskLim"),prefix + missg + "`ringMaskLim`. " + msg);
assert(numel(params.ringMaskLim) == 2,prefix + invld + "`ringMaskLim`. " + msg);
assert(all(params.ringMaskLim<= 2) && all(params.ringMaskLim>=0) ,prefix + invld + "`ringMaskLim`. " + msg);
assert(max(params.limits(1)-0.5,0)>=params.ringMaskLim(1), prefix +"The WF masking `ringMaskLim(1)` is too large compared to the search region given in `limits`."+...
    "`ringMaskLim(1)` should not exceed " +num2str(max(params.limits(1)-0.5,0)));
% nMinima, FilterRefinement, nPoints
msg="Should be a positive integer.";
assert(isfield(params, "nMinima"),prefix + missg + "`nMinima`. " + msg);
assert(mod(params.nMinima, 1) == 0 && params.nMinima> 0,prefix + invld + "`nMinima`. " + msg);
assert(isfield(params, "FilterRefinement"),prefix + missg + "`FilterRefinement`. " + msg);
assert(mod(params.FilterRefinement, 1) == 0 && params.FilterRefinement> 0,prefix + invld + "`FilterRefinement`. " + msg);
assert(isfield(params, "nPoints"),prefix + missg + "`nPoints`. " + msg);
assert(mod(params.nPoints, 1) == 0 && params.nPoints> 0,prefix + invld + "`nPoints`. " + msg);
if params.nPoints < 50
    warning(prefix + "Parameter nPoints is low. Results might be inaccurate");
elseif params.nPoints > 300
    warning(prefix + "Parameter nPoints is large. Computation might be slow");
end
% method
msg="Should be an integer within {0, 1, 2}.";
assert(isfield(params, "method"),prefix + missg + "`method`. " + msg);
assert(mod(params.method, 1) == 0 && params.method > -1 && params.method < 3,prefix + invld + "`method`. " + msg);
if params.nbPh == 1              % Other check related to the choice of method
    assert(params.method == 0, "When providing only one phase, set parameter `method` should be set to 0.");
end
% estiPattLowFreq
msg="Should be a boolean.";
assert(isfield(params, "estiPattLowFreq"),prefix + missg + "`estiPattLowFreq`. " + msg);
assert(params.estiPattLowFreq==0 || params.estiPattLowFreq==1, prefix + invld + "`estiPattLowFreq'. " + msg);  

                                  
% -- Parameters for image Reconstruction 
prefix="[Image Reconstruction] ";
% mu, stepTol
msg="Should be a positive real.";
assert(isfield(params, "mu"),prefix + missg + "`mu`. " + msg);
assert(params.mu> 0,prefix + invld + "`mu`. " + msg);
assert(isfield(params, "stepTol"),prefix + missg + "`stepTol`. " + msg);
assert( params.stepTol > 0,prefix + invld + "`stepTol`. " + msg);
% padSZ, maxIt
msg="Should be a non-negative integer.";
assert(isfield(params, "padSz"),prefix + missg + "`padSz`. " + msg);
assert(params.padSz>= 0,prefix + invld + "`padSz`. " + msg);
assert(isfield(params, "maxIt"),prefix + missg + "`maxIt`. " + msg);
assert(params.maxIt>= 0,prefix + invld + "`maxIt`. " + msg);
% sepOrr
msg="Should be a boolean.";
assert(isfield(params, "sepOrr"),prefix + missg + "`sepOrr`. " + msg);
assert(params.sepOrr==0 || params.sepOrr==1, prefix + invld + "`sepOrr'. " + msg);  
% regType
msg="Should be an integer within {1, 2, 3}.";
assert(isfield(params, "regType"),prefix + missg + "`regType`. " + msg);
assert(mod(params.regType, 1) == 0 && params.regType > 0 && params.regType < 4,prefix + invld + "`regType`. " + msg);














