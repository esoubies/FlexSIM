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


% -- Display, saving and GPU acceleration
prefix="[Display, saving] ";
% dipsl
msg="Should be an integer within {0, 1, 2, 3}.";
assert(isfield(params, "displ"),prefix + missg + "`displ`. " + msg);
assert(mod(params.displ, 1) == 0 && params.displ > -1 && params.displ < 4,prefix + invld + "`displ`. " + msg);
% verbose
msg="Should be an integer within {0, 1, 2}";
assert(isfield(params, "verbose"),prefix + missg + "`verbose`. " + msg);
assert(params.verbose==0 || params.verbose==1 || params.verbose==2, prefix + invld + "`verbose'. " + msg);  
% sav
msg="Should be a boolean.";
assert(isfield(params, "sav"),prefix + missg + "`sav`. " + msg);
assert(params.sav==0 || params.sav==1, prefix + invld + "`sav'. " + msg);  
% GPU usage
%prefix="[GPU]";
%msg="Should be a boolean";
%assert(isfield(params, "GPU"),prefix + missg + "`GPU`. " + msg);
%assert(ismember(params.GPU, [0, 1]), prefix + invld + "`GPU`. " + msg);
%if params.GPU
%    assert(logical(license('test','Distrib_Computing_Toolbox')),prefix  + "Parallel Computing Toolbox not installed. Set parameter GPU to 0.")
%    assert(logical(gpuDeviceCount), prefix + "No GPU was found. Set parameter GPU to `0`.");
%end
% parallelProcess
prefix="[Parallel]";
msg="Should be a boolean.";
assert(isfield(params, "parallelProcess"),prefix + missg + "`parallelProcess`. " + msg);
assert(ismember(params.parallelProcess, [0, 1]), prefix + invld + "`parallelProcess'. " + msg);  
if params.parallelProcess % Check is parallel tolbox is available 
  assert(logical(license('test','Distrib_Computing_Toolbox')),prefix + "Parallel Computing Toolbox not installed. Set parameter parallelProcess to 0.")
  %assert(params.verbose==0,"prefix + parameter verbose should be set to 0 when parallelProcess is activated.");
end


%% Data related parameters
prefix="[Background estimation] ";
% -- Background estimation
% SzRoiBack
msg="Should be either empty or an odd number.";
assert(isfield(params, "SzRoiBack"),prefix + missg + "`SzRoiBack`. " + msg);
assert(isempty(params.SzRoiBack) || mod(params.SzRoiBack,2)==1,prefix + invld + "`SzRoiBack`. " + msg);

% -- Patterns
prefix="[Patterns] ";
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
% pattAmp
msg="Should be a real in (0,1].";
assert(isfield(params, "pattAmp"),prefix + missg + "`pattAmp`. " + msg);
assert(params.pattAmp> 0 && params.pattAmp <=1 ,prefix + invld + "`pattAmp`. " + msg);


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

%% Parameters for patterns estimation
prefix="[Patterns Estimation] ";
% SzRoiPatt
msg="Should be either empty or an odd number.";
assert(isfield(params, "SzRoiPatt"),prefix + missg + "`SzRoiPatt`. " + msg);
assert(isempty(params.SzRoiPatt) || mod(params.SzRoiPatt,2)==1,prefix + invld + "`SzRoiPatt`. " + msg);
% maskWF
msg="Should be a real in [0, 0.5].";
assert(isfield(params, "maskWF"),prefix + missg + "`maskWF`. " + msg);
assert(params.maskWF>= 0 && params.maskWF <=0.5 ,prefix + invld + "`maskWF`. " + msg);
% ringRegionSearch
msg="Should be a nonnegative 2D vector.";
assert(isfield(params, "ringRegionSearch"),prefix + missg + "`ringRegionSearch`. " + msg);
assert(numel(params.ringRegionSearch) == 2,prefix + invld + "`ringRegionSearch`. " + msg);
assert(all(params.ringRegionSearch>=0) ,prefix + invld + "`ringRegionSearch`. " + msg);
% eqPh
msg="Should be a boolean.";
assert(isfield(params, "eqPh"),prefix + missg + "`eqPh`. " + msg);
assert(params.eqPh==0 || params.eqPh==1, prefix + invld + "`eqPh`. " + msg);  
% estiPattLowFreq
msg="Should be a boolean.";
assert(isfield(params, "estiPattLowFreq"),prefix + missg + "`estiPattLowFreq`. " + msg);
assert(params.estiPattLowFreq==0 || params.estiPattLowFreq==1, prefix + invld + "`estiPattLowFreq'. " + msg);  

                                  
%% Parameters for image Reconstruction 
% -- OTF Attenuation
prefix="[OTF Attenuation] ";
% OTFAttStr
msg="Should be a real in [0,1].";
assert(isfield(params, "OTFAttStr"),prefix + missg + "`OTFAttStr`. " + msg);
assert(params.OTFAttStr>= 0 && params.OTFAttStr <=1 ,prefix + invld + "`OTFAttStr`. " + msg);
% OTFAttwdth
msg="Should be a nonnegative real.";
assert(isfield(params, "OTFAttwdth"),prefix + missg + "`OTFAttwdth`. " + msg);
assert(params.OTFAttwdth>=0,prefix + invld + "`OTFAttwdth`. " + msg);

% -- Operators, costs, and optim
prefix="[Image Reconstruction] ";
% mu, stepTol
msg="Should be a positive real.";
assert(isfield(params, "mu"),prefix + missg + "`mu`. " + msg);
assert(params.mu> 0,prefix + invld + "`mu`. " + msg);
assert(isfield(params, "stepTol"),prefix + missg + "`stepTol`. " + msg);
assert( params.stepTol > 0,prefix + invld + "`stepTol`. " + msg);
% padSZ, maxIt
msg="Should be a non-negative integer.";
assert(isfield(params, "maxIt"),prefix + missg + "`maxIt`. " + msg);
assert(params.maxIt>= 0,prefix + invld + "`maxIt`. " + msg);
assert(isfield(params, "padSz"),prefix + missg + "`padSz`. " + msg);
assert(params.padSz>= 0,prefix + invld + "`padSz`. " + msg);
% sepOrr apodize
msg="Should be a boolean.";
assert(isfield(params, "sepOrr"),prefix + missg + "`sepOrr`. " + msg);
assert(params.sepOrr==0 || params.sepOrr==1, prefix + invld + "`sepOrr'. " + msg);  
assert(isfield(params, "apodize"),prefix + missg + " apodize. " + msg);
assert(params.apodize==0 || params.apodize==1, prefix + invld + "`apodize'. " + msg);  
% regType
msg="Should be an integer within {1, 2, 3}.";
assert(isfield(params, "regType"),prefix + missg + "`regType`. " + msg);
assert(mod(params.regType, 1) == 0 && params.regType > 0 && params.regType < 4,prefix + invld + "`regType`. " + msg);














