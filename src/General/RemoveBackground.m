function [y]  = RemoveBackground(y, params)
%--------------------------------------------------------------------------
% Function y = RemoveBackground(y, params)
% 
% Estimate a constant background from an ROI that does not contains objects
% and subtract it to the stack. 
%
% Inputs : y       -> Raw SIM data
%          params  -> Structures with acquisition fields (see EstimatePatterns 
%                     for details). 
%
% Outputs: y       -> Raw SIM stack without background
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

if isfield(params,'roiBack') && ~isempty(params.roiBack)
    roi=y(params.roiBack(1):params.roiBack(1)+params.roiBack(3)-1,params.roiBack(2):params.roiBack(2)+params.roiBack(3)-1,:);
    y=max(y-mean(roi,[1 2]),0);
end

end