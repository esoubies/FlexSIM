function y = RemoveBackground(y,PosRoiBack,SzRoi)
%--------------------------------------------------------------------------
% Function y = RemoveBackground(y,PosRoiBack,SzRoi)
% 
% Estimate a constant background from an ROI that does not contains objects
% and subtract it to the stack. 
%
% Inputs : y          -> Raw SIM data
%          PosRoiBack -> Top left corner of the ROI
%          SzRoi      -> Size of the ROI
%
% Outputs: y       -> Raw SIM stack without background
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

roi=y(PosRoiBack(1):PosRoiBack(1)+SzRoi-1,PosRoiBack(2):PosRoiBack(2)+SzRoi-1,:);
for ii = 1: size(y,3)
    y(:,:,ii) = max(y(:,:,ii) - mean(mean(roi(:,:,ii),1),2),0);
end

end