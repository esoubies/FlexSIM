function maskedFT = MaskFT(ft, fc, radii)
%--------------------------------------------------------------------------
% Function masked_ft = MaskFT(ft, fc, res, radii)
%
% Takes as input a FT (propoerly shifted) and masks it by multiplying it 
% (element-wise) with a boolean ring. The ring is parametrized in terms of
% multiples of the maximum cutoff frequency of the system (provided as input)
%
% Inputs : ft        → Fourier transform (of SIM image) to mask 
%          fc        → Cutoff frequency, normalized by resolution ([0, 0.5])
%          limits    → 1x2 matrix. Limits of the ring as factors of the
%                      maximum cutoff frequency. 
% Outputs : maskedFT → Self-explanatory 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

sz = size(ft);            % Get size of the image. 
maxp=fc*min(sz);          % Radius (in pxls) of the largest possible wavevector
[I, J] = meshgrid(0:sz(1)-1,0:sz(2)-1);    % Create grid and boolean circles
circ1 = sqrt((I-floor(sz(1)/2)-1).^2 + (J-floor(sz(2)/2)-1).^2) < maxp*radii(1);
circ2 = sqrt((I-floor(sz(1)/2)-1).^2 + (J-floor(sz(2)/2)-1).^2) < maxp*radii(2);
ring = circ2 - circ1;     % Define boolean mask
maskedFT = ft.*ring;      % Mask
end