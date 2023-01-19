function maskedFT = MaskFT(ft, fc, radii,shift)
%--------------------------------------------------------------------------
% Function masked_ft = MaskFT(ft, fc, radii,shift)
%
% Takes as input a FT (propoerly shifted) and masks it by multiplying it 
% (element-wise) with a boolean ring. The ring is parametrized in terms of
% multiples of the maximum cutoff frequency of the system (provided as input)
%
% Inputs : ft        -> Fourier transform (of SIM image) to mask 
%          fc        -> Cutoff frequency, normalized by resolution ([0, 0.5])
%          limits    -> 1x2 matrix. Limits of the ring as factors of the
%                      maximum cutoff frequency. 
%          shift     -> optional parameter to shift the center (0,0) of the mask
% Outputs : maskedFT -> Self-explanatory 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

if nargin <4
    shift=[0,0];
end

sz = size(ft);            % Get size of the image. 
maxp=fc*sz;          % Radius (in pxls) of the largest possible wavevector
[I, J] = meshgrid(1:sz(2),1:sz(1));    % Create grid and boolean circles
I=I-shift(1); J=J-shift(2);
ellips1 = sqrt(((I-floor(sz(2)/2)-1).^2)./(maxp(2).^2) + ((J-floor(sz(1)/2)-1).^2)./(maxp(1).^2)) < radii(1); 
ellips2 = sqrt(((I-floor(sz(2)/2)-1).^2)./(maxp(2).^2) + ((J-floor(sz(1)/2)-1).^2)./(maxp(1).^2)) < radii(2); 
mask = ellips2 - ellips1;     % Define boolean mask
maskedFT = ft.*mask;      % Mask
end