function masked_ft = MaskFT(ft, fc, res, limits, display)
%--------------------------------------------------------------------------
%   k_est = LandscapeJ(G,wf,fc, OTF, n_pts, res, zoom, grad, display)
%
% Inputs : ft   : Fourier transform (of SIM image) to mask 
%          fc   : Cutoff frequency of OTF in [nm/rad]
%          res  : Resolution of the image in [nm]
%          limits : 1x2 matrix. Limits of the ring as factors of the system fc
%          display : Whether to display findings or not
%
% MaskFT takes as input the FT (shifted) of a SIM image, the cutoff
% frequency of the system, and the limits of the ring that will mask the
% image as a product of the fc 
%--------------------------------------------------------------------------

sz = size(ft);           % Size of the image. 
maxp=fc*sz(1);         % Radius (in pxls) of the largest possible wavevector
if nargin < 3
    limits = [0.8, 1.1];
    display = 0;
else if nargin < 4
    display = 0;
end

[I, J] = meshgrid(0:sz(1)-1,0:sz(2)-1);
circ1 = sqrt((I-floor(sz(1)/2)-1).^2 + (J-floor(sz(2)/2)-1).^2) < maxp*limits(1);
circ2 = sqrt((I-floor(sz(1)/2)-1).^2 + (J-floor(sz(2)/2)-1).^2) < maxp*limits(2);
ring = circ2 - circ1;

masked_ft = ft.*ring; 


end