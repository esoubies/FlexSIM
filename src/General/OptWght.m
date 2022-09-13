function a = OptWght(x,xref,mask)
%--------------------------------------------------------------------------
% Function a = optWght(x,xref,mask)
%
% Takas as inputs two arrays and a mask, and calculates the parameter `a` 
% that will minimize the objective |x - xref|^2
%
% Inputs : ft        → Fourier transform (of SIM image) to mask 
%          fc        → Cutoff frequency of OTF in [nm/rad]
%          limits    → 1x2 matrix. Limits of the ring as factors of the
%                      maximum cutoff frequency. 
% Outputs : maskedFT → Self-explanatory 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% If no mask was provided, define one that covers the whole image
if nargin < 3
    sz = size(x);
    mask = ones(sz); 
end

x=x.*mask; xref=xref.*mask;         % Multiply both arrays by mask
a=(x(:)'*xref(:))/norm(x(:))^2;     % Solve system to get parameter `a`  
end