function a = OptWght(x,xref,mask)
%--------------------------------------------------------------------------
% Function a = OptWght(x,xref,mask)
%
% Takes as inputs two arrays and a mask, and calculates the parameter `a` 
% that will minimize the objective `|a*x - xref|^2`
%
% Inputs : x      -> Array to multiply by `a`
%          xref   -> Reference array
%          limits -> 1x2 matrix. Limits of the ring as factors of the
%                      maximum cutoff frequency. 
% Outputs : a     -> Factor that minimizes least squares 
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

sz = size(x);
% If no mask was provided, define one that covers the whole image
if nargin < 3, mask = ones(sz);  end
% Solution to |a*x - xref|^2
x=x.*mask;
xref=xref.*mask;
a=sum(conj(x).*xref,[1,2])./sum(abs(x).^2,[1,2]);

end