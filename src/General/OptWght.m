function a = OptWght(x,xref,p,r)
%--------------------------------------------------------------------------
% function a = OptWght(x,xref,p,r)
%
% Computes a = argmin ||mask*(a*x-xref)||^2
% where mask is a disk of radius r centered at p
%
% Input : x    : image to be scaled
%         xref : reference image
%         p    : mask center
%         r    : mask radius
%
% Output: a  : Regression coefficient
%--------------------------------------------------------------------------

% -- Pre computations
sz = size(x);

% -- Create mask
[I,J]=meshgrid(1:sz(1),1:sz(2));
mask=double(((I-p(1)).^2+(J-p(2)).^2) < r^2); 

% -- Compute weight
x=x.*mask;
xref=xref.*mask;
a=(x(:)'*xref(:))/norm(x(:))^2;

end