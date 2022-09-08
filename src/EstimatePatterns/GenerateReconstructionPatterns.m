function [patt, params] = GenerateReconstructionPatterns(params, displ)
%--------------------------------------------------------------------------
%   [patterns, params] = GenerateReconstructionPatterns(filename, params, displ)
%
% Inputs : filename : Raw SIM data 
%          params   : parameters with stored data. Should contain at least
%                     the fields k and phase.
%          displ    : Whereas to display a sample pattern
%
% Function that from the etimated parameters of a raw set of SIM images,
% returns the estimated patterns
% 
%--------------------------------------------------------------------------
% Extracting variables
if isfield(params,'roi')
    lx=-(params.roi(2)-1)*params.res:params.res/2:(params.sz(2)-params.roi(2)+1)*params.res-params.res/2;
    ly=-(params.roi(1)-1)*params.res:params.res/2:(params.sz(1)-params.roi(1)+1)*params.res-params.res/2;
    [X,Y]=meshgrid(lx,ly); 
else
    [X,Y]=meshgrid(0:2*params.sz(2)-1,0:2*params.sz(1)-1); X=X*params.res/2; Y=Y*params.res/2;
end
patt = zeros(2*params.sz(1), 2*params.sz(2), params.nbOr*params.nbPh); 
for i = 1:params.nbOr
    k = params.k(i, :); 
    for j = 1:params.nbPh
        if params.method == 2
           phoff = params.phase(i);
           a = params.a(i);
           patt(:,:,(i-1)*params.nbPh+j) = 1 + a*(cos(2*(phoff+(j-1)*pi/params.nbPh))*cos(2*(k(1)*X+k(2)*Y)) ...
                                                             - sin(2*(phoff+(j-1)*pi/params.nbPh))*sin(2*(k(1)*X+k(2)*Y)));
        else
        a = params.a(i, j); 
        patt(:,:,(i-1)*params.nbPh+j) = 1 + a*(cos(2*ph)*cos(2*(k(1)*X+k(2)*Y)) ...
                                                          - sin(2*ph)*sin(2*(k(1)*X+k(2)*Y)));
        end
    end
end
patt=patt/(mean(patt(:))*size(patt,3));

if displ
    figure; 
    for i = 1:params.nbOr
        subplot(2, ceil(params.nbOr/2), i); imshow(patt(:,:,i*params.nbPh), []); 
        title("Sample pattern of orientation " + num2str(i));
    end
    impixelinfo;
end
end