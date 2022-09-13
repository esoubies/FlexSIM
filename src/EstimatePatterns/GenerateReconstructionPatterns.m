function [patt, params] = GenerateReconstructionPatterns(params, res)
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
sz = size(params.y); 
if isfield(params,'roi') && ~isempty(params.roi)
    lx=-(params.roi(2)-1)*params.res:params.res/2:(sz(2)-params.roi(2)+1)*params.res-params.res/2;
    ly=-(params.roi(1)-1)*params.res:params.res/2:(sz(1)-params.roi(1)+1)*params.res-params.res/2;
    [X,Y]=meshgrid(lx,ly); 
else
    [X,Y]=meshgrid(0:2*sz(2)-1,0:2*sz(1)-1); X=X*params.res/2; Y=Y*params.res/2;
end
patt = zeros(2*sz(1), 2*sz(2), params.nbOr*params.nbPh); 
for i = 1:params.nbOr
    k = res.k(i, :); 
    for j = 1:params.nbPh
        if params.method == 2
           phoff = res.phase(i);
           a = res.a(i);
           patt(:,:,(i-1)*params.nbPh+j) = 1 + a*(cos(2*(phoff+(j-1)*pi/params.nbPh))*cos(2*(k(1)*X+k(2)*Y)) ...
                                                 - sin(2*(phoff+(j-1)*pi/params.nbPh))*sin(2*(k(1)*X+k(2)*Y)));
        else
            a = res.a(i, j); ph = res.phase(i, j);
            patt(:,:,(i-1)*params.nbPh+j) = 1 + a*(cos(2*ph)*cos(2*(k(1)*X+k(2)*Y)) ...
                                                  - sin(2*ph)*sin(2*(k(1)*X+k(2)*Y)));
        end
    end
end
patt=patt/(mean(patt(:))*size(patt,3));
end