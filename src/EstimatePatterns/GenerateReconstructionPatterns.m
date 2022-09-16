function [patt] = GenerateReconstructionPatterns(params, res, y)
%--------------------------------------------------------------------------
%  function [patt] = GenerateReconstructionPatterns(params, res)
%
% Takes as input the etimated pattern parameters (amplitude, phase and 
% wavevector) of a SIM acquisition and returns the corresponding patterns
%
% Inputs :  params   -> Parameters with stored data. Should contain at least
%                       the fields k, phase and a, as well as the raw SIM data.
%           res      -> Structure containing the wave vectors, phases, and 
%                       amplitudes estimated
%           y        -> Raw SIM data
% 
% Outputs : patt     -> Estimated patterns of the SIM acquisition, of the
%                       same size as the data
%
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Extracting variables
sz = size(y); 
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