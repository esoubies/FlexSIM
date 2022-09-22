function patt = GenerateReconstructionPatterns(params,k,phase,a,sz)
%--------------------------------------------------------------------------
% function patt = GenerateReconstructionPatterns(params,k,phase,a,sz)
%
% Sample a 2D sinusoidal pattern
%    w(x) = 1 + a cos(<k,x> + phase)
% from its parameters (amplitude, phase and  wavevector).
%
% Inputs :  params -> Structures with fields:
%                         - roi: region on interest on which the parameters  have been estimated (see function EstimatePatterns)
%                         - res: resolution of the SIM data stack (patterns will be generated on a twice finer grid)
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%                         - method: method used to estimate the parameters (see function EstimatePatterns)
%           k      -> Array (nbOr x 2) containing the wavevector for each orientation
%           phase  -> Array containing the absolute phases (method = 0 or 1) or phase offsets (method =2)
%           a      -> Array containing the amplitudes of the patterns
%           sz     -> Size of the raw SIM stack
% 
% Outputs : patt   -> Sampled patterns
%
% See also EstimatePatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Extracting variables
if isfield(params,'roi') && ~isempty(params.roi)
    lx=-(params.roi(2)-1)*params.res:params.res/2:(sz(2)-params.roi(2)+1)*params.res-params.res/2;
    ly=-(params.roi(1)-1)*params.res:params.res/2:(sz(1)-params.roi(1)+1)*params.res-params.res/2;
    [X,Y]=meshgrid(lx,ly); 
else
    [X,Y]=meshgrid(0:2*sz(2)-1,0:2*sz(1)-1); X=X*params.res/2; Y=Y*params.res/2;
end
patt = zeros(2*sz(1), 2*sz(2), params.nbOr*params.nbPh); 
if params.method
    for i = 1:params.nbOr
        ki = k(i, :); 
        for j = 1:params.nbPh
            if params.method == 2
               phoff = phase(i);
               ai = a(i);
               patt(:,:,(i-1)*params.nbPh+j) = 1 + ai*(cos(2*(phoff+(j-1)*pi/params.nbPh))*cos(2*(ki(1)*X+ki(2)*Y)) ...
                                                     - sin(2*(phoff+(j-1)*pi/params.nbPh))*sin(2*(ki(1)*X+ki(2)*Y)));
            else
                aij = a(i, j); ph = phase(i, j);
                patt(:,:,(i-1)*params.nbPh+j) = 1 + aij*(cos(2*ph)*cos(2*(ki(1)*X+ki(2)*Y)) ...
                                                      - sin(2*ph)*sin(2*(ki(1)*X+ki(2)*Y)));
            end
        end
    end
else
    for i = 1:params.nbOr*params.nbPh
        ki = k(i, :);        
        phoff = phase(i);
        ai = a(i);
        patt(:,:,i) = 1 + ai*(cos(2*(phoff*pi/params.nbPh))*cos(2*(ki(1)*X+ki(2)*Y)) ...
                            - sin(2*(phoff*pi/params.nbPh))*sin(2*(ki(1)*X+ki(2)*Y)));
       
    end
end
patt=patt/(mean(patt(:))*size(patt,3));
end