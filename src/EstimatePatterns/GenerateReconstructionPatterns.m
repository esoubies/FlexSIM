function patt = GenerateReconstructionPatterns(params,PosRoiPatt,k,phase,a,sz,Lf)
%--------------------------------------------------------------------------
% function patt = GenerateReconstructionPatterns(params,PosRoiPatt,k,phase,a,sz,Lf)
%
% Sample a 2D sinusoidal pattern
%    w(x) = Lf(x) + a cos(<k,x> + phase)
% from its parameters (amplitude, phase and  wavevector) and the
% low-frequency component Lf.
%
% Inputs :  params    -> Structures with fields:
%                         - SzRoiPatt: Size of the ROI used for patterns estimation
%                         - res: resolution of the SIM data stack (patterns will be generated on a twice finer grid)
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%                         - method: method used to estimate the parameters (see function EstimatePatterns)
%           SzRoiPatt -> Top-left position of the ROI used for patterns estimation
%           k         -> Array (nbOr x 2) containing the wavevector for each orientation
%           phase     -> Array containing the absolute phases (method = 0 or 1) or phase offsets (method =2)
%           a         -> Array containing the amplitudes of the patterns
%           sz        -> Size of the raw SIM stack
%           Lf_comp   -> Low-freq component of the patterns, default 1 (see EstimateLowFreqPatterns.m)
% 
% Outputs : patt   -> Sampled patterns
%
% See also EstimatePatterns.m EstimateLowFreqPatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Extracting variables
if isfield(params,'SzRoiPatt') && ~isempty(params.SzRoiPatt)
    lx=-(PosRoiPatt(2)-1)*params.res:params.res/2:(sz(2)-PosRoiPatt(2)+1)*params.res-params.res/2;
    ly=-(PosRoiPatt(1)-1)*params.res:params.res/2:(sz(1)-PosRoiPatt(1)+1)*params.res-params.res/2;
    [X,Y]=meshgrid(lx,ly); 
else
    [X,Y]=meshgrid(0:2*sz(2)-1,0:2*sz(1)-1); X=X*params.res/2; Y=Y*params.res/2;
end
patt = zeros(2*sz(1), 2*sz(2), params.nbOr*params.nbPh,size(Lf,4)); 
X = gpuCpuConverter(X); 
Y = gpuCpuConverter(Y);
patt = gpuCpuConverter(patt);


for it=1:size(Lf,4)
    for i = 1:params.nbOr
        ki = k(i, :);
        for j = 1:params.nbPh
            if params.eqPh
                phoff = phase(i);
                patt(:,:,(i-1)*params.nbPh+j,it) = 1 + a*(cos(2*(phoff+(j-1)*pi/params.nbPh))*cos(2*(ki(1)*X+ki(2)*Y)) ...
                    - sin(2*(phoff+(j-1)*pi/params.nbPh))*sin(2*(ki(1)*X+ki(2)*Y)));
            else
                ph = phase(i, j);
                patt(:,:,(i-1)*params.nbPh+j,it) = 1 + a*(cos(2*ph)*cos(2*(ki(1)*X+ki(2)*Y)) ...
                    - sin(2*ph)*sin(2*(ki(1)*X+ki(2)*Y)));
            end
        end
    end
end

patt=patt-mean(patt,[1,2]) +Lf.*mean(patt,[1,2])./mean(Lf,[1,2]);
patt=patt/(mean(patt(:))*size(patt,3));
end
