function  [Jp,K1,K2] = GridEvalJ(params,wf,G,grid_data)
%--------------------------------------------------------------------------
% function  [Jp,K1,K2] = GridEvalJ(params,wf,G,grid_data)
%
% Evaluate the function J over a grid
%
% Inputs  : params    -> Structures with fields:
%                         - lamb: Emission wavelength
%                         - Na: Objective numerica aperture
%                         - res: resolution of the SIM data stac
%                         - limits: 1x2 array (eg. [0.9, 1.1]) defining  the ring over which the J function is evaluated for initializing, givien as factor of fc=1
%                         - nPoints: Number of points in the J evaluation grid
%                         - method: method used to estimate the parameter (see EstimatePatterns.m)
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%           wf        -> widefield image
%           G         -> pre-processed data (output of RemoveWFandMask  function)
%           grid_data -> grid corresponding to the resolution of the data
%      
% Outputs:  Jp      -> grid with evaluation of J
%           K1      -> corresponding first component of wavevector 
%           K2      -> corresponding second component of wavevector 
%
% See also EstimatePatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Precomputations
FCut = 2*params.Na/params.lamb*params.res;            % Cut-off frequency
FCutN=FCut*pi/params.res;                             % Scaled cutoff frequency
maxp=FCut*pi/params.res*max(1,params.limits(2));      % Maximum value of one component of k
ll1=linspace(0,maxp,params.nPoints/2);                % Declare grid of wavevectors...
ll2=linspace(-maxp,maxp,params.nPoints);              % with the x direction halved 
[K1,K2]=meshgrid(ll1,ll2);                            
K = cat(3, K1, K2); 
Knorm = vecnorm(K, 2, 3);                             % Calculate norms in advance

% -- Initializations
if params.GPU
    Jp = zeros(params.nPoints, params.nPoints/2,'double','gpuArray');  % Initialize Jp for batch
else
    Jp = zeros(params.nPoints, params.nPoints/2);         % Initialize Jp for batch
end

% -- Evaluation of J over the grid
for ii=1:params.nPoints                           % Iterate the grid and if the wavevector
    for jj=1:params.nPoints/2                     % is within the limits, evaluate
        if Knorm(ii, jj)>params.limits(1)*FCutN && Knorm(ii, jj)<params.limits(2)*FCutN
            ktest = squeeze(K(ii, jj, :));        % Evaluate cost without filtering
            Jp(ii, jj) = EvalJ(ktest, wf, G, params, grid_data, 0, 0, 0);
        else
            Jp(ii,jj)=NaN;
        end
    end
end

end