function k = ExtractLocMin(params,Jp,K1,K2)
%--------------------------------------------------------------------------
% function k = ExtractLocMin(params,Jp,K1,K2)
%
% Extract the params.nMinima smallest local minima of Jp (evaluation of the
% function J over the grid of wavevectors defined by K1 and K2, see
% GridEvalJ.m)
%
% See also EstimatePatterns.m GridEvalJ.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Initializations
k = zeros(params.nMinima, 2);   % Give user info and initialize...

% -- Extract n local mins
for nth = 1:params.nMinima
    [~,idxMin]=min(Jp(:));                     % Extract and store minima
    [ii,jj]=ind2sub([params.nPoints,params.nPoints/2],idxMin);
    k(nth, :) = [K1(ii,jj);K2(ii,jj)];
    
    % Temporarily erase the vecinity of the found minima
    idx1 = ii - round(params.nPoints/50); idx1(idx1<1) = 1;
    idx2 = ii + round(params.nPoints/50); idx2(idx2<params.nPoints) = params.nPoints;
    idx3 = jj - round(params.nPoints/50); idx3(idx3<1) = 1;
    idx4 = jj + round(params.nPoints/50); idx4(idx4>params.nPoints/2) = params.nPoints/2;
    Jp(idx1:idx2, idx3:idx4)  = max(Jp(:));
end

end