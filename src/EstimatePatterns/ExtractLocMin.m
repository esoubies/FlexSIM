function [k,ii,jj] = ExtractLocMin(n,map,K1,K2)
%--------------------------------------------------------------------------
% function k = ExtractLocMin(n,map,K1,K2)
%
% Extract the n smallest local minima of Jp (evaluation of the
% function map over the grid of wavevectors defined by K1 and K2)
%
% See also EstimatePatterns.m GridEvalJ.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Initializations
k = zeros(n, 2);   % Give user info and initialize...
npts1=size(K1,1);
npts2=size(K1,2);
szroi=max(npts1,npts2);

% -- Extract n local mins
for nth = 1:n
    [~,idxMin]=min(map(:));                     % Extract and store minima
    [ii,jj]=ind2sub([npts1,npts2],idxMin);
    k(nth, :) = [K1(ii,jj);K2(ii,jj)];
    
    % Temporarily erase the vecinity of the found minima
    idx1 = ii - round(szroi/50); idx1(idx1<1) = 1;
    idx2 = ii + round(szroi/50); idx2(idx2<npts1) = npts1;
    idx3 = jj - round(szroi/50); idx3(idx3<1) = 1;
    idx4 = jj + round(szroi/50); idx4(idx4>npts2) = npts2;
    map(idx1:idx2, idx3:idx4)  = max(map(:));
end

end