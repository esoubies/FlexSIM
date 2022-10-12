function Lf = EstimateLowFreqPatterns(y,filt_sz)
%--------------------------------------------------------------------------
% function Lf = EstimateLowFreqPatterns(y,filt_sz)
%
% Estimate a non-constant low-frequency componnent (Lf) of a 2D sinusoidal
% pattern of the form:
%    w(x) = Lf(x) + a cos(<k,x> + phase)
%
% Inputs: y       -> SIM data stack
%         filt_sz -> size (in px) of the Gaussian filter used to extract low-freqs
%
% Output: Lf      -> Extracted low frequency component of the pattern
%
% See also GenerateReconstructionPatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

sz=size(y);
h = fspecial('gaussian',min(min(sz(1:2)),6*filt_sz),filt_sz);
% y_filt = imfilter(y./mean(y,3),h,'symmetric');
y_filt = imfilter(y,h,'symmetric')./imfilter(mean(y,3),h,'symmetric');
 
Lf=zeros(size(y_filt(:,:,1))*2);
for kk=1:size(y_filt,3)
    Lf(:,:,kk)=imresize(y_filt(:,:,kk),size(y_filt(:,:,kk))*2);
end

end