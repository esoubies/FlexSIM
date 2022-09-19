function Lf = EstimateLowFreqPatterns(y,filt_sz)
%--------------------------------------------------------------------------
% function Lf = EstimateLowFreqPatterns(y,filt_sz)
%
%
%
% Inputs: y       -> SIM data stack
%         filt_sz -> size (in px) of the Gaussian filter used to extract low-freqs
%
% Output: Lf      -> Extracted low frequency component of the pattern
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

sz=size(y);
h = fspecial('gaussian',min(min(sz(1:2)),128),filt_sz);
% y_filt = imfilter(y./mean(y,3),h,'symmetric');
y_filt = imfilter(y,h,'symmetric')./imfilter(mean(y,3),h,'symmetric');
Lf=zeros(size(y_filt(:,:,1))*2);
for kk=1:size(y_filt,3)
    Lf(:,:,kk)=imresize(y_filt(:,:,kk),size(y_filt(:,:,kk))*2);
end
%                 patterns=patterns-mean(patterns,[1,2]) +tt.*mean(patterns,[1,2])./mean(tt,[1,2]);
%                 patterns=patterns/(mean(patterns(:))*size(patterns,3));
end