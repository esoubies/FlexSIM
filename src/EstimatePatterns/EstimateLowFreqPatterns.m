function Lf = EstimateLowFreqPatterns(params,y,wf_stack,filt_sz)
%--------------------------------------------------------------------------
% function Lf = EstimateLowFreqPatterns(params,y,wf_stack,filt_sz)
%
% Estimate a non-constant low-frequency componnent (Lf) of a 2D sinusoidal
% pattern of the form:
%    w(x) = Lf(x) + a cos(<k,x> + phase)
%
% Inputs: params   -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)
%         y        -> SIM data stack
%         wf_stack -> stack of wf images (either 1 per orientation or only 1)
%         filt_sz  -> size (in px) of the Gaussian filter used to extract low-freqs
%
% Output: Lf      -> Extracted low frequency component of the pattern
%
% See also GenerateReconstructionPatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
sz=size(y);
h = fspecial('gaussian',min(min(sz(1:2)),6*filt_sz),filt_sz);

% -- Loop over orientations and phases
Lf=zeros(sz(1:2)*2);
for ii=1:params.nbOr
    for jj=1:params.nbPh
        id=(ii-1)*params.nbPh+jj;
        y_filt = imfilter(y(:,:,id),h,'symmetric')./imfilter(wf_stack(:,:,min(ii,size(wf_stack,3))),h,'symmetric');
        Lf(:,:,id)=imresize(y_filt,sz(1:2)*2);
    end
end


end