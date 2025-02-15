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
% [1] Handling Challenging Structured Illumination Microscopy Data with FlexSIM
%     E. Soubies et al, Preprint, 2023
%
% See also GenerateReconstructionPatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
sz=size(y);
h = gpuCpuConverter(fspecial('gaussian',min(min(sz(1:2)),6*filt_sz),filt_sz));
y = gpuCpuConverter(y);
wf_stack = gpuCpuConverter(wf_stack);

% -- Loop over orientations and phases
Lf=zeros_([sz(1:2)*2,sz(3:end)]);
for it=1:size(y,4)
    for ii=1:params.nbOr
        for jj=1:params.nbPh
            id=(ii-1)*params.nbPh+jj;
            wf_i=wf_stack(:,:,min(ii,size(wf_stack,3)),it);
            y_i=y(:,:,id,it);
            y_filt = imfilter(y_i,h,'symmetric')./max(imfilter(wf_i,h,'symmetric'),eps);
            Lf(:,:,id,it)=imresize(y_filt,sz(1:2)*2);
        end
    end
end


end
