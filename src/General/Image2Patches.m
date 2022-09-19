function patches = Image2Patches(im,sz_patch,overlap)
%--------------------------------------------------------------------------
% function patches = Image2Patches(im,sz_patch,overlap)
%
% Extract 2D patches from an image
%
% Inputs: im       -> input image (can be a stack)
%         sz_patch -> size [in px] of (square) patchs
%         overlap  -> overlap [in px] between patches
%
% Output: patches  -> cell array containing extracted patches
%
% See also Patches2Image.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
sz=size(im);
I=ceil((sz(1)-overlap)/(sz_patch-overlap));
J=ceil((sz(2)-overlap)/(sz_patch-overlap));
step=sz_patch-overlap;

% -- Extract patches
patches=cell(I,J);
for i=1:I
    for j=1:J
        idx_b=1 + (i-1)*step; 
        idx_e=min(sz_patch + (i-1)*step,sz(1));
        jdx_b=1 + (j-1)*step; 
        jdx_e=min(sz_patch + (j-1)*step,sz(2));
        patches{i,j}=im(idx_b:idx_e,jdx_b:jdx_e,:);
    end
end

end
