function im = Patches2Image(patches,overlap)
%--------------------------------------------------------------------------
% function im = Patches2Image(patches,overlap)
%
% Reconstruct image from 2D patches
%
% Inputs: patches  -> cell array of patches (see function Image2Patches)
%         overlap  -> overlap [in px] between patches
%
% Output: im       -> reconstructed image
%
% See also Image2Patches.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Pre-computations
[I,J]=size(patches);
sz_patch=size(patches{1},1);
step=sz_patch-overlap;
% Size of the output image
sz_1=(I-1)*step+size(patches{I,1},1);
sz_2=(J-1)*step+size(patches{1,J},2);
sz_3=size(patches{1},3);
% Apodisation function
[X,Y]=meshgrid(linspace(-1,1,sz_patch),linspace(-1,1,sz_patch));
Apo=1./(1+exp(30*(abs(X)-0.8)))./(1+exp(30*(abs(Y)-0.8)));

% -- Reconstruct image
im=zeros([sz_1,sz_2,sz_3]);
ww=zeros([sz_1,sz_2]);
for i=1:I
    for j=1:J
        if norm(patches{i,j}(:))>0
            idx_b=1 + (i-1)*step;
            idx_e=min(sz_patch + (i-1)*step,sz_1);
            jdx_b=1 + (j-1)*step;
            jdx_e=min(sz_patch + (j-1)*step,sz_2);
            if isequal(size(patches{i,j}),size(Apo))
                im(idx_b:idx_e,jdx_b:jdx_e,:) = im(idx_b:idx_e,jdx_b:jdx_e,:) + patches{i,j}.*Apo;
                ww(idx_b:idx_e,jdx_b:jdx_e)=ww(idx_b:idx_e,jdx_b:jdx_e)+Apo;
            else
                [X,Y]=meshgrid(linspace(-1,1,size(patches{i,j},2)),linspace(-1,1,size(patches{i,j},1)));
                Apo_bis=1./(1+exp(30*(abs(X)-0.8)))./(1+exp(30*(abs(Y)-0.8)));
                im(idx_b:idx_e,jdx_b:jdx_e,:) = im(idx_b:idx_e,jdx_b:jdx_e,:) + patches{i,j}.*Apo_bis;
                ww(idx_b:idx_e,jdx_b:jdx_e)=ww(idx_b:idx_e,jdx_b:jdx_e)+Apo_bis;
            end
        end
    end
end
ww(ww==0)=1;
im=im./ww;

end


