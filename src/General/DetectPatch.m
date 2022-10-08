function pos=DetectPatch(im,sz_patch,minmax)
%--------------------------------------------------------------------------
% function pos=DetectPatch(im,sz_patch,minmax)
%
% Search the patch of size sz_patch (odd number) with minimal (minmax=-1)
% or maximal (minmax=1) intensity within the image im
%
% Retuns the coordinates of the top left corner of the detected patch
%--------------------------------------------------------------------------

k=ones(sz_patch);
c=conv2(im,k,'valid');
if minmax==1
    [~,idx]=max(c(:));
elseif minmax==-1
    [~,idx]=min(c(:));
end
[i,j]=ind2sub(size(c),idx);
pos=[i,j];
end
