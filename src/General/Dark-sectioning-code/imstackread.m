function ImStack = imstackread(img)
% This function can import 3D image stack.
% The bit depth is limited to 8,16,32-bit
% Written by Ethan Zhao, Jan. 2021
% Tutorial: https://zhuanlan.zhihu.com/p/326688179/

imgInfo = imfinfo(img);
imgRow = imgInfo(1).Height;
imgCol = imgInfo(1).Width;
imgDepth = length(imgInfo);
imgBitDepth = ['uint',num2str(imgInfo(1).BitDepth)];
ImStack = zeros(imgRow, imgCol, imgDepth,imgBitDepth);
for ii = 1 : length(imgInfo)
    ImStack(:,:,ii) = imread(img, ii);
end

end
