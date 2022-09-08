function [y, wfs, params] = Preprocess(filename, psfname, params, displ)
%--------------------------------------------------------------------------
% function [y, wfs, params] = Preprocess(filename, params, displ)
%
% Inputs :
%
% 
%--------------------------------------------------------------------------
%% Signal Conditioning & Preprocessing
params.fc = 2*params.Na/params.lamb*params.res;

y = double(loadtiff(filename)); params.origSz=size(y); % Load and keep size before crop

if displ > 1
    figure(); subplot(2,3,1); imshow(y(:,:,1), []); title("Sample raw image"); colorbar; colormap(viridis);
    subplot(2,3,4); imshow(abs(log10(fftshift(fft2(y(:,:,1)))+1)), []); title("Raw image FT"); colorbar; colormap(viridis);
end

% First we get the ROI and/or ensure that the stack is squared
if isfield(params,'roi') && ~isempty(params.roi)   % Crop to ROI
    y=y(params.roi(1):params.roi(1)+params.roi(3)-1,params.roi(2):params.roi(2)+params.roi(3)-1,:);
elseif params.sz(1) ~= params.sz(2) % Make square centered
    dim = min(params.sz(1), params.sz(2));
    center = [floor(params.sz(1)/2)+1, floor(params.sz(2)/2)+1]; 
    y = y(center(1)-floor(dim/2):center(1)+floor(dim/2)-1, center(2)-floor(dim/2):center(2)+floor(dim/2)-1, :);     
end
params.sz = size(y);         % Recalc size if it was modified

%% Process y stack - Calculate and substract WF, mask with rings on ROIs
yTmp = zeros(params.sz);                      % Initialize variables
wfs =zeros(params.sz(1), params.sz(2),params.nbOr);
y = (y - min(y(:))) / (max(y(:)) - min(y(:))); % Normalize images

% Select the indexes to use. If Zeiss convention (zap), flip
params.imgIdxs = 1:params.nbPh*params.nbOr; 
params.imgIdxs = reshape(params.imgIdxs, [params.nbPh, params.nbOr]); 
if strcmp(params.AcqConv, "zap")  
    params.imgIdxs = params.imgIdxs'; 
end
% If there are any assumptions on multiple images, we will use all the phases
if params.method > 0
    params.nbImgs = params.nbPh; 
else 
    % Else we will only use the first image of each orientation
    params.imgIdxs = params.imgIdxs(1,:); 
    params.nbImgs = 1;
end

% Remove widefield and mask
idwf=0;
for idx = params.imgIdxs
    idwf=idwf+1;
    % Compute the wf and fft per batch of images (orientation)
    wfs(:,:,idwf)=sum(y(:,:,idx), 3)/length(idx); fft_wf = fftshift(fft2(wfs(:,:,idwf)));
    % Treat images individually - get FFT, remove scaled fftWF & mask
    for img = idx' 
        y_tmp = y(:,:,img); ffty=fftshift(fft2(y_tmp));   % image and FFT data
        a=real(OptWght(ffty,fft_wf,(params.sz(1:2)+1)/2, params.sz(1)/8)); % fc * 0.2                   
        if isfield(params, "ringMaskLim")
            % Remove the scaled WF + filter freq within a ring of interest
            yTmp(:,:,img) = real(ifft2(ifftshift(MaskFT((a*ffty-fft_wf), params.fc, params.res, params.ringMaskLim, 1)))); 
        else
            % Remove the scaled wf + constrain to OTF supportn
            yTmp(:,:,img) = real(ifft2(ifftshift(params.OTF0.*(a*ffty-fft_wf)))); 
        end
    end
end

% Mask widefield
if isfield(params, "ringMaskLim") 
%     wfs = real(ifft2(ifftshift(MaskFT(fftshift(fft2(wfs)), params.fc, params.res, params.ringMaskLim, 1))));
    for ii=1:size(wfs,3)
        wfs(:,:,ii) = real(ifft2(ifftshift(MaskFT(fftshift(fft2(wfs(:,:,ii))), params.fc, params.res, params.ringMaskLim, 1))));
    end
end
y = yTmp; 
%% Create and save PSF
% [PSF,OTF] = GeneratePSF(params.Na,params.lamb,params.sz ,params.res,1,params.damp,displ);
[PSF,OTF] = GeneratePSF(params.Na,params.lamb,params.sz ,params.res,1,1,displ);
params.OTF=OTF;
params.OTF0=fftshift(double(OTF>0));
saveastiff(fftshift(PSF),psfname);

%% Declare some useful variables 
[params.I, params.J] = meshgrid(0:params.sz(2)-1,0:params.sz(1)-1); % Numerical mesh
params.X = params.I*params.res; params.Y=params.J*params.res;       % Scaled (by resolution) mesh 

if displ > 1
    subplot(2,3,2); imshow(wfs(:,:,1), []); title("WF"); colorbar; colormap(viridis);
    subplot(2,3,5); imshow(abs(log10(fftshift(fft2(wfs(:,:,1)))+1)), []); title("WF FT"); colorbar; colormap(viridis);
    subplot(2,3,3); imshow(y_out(:,:,1), []); title("Sample processed image"); colorbar; colormap(viridis);
    subplot(2,3,6); imshow(abs(log10(fftshift(fft2(y_out(:,:,1)))+1)), []); title("FT (sample processed image)"); colorbar; colormap(viridis);
end
