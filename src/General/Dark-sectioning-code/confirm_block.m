function block_size = confirm_block(params,lp)
PSF = PSF_Generator(params.emwavelength,params.pixelsize,params.NA,params.Nx,params.factor);
PSF_Lo = abs(ifft2(fftshift(fft2(PSF)).*fftshift(lp)));
PSF_Lo = PSF_Lo./max(max(PSF_Lo));
% figure;plot(PSF(floor(params.Nx/2)+1,:))
for count_x = floor(params.Nx/2):params.Nx
    if PSF_Lo(count_x,floor(params.Nx/2)) <0.01
        break;
    end
end
block_size = count_x-floor(params.Nx/2);
end