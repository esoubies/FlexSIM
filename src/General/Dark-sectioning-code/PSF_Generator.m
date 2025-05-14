function PSF = PSF_Generator(lambada,pixelsize,NA,w,factor)

[X,Y]=meshgrid(linspace(0,w-1,w),linspace(0,w-1,w));
scale=2*pi*NA/lambada*pixelsize;
scale=scale*factor;

R=sqrt(min(X,abs(X-w)).^2+min(Y,abs(Y-w)).^2);
PSF=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2;
PSF = PSF/(sum(sum(PSF)));
PSF=fftshift(PSF);

end

