function OTF = OTF_estimate(fft_image)
[Nx,Ny] = size(fft_image);

 
mx = floor(Ny/2)+1;
my = floor(Ny/2)+1;

r = zeros(Nx*Ny,1);
k = zeros(Nx*Ny,1);
for jx = 1:Nx
    for jy = 1:Ny
        r((jx-1)*Nx+jy,1) = sqrt((jx-mx).^2 + (jy-my).^2);
        k((jx-1)*Nx+jy,1) = fft_image(jx,jy);
    end
end



end