function [ radiance ] = dehaze_fast2(  image, omega, win_size,EL,dep,thres )
% Copyright (c) 2014 Stephen Tierney


[Nx,Ny] = size(image);
if ~exist('omega', 'var')
    omega = 0.95;
end

if ~exist('win_size', 'var')
    win_size = 15;
end

r = 15;
res = 0.001;

[m, n, ~] = size(image);

Mask = zeros(Nx,Ny);
Mask(image<thres) = 1;
dark_channel = get_dark_channel(image.*Mask, win_size);
min_atmosphere = get_atmosphere(image.*Mask, dark_channel);

dark_channel = get_dark_channel(image, win_size);
max_atmosphere = get_atmosphere(image, dark_channel);
EL = EL - min(min(EL));
rep_atmosphere_process = EL/max(max(EL))*(max_atmosphere-min_atmosphere)+min_atmosphere;
rep_atmosphere_process = dep*rep_atmosphere_process;
trans_est = get_transmission_estimate(rep_atmosphere_process,image, omega, win_size);
x = guided_filter(image, trans_est, r, res);
transmission = reshape(x, m, n);
radiance = get_radiance(rep_atmosphere_process,image, transmission);


end





