function trans_est = get_transmission_estimate(rep_atmosphere, image, omega, win_size)
% Copyright (c) 2014 Stephen Tierney
[m, n, ~] = size(image);

%rep_atmosphere = repmat(atmosphere, m, n);

trans_est = 1 - omega * get_dark_channel( image ./ rep_atmosphere, win_size);

end