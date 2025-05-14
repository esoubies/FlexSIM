function radiance = get_radiance(rep_atmosphere,image, transmission)
% Copyright (c) 2014 Stephen Tierney

[m, n, ~] = size(image);

max_transmission = max(transmission, 0.1);

radiance = ((image - rep_atmosphere) ./ max_transmission) + rep_atmosphere;

end