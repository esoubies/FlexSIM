function GenerateSIMData(imgType, params)
%--------------------------------------------------------------------------
% Function PeakDetect(x, MEP, params)
%
% Generates a stack of SIM data according to the specifications of params 
% (see below). Takes as input an image type, a noise level, and a structure
% containing the parameters of the (simulated) acquisition.
% Saves the SIM stack in the path specified by the parameter structure
% 
% Inputs :  imgType -> (int) Object to simulate. One of {0-> star, ...?}
%           MEP     -> Maximum expected number of photons 
%           params  -> Structures with fields:
%                         - DataPath (str) : Location to save the SIM stack in
%                         - nbPh (str)     : Number of phases per orientation
%                         - nbOr (str)     : Number of orientations
%                         - sz (list)      : The number of pixels in x and y of the image 
%                         - lamb (float)   : Emission wavelength [nm]
%                         - StackOrder(str): Phase (p), angle (a) and widefield (w)convention. 
%                         - Na             : Objective numerica aperture
%                         - res            : resolution of the SIM data stack
%                         - MEP            : Maximum estimated nb photons
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Generate OTF, Ground Truth
switch imgType
    case 0
        N = min(params.sz) - 1 * mod(min(params.sz), 2);             % GT Generation
        x = StarLikeSample(2,N,2,1,1,0); 
    case 1
        x = loadtiff('phantom.tif');
        x = imresize(x, 0.5);
        params.sz = size(x);
end
OTF = GenerateOTF(params.Na, params.lamb, params.sz, params.res, 1); % OTF Computation
% Pattern initialization according to convention
[X,Y]=meshgrid(0:params.sz(2)-1,0:params.sz(1)-1); X=X*params.res; Y=Y*params.res;
n = params.nbOr*params.nbPh+1*ismember('w',char(params.StackOrder));
patt = ones([params.sz n]);
for or = 1:params.nbOr                                             % Pattern generation
    k  = 2*pi*params.ns/params.lamb*[cos(params.or(or)), sin(params.or(or))]*params.Na/params.nl;
    for ph = 1:params.nbPh                                         % Iterate orientation and phase
        a = params.a;
        ac = a*cos(2*params.ph(ph)); as = a*sin(2*params.ph(ph));  % Calculate phase factors        
        patt(:,:,ph+params.nbPh*(or - 1)) = ...                    
            1 + ac*cos(2*(k(1)*X+k(2)*Y)) - as*sin(2*(k(1)*X+k(2)*Y));
    end
end

y=zeros(size(patt));
for ii=1:size(patt,3)
    y_noNoise = real(ifft2(OTF.*fft2(patt(:,:,ii).*x)));      % Simulate clean acquisitions
    y_noNoise = im2double(y_noNoise);                         % Min/max normalization
    scaling = 1e12/params.MEP;                                % Scaling factor (see imnoise doc)
    y(:,:,ii) = imnoise(y_noNoise/scaling,'poisson')*scaling; % Add poisson noise       
end

saveastiff(y, params.DataPath);
end