function A = BuildA(ktest, wf, filt, params, grids)
%--------------------------------------------------------------------------
% function A = BuildA(ktest, wf, filt, params, grids)
%
% Builds the system A = [a1(:),a2(:)], where 
% a1 =   wf.* cos(2*(k_1*X+k_2*Y));            
% a2 = - wf.* sin(2*(k_1*X+k_2*Y)));
% To avoid the repetition of Cosine, Sine and FFT calculations we take
% advantage of the trigonometric identities of sum of angles, as well as of
% the linearity of the Fourier transform, to only calculate expensive
% cosines and FFTs of each component once. 
%
% Inputs :  ktest   -> The wavevector to test 
%           wf      -> The corresponding wf information
%           filt    -> The filter (in FT domain) to apply to A. If set to 0,
%                      no filter is applied           
%           params  -> Structure containing the input parameters of the system
%           grids   -> Structure that stores all necessary grids
%
% Outputs : A       -> Array of size [N, 2] (N being the # of pixels) that
%                      contains the sin and cos components of the pattern
%
% [1] FlexSIM: ADD REF TO PAPER
% 
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com), 
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Calculate necessary sins and cos first we get the grids corresponding to the test wavevector
k1x = 2*ktest(1)*grids.X; k2y = 2*ktest(2)*grids.Y;
% Now sins and cos of each (the ones correspondong to the X grid multiplied by the wf)
ck1x_wf = wf.*cos(k1x); ck2y = cos(k2y); sk1x_wf = wf.*sin(k1x); sk2y = sin(k2y);      
% Elementwise multiplications of the combinations needed for a1 & a2
c1c2 = ck1x_wf .* ck2y; s1c2 = sk1x_wf .* ck2y; s1s2 = sk1x_wf .* sk2y; c1s2 = ck1x_wf .* sk2y; 

if all(size(filt) == size(wf))
    % Calculate FFTs if necessary
    c1c2_ft = fft2(c1c2); s1c2_ft = fft2(s1c2); s1s2_ft = fft2(s1s2); c1s2_ft = fft2(c1s2);
    % Build system
    a1=real(ifft2( filt.*(c1c2_ft - s1s2_ft)));
    a2=real(ifft2(-filt.*(s1c2_ft + c1s2_ft)));
else
    % Build system
    a1 =    c1c2 - s1s2;
    a2 = - (s1c2 + c1s2); 
end

% If there is the assumption of equally spaced phases, build according system
if params.eqPh
    % Preallocate A for performance
    if params.GPU
        A= zeros([params.nbPh*numel(grids.X(:)), 2],'double','gpuArray');
    else
        A = zeros([params.nbPh*numel(grids.X(:)), 2]);
    end
    
    A(1:numel(grids.X(:)), :) = [a1(:),a2(:)];
    for i = 2:params.nbPh
        deltaph = (i - 1)*pi/params.nbPh;         % Calculate deltaPhi
        c3 = cos(2*deltaph); s3 = sin(2*deltaph); % Calculate additional sines and cos        
        if all(size(filt) == size(wf))            % Calculate additional FFTs and build
            a1=real(ifft2( filt.*(c1c2_ft*c3 - s1s2_ft*c3 - s1c2_ft*s3 - c1s2_ft*s3)));
            a2=real(ifft2(-filt.*(s1c2_ft*c3 + c1s2_ft*c3 + c1c2_ft*s3 - s1s2_ft*s3)));
        else
            a1 =    c1c2*c3 - s1s2*c3 - s1c2*s3 - c1s2*s3;
            a2 = - (s1c2*c3 + c1s2*c3 + c1c2*s3 - s1s2*s3); 
        end
        A((i-1)*numel(grids.X(:))+1:i*numel(grids.X(:)),:) = [a1(:),a2(:)];
    end
else
        A = [a1(:),a2(:)];                          % Simply declare A
end  
end