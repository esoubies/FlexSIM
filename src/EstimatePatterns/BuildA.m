function A = BuildA(ktest, wf, filt, params)
%--------------------------------------------------------------------------
% function A = BuildA(ktest, wf, filt, params)
%
% Builds the system A = [a1(:),a2(:)], where 
% a1 =   wf.* cos(2*(k_1*X+k_2*Y));            
% a2 = - wf.* sin(2*(k_1*X+k_2*Y)));
%
% Inputs :  ktest  → The wavevector to test
%           wf     → The corresponding wf information
%           filt   → The filter (in FT domain) to apply to A
%           params → Structure containing the parameters of the system
%
% Outputs : A      → Array of size [N, 2] (N being the # of pixels) that
%                    contains the sin and cos components of the pattern
%--------------------------------------------------------------------------

% Calculate necessary sins and cos first we get the grids corresponding to the test wavevector
k1x = 2*ktest(1)*params.X; k2y = 2*ktest(2)*params.Y;
% Now sins and cos of each (the ones correspondong to the X grid multiplied by the wf)
ck1x_wf = wf.*cos(k1x); ck2y = cos(k2y); sk1x_wf = wf.*sin(k1x); sk2y = sin(k2y);      
% Elementwise multiplications of the combinations needed for a1 & a2
c1c2 = ck1x_wf .* ck2y; s1c2 = sk1x_wf .* ck2y; s1s2 = sk1x_wf .* sk2y; c1s2 = ck1x_wf .* sk2y; 
if all(size(filt) == size(wf))
    % Get FTs
    c1c2_ft = fft2(c1c2); s1c2_ft = fft2(s1c2); s1s2_ft = fft2(s1s2); c1s2_ft = fft2(c1s2);
    a1=real(ifft2( filt.*(c1c2_ft - s1s2_ft)));
    a2=real(ifft2(-filt.*(s1c2_ft + c1s2_ft)));
else
    a1 =    c1c2 - s1s2;
    a2 = - (s1c2 + c1s2); 
end
A = [a1(:),a2(:)];

if params.method == 2
    % - Build system              
    for i = 2:params.nbPh
        deltaph = (i - 1)*pi/params.nbPh;
        c3 = cos(2*deltaph); s3 = sin(2*deltaph); 
        
        if all(size(filt) == size(wf))
            a1=real(ifft2( filt.*(c1c2_ft*c3 - s1s2_ft*c3 - s1c2_ft*s3 - c1s2_ft*s3)));
            a2=real(ifft2(-filt.*(s1c2_ft*c3 + c1s2_ft*c3 + c1c2_ft*s3 - s1s2_ft*s3)));
        else
            a1 =    c1c2*c3 - s1s2*c3 - s1c2*s3 - c1s2*s3;
            a2 = - (s1c2*c3 + c1s2*c3 + c1c2*s3 - s1s2*s3); 
        end
        A = vertcat(A, [a1(:),a2(:)]);
    end
end      
end