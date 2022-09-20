function [ph,a]=GetPhaseAndAmp(k,wf,G,grids,OTF,sz,params)
%--------------------------------------------------------------------------
% function [ph,a]=GetPhaseAndAmp(k,wf,G,grids,OTF,sz,params)
%
% Given a wavevector k, compute the optimal phase and amplitude through the
% minimization of J (solving of a linear system). See the function EvalJ for
% input parameters.
%
% See also EvalJ.m
% 
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Build matrix A
[att_filt, filt] = BuildFilter(k, sz, OTF, params, grids);
A=BuildA(k, wf, filt, params, grids); AA = A'*A;
% -- Build second term b
G_filt=real(ifft2(fft2(G).*att_filt));
% -- Solve linear system
if params.method == 2
    s = AA\A'*G_filt(:);
else
    G_filt_tmp= G_filt(:,:,1);
    s = AA\A'*G_filt_tmp(:);
end
% -- Extract cos and sin components
ac=s(1); as=s(2);

% -- Convert (cos,sin) --> angle
tmp = atan(as/ac);           % Calculate phase, ac, and as
if ac <0, tmp=pi+tmp; end
ph=mod(tmp,2*pi)/2;   

% - Get amplitude of the pattern
% TODO : To improve...
A=BuildA(k, wf, 0, params, grids);
AA = A'*A;                                   % Recalculate ac and as without...
if params.method == 2                        % filters and extract amplitude
    s = AA\A'*G(:);
else
    G_tmp= G(:,:,1); s = AA\A'*G_tmp(:);
end
a=sqrt(s(1)^2 + s(2)^2);

end