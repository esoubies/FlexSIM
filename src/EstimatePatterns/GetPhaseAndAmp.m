function ph=GetPhaseAndAmp(k,wf,G,grids,OTF,sz,params)
%--------------------------------------------------------------------------
% function ph=GetPhaseAndAmp(k,wf,G,grids,OTF,sz,params)
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
[OTFshiftCrop, OTFCrop] = BuildFilter(k, sz, OTF, params, grids);
A=BuildA(k, wf, OTFCrop, params, grids); AA = A'*A;
% -- Build second term b
G_filt=real(ifft2(fft2(G).*OTFshiftCrop));
% -- Solve linear system
if params.method == 2
    s = AA\A'*G_filt(:);
elseif params.method == 1
    s = zeros(2, params.nbPh);
    for phNb = 1:params.nbPh
        G_filt_tmp= G_filt(:,:,phNb);
        s(:,phNb) = AA\A'*G_filt_tmp(:);
    end
else
    G_filt_tmp= G_filt;
    s = AA\A'*G_filt_tmp(:);
end
% -- Extract cos and sin components
ac=s(1,:); as=s(2,:);

% -- Convert (cos,sin) --> angle
tmp = atan(as./ac);           % Calculate phase, ac, and as
tmp = tmp + pi * (ac < 0);
ph=mod(tmp,2*pi)/2;   

% - Get amplitude of the pattern
% Not used : a is now a parameter that needs to be adjusted by the user
% A=BuildA(k, wf, 0, params, grids);
% AA = A'*A;                                   % Recalculate ac and as without...
% if params.method == 2                        % filters and extract amplitude
%     s = AA\A'*G(:);
% elseif params.method == 1
%     s = zeros(2, 3);
%     for phNb = 1:params.nbPh
%         G_filt_tmp= G_filt(:,:,phNb);
%         s(:,phNb) = AA\A'*G_filt_tmp(:);
%     end
% else
%     G_filt_tmp= G_filt;
%     s = AA\A'*G_filt_tmp(:);
% end
% a=sqrt(s(1,:).^2 + s(2,:).^2);

end