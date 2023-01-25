function [c, g] = EvalJ(ktest, wf, G, params, grids, filt, att_filt, grad)
%--------------------------------------------------------------------------
% function [c, g] = EvalJ(ktest, wf, G, params, grids, filt, att_filt, grad)
%
% Evaluates the cost `|As - G|^2` with `A` and `s` being determined by the 
% wavevector `k`. Depending on the parameters, A and G can be filtered. If 
% requested, also calculates the gradient.
% 
% Inputs : ktest           -> Wavevector 
%          wf              -> widefield component
%          G               -> SIM image(s) without widefield component
%          params          -> Structure with main parameters
%          grids           -> Structure with the grids of system
%          filter          -> 0 (bool, filter not used) or Fourier filter for `A`
%          att_filter      -> 0 (bool) or Fourier filter for the acquired image
%          grad            -> whether or not to calculate the gradient
%
% Output : c  (float)      -> Cost of the given wavevector
%          g  (1x2 float)  -> Gradient
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

if all(size(att_filt) == size(wf)) && (~isfield(params,'useFilter') || params.useFilter)   % If a filter was given, filter data
    G=real(ifft2(fft2(G).*att_filt));
else
    G=real(ifft2(fft2(G).*double(att_filt>0)));
end

if params.eqPh                                   % If no summation needed...
    A = BuildA(ktest, wf, filt, params, grids);  % Build system (inner function takes
    AA = A'*A;                                   % care of the size and delta ph)
    s = AA\A'*G(:);
    c = 0.5*norm(A*s-G(:))^2/numel(G);           % Extract cost and gradient
%     if grad
%         tmp=(A*s-G(:))/numel(G);
%         g(1) = - (A*[s(2);-s(1)].*2.*repmat(grids.X(:), params.nbPh, 1))'*tmp;
%         g(2) = - (A*[s(2);-s(1)].*2.*repmat(grids.Y(:), params.nbPh, 1))'*tmp;
%     end
    if grad
        A_noFilt = BuildA(ktest, wf, 0, params, grids);
        tmp=real(ifft2(conj(filt).*fft2(reshape((A*s-G(:))/numel(G),size(G)))));
        g(1) = - (A_noFilt*[s(2);-s(1)].*2.*repmat(grids.X(:), params.nbPh, 1))'*tmp(:);
        g(2) = - (A_noFilt*[s(2);-s(1)].*2.*repmat(grids.Y(:), params.nbPh, 1))'*tmp(:);
    end
else
    c = 0;                               % Initialize variables and system `A`
    g = [0; 0];
    A = BuildA(ktest, wf, filt, params, grids);
    AA = A'*A;
    for phNb = 1:params.nbPh             
        G_i = G(:,:,phNb);               % Select corresponding image
        s = AA\A'*G_i(:);                % Extract phase of current image
        c = c + 0.5*norm(A*s-G_i(:))^2/numel(G);
%         if grad
%             tmp=(A*s-G_i(:))/numel(G);
%             g(1) = g(1) - (A*[s(2);-s(1)].*2.*grids.X(:))'*tmp(:);
%             g(2) = g(2) - (A*[s(2);-s(1)].*2.*grids.Y(:))'*tmp(:);
%             g = g';
%         end
        if grad
            A_noFilt = BuildA(ktest, wf, 0, params, grids);
            tmp=real(ifft2(conj(filt).*fft2(reshape((A*s-G_i(:))/numel(G),size(G_i)))));
            g(1) = g(1) - (A_noFilt*[s(2);-s(1)].*2.*grids.X(:))'*tmp(:);
            g(2) = g(2) - (A_noFilt*[s(2);-s(1)].*2.*grids.Y(:))'*tmp(:);
            g = g';
        end
    end
end
end