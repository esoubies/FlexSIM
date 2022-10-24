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
nb_imgs = size(G,3)  ;

if all(size(att_filt) == size(wf))      % If a filter was given, filter data             
    G=real(ifft2(fft2(G).*att_filt));
end

if ismember(params.method, [0, 2])               % If no summation needed...
    A = BuildA(ktest, wf, filt, params, grids);  % Build system (inner function takes 
    AA = A'*A;                                   % care of the size and delta ph)
    if params.method == 0
        G = G(:,:,1);                            % If we're only using one image select it
    end
    s = AA\A'*G(:);
    c = 0.5*norm(A*s-G(:))^2/numel(G);           % Extract cost and gradient
    if grad                                        
        g(1) = - (A*[s(2);-s(1)].*2.*repmat(grids.X(:), nb_imgs, 1))'*(A*s-G(:))/numel(G);
        g(2) = - (A*[s(2);-s(1)].*2.*repmat(grids.Y(:), nb_imgs, 1))'*(A*s-G(:))/numel(G);
    end 
else
    c = 0;                               % Initialize variables and system `A`
    g = [0; 0];   
    A = BuildA(ktest, wf, filt, params, grids);
    AA = A'*A; 
    
    for ith_img = 1:nb_imgs              % Iterate the images for summation        
        G_i = G(:,:,ith_img);            % Select corresponding image               
        s = AA\A'*G_i(:);                % Extract phase of current image
        if grad                                    
            g(1) = g(1) - (A*[s(2);-s(1)].*2.*grids.X(:))'*(A*s-G_i(:))/numel(G_i);
            g(2) = g(2) - (A*[s(2);-s(1)].*2.*grids.Y(:))'*(A*s-G_i(:))/numel(G_i);
            g = g';
        end 
        c = c + 0.5*norm(A*s-G_i(:))^2/numel(G_i);
    end    
end
end