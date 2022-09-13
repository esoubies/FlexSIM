function [c, g] = evalJ(ktest, wf, G, params, results, filt, att_filt, grad)
% function [c, g] = evalJ(ktest, wf, G, eq_ph, filter, OTF, res, X, Y, nb_imgs)
%--------------------------------------------------------------------------
% [c, g] = evalJ(ktest, wf, G, filter, att_filter, params)
%
% Evaluate the cost of the wavevector k, with or without filters
% (attenuated for the actual image). Also calculates gradient on demand.
% 
% Input: ktest   : Wavevector 
%        wf      : widefield component
%        G       : SIM image(s) without widefield component
%        params  : Structure with main parameters
%        filter  : 0 (bool) or constant filter for the simulated image
%        att_filter  : 0 (bool) or attenuated filter for the acquired image
%        grad        : whether or not to calculate the gradient
%        displ   : Level of information to give to the user.
%
% Output: c  (float)     : Cost of the given wavevector
%         g  (1x2 float) : Gradient
%--------------------------------------------------------------------------
    nb_imgs = params.nbImgs;

    if all(size(att_filt) == size(wf))      % Already filter data             
        G=real(ifft2(fft2(G).*att_filt));
    end

%     if ~params.UseAllImages | (params.UseAllImages & params.eq_ph)
    if ismember(params.method, [0, 2])
        A = BuildA(ktest, wf, filt, params, results); 
        AA = A'*A;
        if params.method == 0
            G = G(:,:,1); 
        end
        s = AA\A'*G(:);
        c = 0.5*norm(A*s-G(:))^2/numel(G);
        if grad           
            g(1) = - (A*[s(2);-s(1)].*2.*repmat(results.X(:), nb_imgs, 1))'*(A*s-G(:))/numel(G);
            g(2) = - (A*[s(2);-s(1)].*2.*repmat(results.Y(:), nb_imgs, 1))'*(A*s-G(:))/numel(G);
        end 
    else

    c = 0;
    g = [0; 0];

    % - Build system    
    A = BuildA(ktest, wf, filt, params);
    AA = A'*A; 

    for ith_img = 1:nb_imgs
        % Select image (G if one image supplied, G_i in any other case)
        G_i = G(:,:,ith_img);
        
        % If grad was required, solve system, and calculate it.
        s = AA\A'*G_i(:);
        if grad                                    
            g(1) = g(1) - (A*[s(2);-s(1)].*2.*results.X(:))'*(A*s-G_i(:))/numel(G_i);
            g(2) = g(2) - (A*[s(2);-s(1)].*2.*results.Y(:))'*(A*s-G_i(:))/numel(G_i);
            g = g';
        end 
        c = c + 0.5*norm(A*s-G_i(:))^2/numel(G_i);
    end    
    end
end