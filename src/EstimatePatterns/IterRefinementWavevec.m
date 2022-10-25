function k = IterRefinementWavevec(k_init,wf,G,grids,OTF,sz,params)
%--------------------------------------------------------------------------
% function k = IterRefinementWavevec(k_init,wf,G,grids,OTF,sz,params)
%
% Given an initial wavevector k_init, minimize locally the function J
% through gradient descent with line search. See the function EvalJ for
% input parameters.
%
% See also EvalJ.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% -- Parameters
tau_fact=1.5;        % Factor to update the descent step-size during line search
tau_min=1e-10;       % Minimal value of step-size
nit_tot=100;         % Max number of iterations
cf_tolerance = 1e-6; % Tolerance on the relative difference of the cost function btw 2 iterates
tau_init=1e-2;       % Initial descent step-size

% -- Initializations
if params.displ>2
    fig=figure;
end
k=k_init;
count_it = 1;  

% -- Outer loop (each time we update the filters)
for nit_filt = 1:params.FilterRefinement   % Iterate filter refinement
    % Update filters and (re)initialization of step-size
    [att_filt, filt] = BuildFilter(k, sz, OTF, params, grids);                       
    tau=tau_init; 
    
    % Calculate current cost
    cf(count_it) = EvalJ(k, wf, G, params, grids, filt, att_filt, 0); %#ok<AGROW>
    
    % Gradient descent w.r.t. wavevector
    for jj=1:nit_tot                      
        [~, g] = EvalJ(k, wf, G, params, grids, filt, att_filt, 1);    % Get current gradient
        ktmp = k - tau * g';                                           % Update k
        cf_new = EvalJ(ktmp, wf, G, params, grids, filt, att_filt, 0); % Calculate new cost
        
        % Line search by backtracking
        if cf(count_it) > cf_new
            while cf(count_it) > cf_new
                k = ktmp;
                tau = tau * tau_fact;
                ktmp = k - tau * g';
                cf_new = EvalJ(ktmp, wf, G, params, grids, filt, att_filt, 0);
            end
            ktmp = k; tau = tau/tau_fact;
        else
            while cf(count_it)<=cf_new && (tau>tau_min)
                tau = tau/tau_fact;
                ktmp = k - tau * g';
                cf_new = EvalJ(ktmp, wf, G, params, grids, filt, att_filt, 0);
            end
            if (tau>tau_min)
                k=ktmp;
            end
        end
        
        % Stop grad descent if step size becomes negligible
        if (tau<tau_min) 
            break;
        end
        
        % Calculate and store new cost
        count_it=count_it+1;                  
        cf(count_it) = EvalJ(ktmp, wf, G, params, grids, filt, att_filt, 0);
        if params.displ>2
            figure(fig);
            plot(cf,'linewidth',2); xlim([0 nit_tot]);grid;
            title('Freq/Phase Opti: Cv curve')
            set(gca,'FontSIze',14);
            drawnow;
        end
        
        % Stop grad descent if cost improvement becomes negligible
        if abs(cf(count_it)-cf(count_it-1))/abs(cf(count_it-1))<cf_tolerance
            break;
        end

    end
    
end