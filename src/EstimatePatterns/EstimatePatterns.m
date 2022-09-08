function [k_final, phase, a, params] = EstimatePatterns(y, wfs, params, displ)
%--------------------------------------------------------------------------
% function [k_final, phase, params.a, params] = EstimateParameters(y, wfs, params, displ)
%
% Inputs : y, 
%
% 
%--------------------------------------------------------------------------

% - Pre computations of J evaluation
cmin=params.limits(1);        % Ring used for peaks search
cmax=params.limits(2);
fc_n=params.fc*pi/params.res;
maxp=params.fc*pi/params.res*max(1,cmax); % Maximum value of one component of k
cen = [0,0];
ll1=linspace(cen(1),cen(1)+maxp,params.n_pts/2);           % For half grid
ll2=linspace(cen(2)-maxp,cen(2)+maxp,params.n_pts);
[K1,K2]=meshgrid(ll1,ll2);                          % Grid of wavevectors
K = cat(3, K1, K2); 
Knorm = vecnorm(K, 2, 3);
directionCounter = 1;

% Loop through images to calculate J
idwf=0;
for idx = params.imgIdxs
    disp([' Batch of images: ',num2str(idx')]);
    % Get the corresponding wf for the current batch of images
    idwf=idwf+1;
    wf =  wfs(:,:,idwf);
    disp('   - Grid-based evaluation of J landscape...');
    Jp = zeros(params.n_pts, params.n_pts/2);
    % nb_imgs = length(idx);
    G = y(:,:,idx);
    if displ > 1
        figure; subplot(1, 2, 2); imshow(G(:,:,1), []); title("Filtered Image"); colorbar; colormap(fire);
        subplot(1, 2, 2); imshow(abs(log10(fftshift(fft2(G(:,:,1)))+1)), []); title("Filtered Image FT"); colorbar; colormap(fire);
    end

    for ii=1:params.n_pts
        for jj=1:params.n_pts/2
            if Knorm(ii, jj)>cmin*fc_n && Knorm(ii, jj)<cmax*fc_n
                ktest = squeeze(K(ii, jj, :));
                % Evaluate cost without filtering 
                Jp(ii, jj) = EvalJ(ktest, wf, G, params, 0, 0, 0);
            else
                Jp(ii,jj)=NaN;
            end
        end
    end

    % Estimate the minima of J   
    if ~isfield(params, "n_minima")
        params.n_minima = 1;
    end
    disp(['   - Extracting the ',num2str(params.n_minima),' smallest local minima...']);
    k_est_landscape = zeros(params.n_minima, 2);
    Jp_tmp = Jp; 
    for nth = 1:params.n_minima
        [~,idxMin]=min(Jp_tmp(:));                     % Extract minima
        [ii,jj]=ind2sub([params.n_pts,params.n_pts/2],idxMin);
        k_est_tmp=[K1(ii,jj);K2(ii,jj)];               
%         if k_est_tmp(1) == 0                           % Is this necessary? We're using half anyway
%             k_est_tmp(1) = 1e-10;
%         end
%         k_est_tmp=k_est_tmp*sign(k_est_tmp(1)); 
        k_est_landscape(nth, :) = k_est_tmp;
        % Erase temporarily the vecinity of the found minima
        idx1 = ii - round(params.n_pts/50); idx1(idx1<1) = 1;
        idx2 = ii + round(params.n_pts/50); idx2(idx2<params.n_pts) = params.n_pts;
        idx3 = jj - round(params.n_pts/50); idx3(idx3<1) = 1;
        idx4 = jj + round(params.n_pts/50); idx4(idx4>params.n_pts/2) = params.n_pts/2;
        Jp_tmp(idx1:idx2, idx3:idx4)  = max(Jp(:));
    end
    
    if displ
        mJp=max(Jp(:));
        fg=figure; subplot(1,2,1); axis xy;hold on; view(45,45);
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14); axis([0 maxp -maxp maxp]);grid;
        for nth = 1:params.n_minima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20); %#ok<AGROW> 
           % line([k_est_landscape(nth, 1) k_est_landscape(nth, 1)],[k_est_landscape(nth, 2) k_est_landscape(nth, 2)],[min(Jp(:)),max(Jp(:))],'Color', color{nth},'Linewidth',2);
        end        
        title('J landscape and initial local minima');
        subplot(1,2,2); axis xy;hold on; title('Zoom');
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14)
        for nth = 1:params.n_minima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20);
           % line([k_est_landscape(nth, 1) k_est_landscape(nth, 1)],[k_est_landscape(nth, 2) k_est_landscape(nth, 2)],[min(Jp(:)),max(Jp(:))],'Color', color{nth},'Linewidth',2);
        end  
        drawnow;
    end

    % Refine each of the J with an iterative approach    
    disp('   - Refine position of extracted local min...');
    fprintf('%s','     - local min #');
    for ithk = 1:params.n_minima
        if mod(ithk,ceil(params.n_minima/10))==0
            if ithk==params.n_minima
                fprintf('%i\n',ithk);
            else
                fprintf('%i, ',ithk);
            end
        end

        G = y(:,:,idx);                       % Select image and minima
        k=k_est_landscape(ithk, :); k = k(:);
        if displ>1
            fig=figure;
        end
        count_it = 1; 
        clear cf
        if ~isfield(params, "FilterRefinement")
            params.FilterRefinement = 1;
        end

        for nit_filt = 1:params.FilterRefinement
            [att_filt, filt] = BuildFilter(k, params, displ);
    
            % Gradient descent parameters
            tau=1e-2;
            tau_fact=1.5;
            tau_min=1e-10;            
            nit_tot=100;
            cf_tolerance = 1e-6;
            %grad_tol=1e-12;
            
            cf(count_it) = EvalJ(k, wf, G, params, filt, att_filt, 0); %#ok<AGROW> 
            
            for jj=1:nit_tot
            % -- Update Wavevector with gradient descent steps
                [~, g] = EvalJ(k, wf, G, params, filt, att_filt, 1);
                ktmp = k - tau * g';
                cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                % Line search by Backtracking
                if cf(count_it) > cf_new
                    while cf(count_it) > cf_new
                        k = ktmp;
                        tau = tau * tau_fact;
                        ktmp = k - tau * g';
                        cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                    end
                    % Reverse to previous ktmp (as new cf is higher)
                    ktmp = k; 
                    tau = tau/tau_fact;
                else
                    while cf(count_it)<=cf_new && (tau>tau_min)
                        tau=tau/tau_fact;
                        ktmp = k -tau*g';
                        cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                    end
                    if (tau>tau_min)
                        k=ktmp;
                    end
                end
                if (tau<tau_min) % Stop if the descent step becomes too small
                    break;
                end
                count_it=count_it+1;
                cf(count_it) = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);        
            if displ>1
                figure(fig);
                plot(cf,'linewidth',2); xlim([0 nit_tot*params.FilterRefinement]);grid;
                title('Freq/Phase Opti: Cv curve')
                set(gca,'FontSIze',14);
                drawnow;
            end
            if (tau<tau_min) || abs(cf(count_it)-cf(count_it-1))/abs(cf(count_it-1))<cf_tolerance %
                break;
            end        
            end

        end
        % Store final cost and wavevector
        k_alt_est(ithk, :) = k; %#ok<AGROW> 
        [att_filt, filt] = BuildFilter(k, params, displ);
        A=BuildA(k, wf, filt, params); 
        AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));       
        if params.method == 2
            s = AA\A'*G_filt(:);
        else
            G_filt_tmp= G_filt(:,:,1); 
            s = AA\A'*G_filt_tmp(:);
        end           
        ac(ithk)=s(1); as(ithk)=s(2);        %#ok<AGROW> 

        
    end

    disp('   - Choosing the best...'); % By reevaluating without filters
    lowestJ = realmax;
    for iii = 1:params.n_minima        
        Jp = EvalJ(k_alt_est(iii,:), wf, G, params, 0, 0, 0);
        if Jp < lowestJ  
            lowestJ = Jp;
            optIdx = iii;
        end
    end

    if displ
        figure(fg);subplot(1,2,1);
        plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        subplot(1,2,2);
        plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        axis([k_alt_est(optIdx, 1) - maxp*0.1 k_alt_est(optIdx, 1) + maxp*0.1 k_alt_est(optIdx, 2) - maxp*0.1 k_alt_est(optIdx, 2) + maxp*0.1 ]);grid;
        pause(0.1);
    end

    k_final(directionCounter, :) = k_alt_est(optIdx, :); %#ok<AGROW> 

    % Calculate phase and stick to convention
    tmp = atan(as(optIdx)/ac(optIdx));
    if ac(optIdx) <0, tmp=pi+tmp; end
    ph_init=mod(tmp,2*pi)/2; 
    phase(directionCounter, 1) = ph_init; %#ok<AGROW> 
    params.phase(directionCounter, 1) = ph_init;
    params.ac(directionCounter,1) = ac(optIdx); 
    params.as(directionCounter,1) = as(optIdx);
    
    % Extract amplitude factor (no filter used)
    A=BuildA(k_final(directionCounter, :), wf, filt, params); AA = A'*A;           
    if params.method == 2
        s = AA\A'*G(:);
    else
        G_tmp= G(:,:,1); s = AA\A'*G_tmp(:);
    end           
    ac(ithk)=s(1); as(ithk)=s(2); params.a(directionCounter, 1) = as(ithk)./sin(ph_init); 

    if ~isfield(params, "store_all_ph")
        params.store_all_ph = 0; 
    end
    if params.store_all_ph && params.method < 2
        k = k_alt_est(optIdx, :); 
        [att_filt, filt] = BuildFilter(k, params, displ);
        A=BuildA(k, wf, filt, params); AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));

        % For amplitude estimation
        ANoFilt = BuildA(k_final(directionCounter, :), wf, filt, params); 
        AANoFilt = ANoFilt'*ANoFilt;

        for i = 1:params.nbPh - 1
            Gtmp = G_filt(:,:,i + 1); 
            s = AA\A'*Gtmp(:);
            ac_tmp=s(1); as_tmp=s(2);

            % For amplitude estimation
            GtmpNoFilt = G(:,:,i + 1); 
            sNoFilt = AANoFilt\ANoFilt'*GtmpNoFilt(:);
            acNoFilt=sNoFilt(1); asNoFilt=sNoFilt(2); %#ok<NASGU> 
            
            tmp=atan(as_tmp/ac_tmp);
            if ac_tmp < 0, tmp=pi+tmp; end
            tmp=mod(tmp,2*pi)/2;
            phase(directionCounter, i+1) = tmp;
            params.phase(directionCounter, i+1) = tmp;
            params.ac(directionCounter,i+1) = ac_tmp; 
            params.as(directionCounter,i+1) = as_tmp;
            params.a(directionCounter, 1) = asNoFilt/sin(tmp);
        end
    end
%     fprintf("Parameters found for orientation %d:\n Kx = %5.4f | Ky = %5.4f | Ph = %5.4f\nk", ...
%         directionCounter, k_alt_est(optIdx, 1), k_alt_est(optIdx, 2), phase)
    directionCounter = directionCounter + 1;
%     display(["Parameters found for orientation", num2])
end
if displ
    pause(0.1);
    figure; imshow(log10(sum(abs(fftshift(fft2(y))),3)+1), []); colormap(viridis); 
    colorbar; axis on; hold on;
    % Plot cross at row 100, column 50
    for i = 1:params.nbOr
        tmp = k_final(i, :) * params.sz(1) * params.res / pi + params.sz(1)/2+1; 
        plot(tmp(1), tmp(2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    end
    title('Sum of data Fourier amplitudes + detected wavevectors'); drawnow;
end
% After parameter estimation, return to original size
params.sz=params.origSz; params.k = k_final; a = params.a;
disp('=== Patterns parameter estimation END');
end