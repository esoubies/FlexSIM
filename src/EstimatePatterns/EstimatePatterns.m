function [k_final, phase, a, params] = EstimatePatterns(params, displ)
%--------------------------------------------------------------------------
% Function [k_final, phase, a, params] = EstimatePatterns(params)
%
% EstimatePatterns takes as input the structure defined in Main.m (now
% containing the SIM data) and returns the parameters (wavevector, phase
% and amplitude) of the SIM patterns
%
% Inputs  : params  → Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, SIM image, etc.)  
%
% Outputs : k_final → Structure with the final image (in the field `res`)
%                    and other intermediate results, like patterns. 
%           phase   → Depending on the method, of size (nbOr, nbPh) or just
%                     (nbPh). Contains the phase offset of the images 
%           a       → Same size as phase, contains the amplitude of the SIM patterns
%           params  → Same input structure, now with the pattern parameter info
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com) , E. Soubies 
% (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Preprocessing of stack 

% If requested, we crop to the ROI, keeping the original size before
if isfield(params,'roi') && ~isempty(params.roi)   % Crop to ROI
    y = params.y(params.roi(1):params.roi(1)+params.roi(3)-1,params.roi(2):params.roi(2)+params.roi(3)-1,:);
end
sz = size(y); % Recalc size if it was modified

% Computation of the OTF
params.OTF = GenerateOTF(params.Na, params.lamb, sz ,params.res, 1);
params.OTF0=fftshift(double(params.OTF>0));
% saveastiff(fftshift(PSF),psfname); Is it necessary?

% Calculate and substract WF, mask with rings on ROIs
% yTmp = zeros(params.sz); ERASE - deemed unnecessary                     
wfs =zeros(sz(1), sz(2),params.nbOr); % Initialize variables
y = (y - min(y(:))) / (max(y(:)) - min(y(:))); % Normalize images

% Select the indexes to use with array imgIdxs. If Zeiss convention (zap), flip
imgIdxs = 1:params.nbPh*params.nbOr; 
imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]); 
if strcmp(params.AcqConv, "zap")  
    imgIdxs = imgIdxs'; 
end
% If there are any assumptions on multiple images, we will use all the phases
if params.method > 0
    params.nbImgs = params.nbPh; 
else 
    % Else we will only use the first image of each orientation
    imgIdxs = imgIdxs(1,:); 
    params.nbImgs = 1;
end

% Remove widefield and mask
idwf=0;
FCut = 2*params.Na/params.lamb*params.res;
for idx = imgIdxs
    idwf=idwf+1;
    % Compute the wf and fft per batch of images (orientation)
    wfs(:,:,idwf)=sum(y(:,:,idx), 3)/length(idx); fft_wf = fftshift(fft2(wfs(:,:,idwf)));
    % Treat images individually - get FFT, remove scaled fftWF & mask
    for img = idx' 
        ffty=fftshift(fft2(y(:,:,img)));   % image and FFT data
        a=real(OptWght(ffty,fft_wf,(sz(1:2)+1)/2, sz(1)/8)); % FCut * 0.2   
        % Remove the scaled WF + filter freq within a ring of interest
        y(:,:,img) = real(ifft2(ifftshift(MaskFT((a*ffty-fft_wf), FCut, params.ringMaskLim)))); 
    end
end

% Mask widefield
for ii=1:size(wfs,3)
    wfs(:,:,ii) = real(ifft2(ifftshift(MaskFT(fftshift(fft2(wfs(:,:,ii))), FCut, params.ringMaskLim))));
end

%% Precomputations for the evaluation of the cost function J 
[params.I, params.J] = meshgrid(0:sz(2)-1,0:sz(1)-1); % Numerical mesh
params.X = params.I*params.res; params.Y=params.J*params.res;       % Scaled (by resolution) mesh 

cmin=params.limits(1);  cmax=params.limits(2); % Ring used for peaks search
FCutN=FCut*pi/params.res;             % Scaled cutoff frequency
maxp=FCut*pi/params.res*max(1,cmax);  % Maximum value of one component of k
cen = [0,0];                          % Grid center
ll1=linspace(cen(1),cen(1)+maxp,params.nPoints/2);    % Declare grid of 
ll2=linspace(cen(2)-maxp,cen(2)+maxp,params.nPoints); % wavevectors with the 
[K1,K2]=meshgrid(ll1,ll2);                          % x direction halved
K = cat(3, K1, K2); 
Knorm = vecnorm(K, 2, 3);                           % Calculate norms

OrientCount = 1;                                    % Dummy var to count orientations
% Loop through images to calculate J
for idx = imgIdxs
    disp([' Batch of images: ', num2str(idx')]);   
    wf =  wfs(:,:,OrientCount);                  % Get corresponding wf 
    disp('   - Grid-based evaluation of J landscape...');
    Jp = zeros(params.nPoints, params.nPoints/2);    % Initialize Jp for batch
    G = y(:,:,idx);                              % Select current batch    

    for ii=1:params.nPoints
        for jj=1:params.nPoints/2
            if Knorm(ii, jj)>cmin*FCutN && Knorm(ii, jj)<cmax*FCutN
                ktest = squeeze(K(ii, jj, :));
                % Evaluate cost without filtering 
                Jp(ii, jj) = EvalJ(ktest, wf, G, params, 0, 0, 0);
            else
                Jp(ii,jj)=NaN;
            end
        end
    end
    
    disp(['   - Extracting the ',num2str(params.nMinima),' smallest local minima...']);
    k_est_landscape = zeros(params.nMinima, 2);
    Jp_tmp = Jp;

    % Analyze Jp (through Jp_tmp) to extract the n local minima
    for nth = 1:params.nMinima
        [~,idxMin]=min(Jp_tmp(:));                     % Extract minima
        [ii,jj]=ind2sub([params.nPoints,params.nPoints/2],idxMin);
        k_est_tmp=[K1(ii,jj);K2(ii,jj)];               
        k_est_landscape(nth, :) = k_est_tmp;

        % Erase temporarily the vecinity of the found minima
        idx1 = ii - round(params.nPoints/50); idx1(idx1<1) = 1;
        idx2 = ii + round(params.nPoints/50); idx2(idx2<params.nPoints) = params.nPoints;
        idx3 = jj - round(params.nPoints/50); idx3(idx3<1) = 1;
        idx4 = jj + round(params.nPoints/50); idx4(idx4>params.nPoints/2) = params.nPoints/2;
        Jp_tmp(idx1:idx2, idx3:idx4)  = max(Jp(:));
    end
    
    if displ
        mJp=max(Jp(:));
        fg=figure; subplot(1,2,1); axis xy;hold on; view(45,45);
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14); axis([0 maxp -maxp maxp]);grid;
        for nth = 1:params.nMinima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20); %#ok<AGROW> 
           % line([k_est_landscape(nth, 1) k_est_landscape(nth, 1)],[k_est_landscape(nth, 2) k_est_landscape(nth, 2)],[min(Jp(:)),max(Jp(:))],'Color', color{nth},'Linewidth',2);
        end        
        title('J landscape and initial local minima');
        subplot(1,2,2); axis xy;hold on; title('Zoom');
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14)
        for nth = 1:params.nMinima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20);
           % line([k_est_landscape(nth, 1) k_est_landscape(nth, 1)],[k_est_landscape(nth, 2) k_est_landscape(nth, 2)],[min(Jp(:)),max(Jp(:))],'Color', color{nth},'Linewidth',2);
        end  
        drawnow;
    end

    % Refine each of the found minima with  an iterative approach    
    disp('   - Refine position of extracted local min...'); 
    fprintf('%s','     - local min #');
    for ithk = 1:params.nMinima                  % Iterate minima
        if mod(ithk,ceil(params.nMinima/10))==0
            if ithk==params.nMinima
                fprintf('%i\n',ithk);
            else
                fprintf('%i, ',ithk);
            end
        end

%         G = y(:,:,idx);                        % Don't needed - we never modify G
        k=k_est_landscape(ithk, :); k = k(:);
        if displ>1
            fig=figure;
        end
        count_it = 1; 
        clear cf        

        for nit_filt = 1:params.FilterRefinement % Iterate filter refinement
            [att_filt, filt] = BuildFilter(k, sz, params);
    
            % Gradient descent parameters
            tau=1e-2;
            tau_fact=1.5;
            tau_min=1e-10;            
            nit_tot=100;
            cf_tolerance = 1e-6;
            
            % Calculate current cost
            cf(count_it) = EvalJ(k, wf, G, params, filt, att_filt, 0); %#ok<AGROW> 
            
            for jj=1:nit_tot   % Update k with nit_tot gradient descent steps            
                [~, g] = EvalJ(k, wf, G, params, filt, att_filt, 1);
                ktmp = k - tau * g';
                cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                
                % Line search by backtracking
                if cf(count_it) > cf_new % Increasing tau until minima, or 
                    while cf(count_it) > cf_new
                        k = ktmp;
                        tau = tau * tau_fact;
                        ktmp = k - tau * g';
                        cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                    end
                    % Reverse to previous ktmp and tau (as new cf is higher)
                    ktmp = k; tau = tau/tau_fact;
                else % Decreasing tau until minima
                    while cf(count_it)<=cf_new && (tau>tau_min)
                        tau = tau/tau_fact;
                        ktmp = k - tau * g';
                        cf_new = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);
                    end
                    if (tau>tau_min)
                        k=ktmp;
                    end
                end

                % Calculate new cost 
                count_it=count_it+1;
                cf(count_it) = EvalJ(ktmp, wf, G, params, filt, att_filt, 0);        
                
                % Stop grad descent if step size or cost improvement become negligible
                if (tau<tau_min) || abs(cf(count_it)-cf(count_it-1))/abs(cf(count_it-1))<cf_tolerance 
                    break;
                end        
            end
            
            % Draw curve before erasing cf
            if displ>1
                figure(fig); plot(cf,'linewidth',2); 
                xlim([0 nit_tot*params.FilterRefinement]);grid;
                title('Freq/Phase Opti: Cv curve'); set(gca,'FontSIze',14); drawnow;
            end
        end

        % Store final wavevector and phase (with updated filter)
        k_alt_est(ithk, :) = k; %#ok<AGROW> 
        [att_filt, filt] = BuildFilter(k, sz, params);
        A=BuildA(k, wf, filt, params); AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));       
        if params.method == 2
            s = AA\A'*G_filt(:);
        else
            G_filt_tmp= G_filt(:,:,1); 
            s = AA\A'*G_filt_tmp(:);
        end           
        ac(ithk)=s(1); as(ithk)=s(2);        %#ok<AGROW>         
    end

    disp('   - Choosing the best wavevector...'); % By reevaluating without filters
    lowestJ = realmax;
    for iii = 1:params.nMinima        
        Jp = EvalJ(k_alt_est(iii,:), wf, G, params, 0, 0, 0);
        if Jp < lowestJ  
            lowestJ = Jp;
            optIdx = iii;
        end
    end
    k_final(OrientCount, :) = k_alt_est(optIdx, :); %#ok<AGROW> Store in final array

    if displ
        figure(fg);subplot(1,2,1);
        plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        subplot(1,2,2);
        plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        axis([k_alt_est(optIdx, 1) - maxp*0.1 k_alt_est(optIdx, 1) + maxp*0.1 k_alt_est(optIdx, 2) - maxp*0.1 k_alt_est(optIdx, 2) + maxp*0.1 ]);grid;
        pause(0.1);
    end

    % Calculate phase and stick to convention
    tmp = atan(as(optIdx)/ac(optIdx));
    if ac(optIdx) <0, tmp=pi+tmp; end
    ph_init=mod(tmp,2*pi)/2; 
    phase(OrientCount, 1) = ph_init; %#ok<AGROW> 
    params.ac(OrientCount,1) = ac(optIdx); 
    params.as(OrientCount,1) = as(optIdx);
    
    % Extract amplitude factor (no filter used)
    A=BuildA(k_final(OrientCount, :), wf, filt, params); AA = A'*A;           
    if params.method == 2
        s = AA\A'*G(:);
    else
        G_tmp= G(:,:,1); s = AA\A'*G_tmp(:);
    end           
    ac(ithk)=s(1); as(ithk)=s(2); a(OrientCount, 1) = as(ithk)./sin(2*ph_init); 
    
    % We only calculate every phase if there is no assumption on the phase being made
    if params.method < 2
        k = k_alt_est(optIdx, :); 
        [att_filt, filt] = BuildFilter(k, sz, params);
        A=BuildA(k, wf, filt, params); AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));

        % For amplitude estimation
        ANoFilt = BuildA(k_final(OrientCount, :), wf, filt, params); 
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
            phase(OrientCount, i+1) = tmp;
            params.ac(OrientCount,i+1) = ac_tmp; 
            params.as(OrientCount,i+1) = as_tmp;
            % Calculate amplitude with the unfiltered as, but the accurate phase
            a(OrientCount, 1) = asNoFilt/sin(2*tmp);
        end
    end
%     fprintf("Parameters found for orientation %d:\n Kx = %5.4f | Ky = %5.4f | Ph = %5.4f\nk", ...
%         OrientCount, k_alt_est(optIdx, 1), k_alt_est(optIdx, 2), phase)
    OrientCount = OrientCount + 1;
%     display(["Parameters found for orientation", num2])
end

params.k = k_final; params.a = a; params.phase = phase; 
end