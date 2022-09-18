function [k_final, phase, a] = EstimatePatterns(params,y)
%--------------------------------------------------------------------------
% Function [k_final, phase, a] = EstimatePatterns(params,y)
%
% Estimate the parameters of 2D sinusoidal patterns of the form 
%    w(x) = 1 + a cos(<k,x> + phase),
% as described in [1].
%
% Inputs  : params  -> Structures with fields:
%                         - AcqConv: order of the SIM stack. Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
%                         - roi: region on interest on which the parameters will be estimated
%                         - lamb: Emission wavelength
%                         - Na: Objective numerica aperture
%                         - res: resolution of the SIM data stac
%                         - nbOr: number of orientations
%                         - nbPh: number of phases
%                         - method: method used to estimate the parameters, 3 choices
%                                    0 - treat all images independently
%                                    1 - use all images with same orientation to estimate a unique wavevector
%                                    2 - 1 + assume equally spaced phases
%                         - ringMaskLim: 1x2 array (eg. [0.3, 1.1]) defining the ring used to mask the WF component and the high-freq of the data, givien as factor of fc=1
%                         - limits: 1x2 array (eg. [0.9, 1.1]) defining  the ring over which the J function is evaluated for initializing, givien as factor of fc=1
%                         - nMinima: Number of starting points for the refinement steps
%                         - nPoints: Number of points in the J evaluation grid
%                         - FilterRefinement: Number of times that the filter is upgraded (gradient descent cycles)
%                         - displ: if >1, displays intermediate results
%           y       -> SIM data stack
%
% Outputs : k_final -> array with the estimated wavevectors (if method=0, of size nbOr*nbPh. Otherwise, of size nbOr)
%           phase   -> array with the estimated phases (if method=0,1, of size nbOr*nbPh. Otherwise, of size nbOr)
%           a       -> Same size as phase, contains the amplitude of the SIM patterns
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also FlexSIM.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Preprocessing of stack & initialization of useful variables
% If requested, we crop to the ROI, keeping the original size before
if isfield(params,'roi') && ~isempty(params.roi)   % Crop to ROI
    y = y(params.roi(1):params.roi(1)+params.roi(3)-1, ...
        params.roi(2):params.roi(2)+params.roi(3)-1,:);
end
sz = size(y);                                      % Calculate size of ROI

[grids.I, grids.J] = meshgrid(0:sz(2)-1,0:sz(1)-1);        % Numerical mesh - multipurpose
grids.X = grids.I*params.res; grids.Y=grids.J*params.res;  % Scaled (by resolution) mesh 

OTF = GenerateOTF(params.Na, params.lamb, sz, params.res, 1); % Computation of the OTF

wfs =zeros(sz(1), sz(2),params.nbOr);              % Initialize variable for widefield
y = (y - min(y(:))) / (max(y(:)) - min(y(:)));     % Normalize stack images

imgIdxs = 1:params.nbPh*params.nbOr;               % Select the indexes to use (through imgIdxs)...
imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]); % according to the selected phase,z angle convention
if strcmp(params.AcqConv, "zap")                   % . If Zeiss convention (zap), flip
    imgIdxs = imgIdxs'; 
end 
                            
if params.method > 0                  % If there are assumptions on multiple...
    params.nbImgs = params.nbPh;      % images, we will use all the phases
else     
    imgIdxs = imgIdxs(1,:);           % Else we will only use the first... 
    params.nbImgs = 1;                % image of each orientation
end

%% Processing of stack - Removal of widefield and masking
idwf=0;
FCut = 2*params.Na/params.lamb*params.res;    % Initialize cut-off frequency
for idx = imgIdxs                             % Iterate through the orientations using the corresponding batch
    idwf=idwf+1;
    % Compute the wf and fft per batch of images (orientation)
    wfs(:,:,idwf)=sum(y(:,:,idx), 3)/length(idx); fft_wf = fftshift(fft2(wfs(:,:,idwf)));
    for img = idx'                            % Iterate images in the batch
        ffty=fftshift(fft2(y(:,:,img)));      % Calculate FFT
        [I,J]=meshgrid(1:sz(1),1:sz(2));      % Create a LP filter in Fourier domain to remove WF
        p = (sz(1:2)+1)/2; r = sz(1)/8; 
        mask=double(((I-p(1)).^2+(J-p(2)).^2) < r^2); 
        a=real(OptWght(ffty,fft_wf, mask));   % Calculate argmin_a |a*ffty - fft_wf|^2
        % Remove the scaled WF + filter freq within a ring of interest
        y(:,:,img) = real(ifft2(ifftshift(MaskFT((a*ffty-fft_wf), FCut, params.ringMaskLim)))); 
    end
end

for ii=1:size(wfs,3)                           % Mask widefields with the same ring as the images
    wfs(:,:,ii) = real(ifft2(ifftshift(MaskFT(fftshift(fft2(wfs(:,:,ii))), FCut, params.ringMaskLim))));
end

if params.displ > 1                               % Debug mode - display conditioned images
    figure; sliceViewer(y, "Colormap",viridis);   % Display masked `b` in [1]
    title('Conditioned Images');                  
    figure; sliceViewer(wfs, "Colormap",viridis); % Display widefields
    title('Conditioned widefields (masked on FT domain)');
    figure; sliceViewer(log10(abs(fftshift(fft2(y))+1)), "Colormap",viridis);   % Display FTs of both
    title('Conditioned Images FT');               
    figure; sliceViewer(log10(abs(fftshift(fft2(wfs))+1)), "Colormap",viridis); % Display widefields
    title('Conditioned widefields FT');
end
%% Precomputations and evaluation of the cost function J
cmin=params.limits(1);  cmax=params.limits(2);        % Ring used to delimit Fourier sin peaks search
FCutN=FCut*pi/params.res;                             % Scaled cutoff frequency
maxp=FCut*pi/params.res*max(1,cmax);                  % Maximum value of one component of k
cen = [0,0];                                          % Grid center
ll1=linspace(cen(1),cen(1)+maxp,params.nPoints/2);    % Declare grid of wavevectors...
ll2=linspace(cen(2)-maxp,cen(2)+maxp,params.nPoints); % with the x direction halved 
[K1,K2]=meshgrid(ll1,ll2);                            
K = cat(3, K1, K2); 
Knorm = vecnorm(K, 2, 3);                             % Calculate norms in advance

OrientCount = 1;                                      % Dummy var to count orientations
for idx = imgIdxs                                     % Loop through images to calculate J
    wf =  wfs(:,:,OrientCount);                       % Get corresponding wf
    disp([' Batch of images: ', num2str(idx')]);      % Display info to the user 
    disp('   - Grid-based evaluation of J landscape...');
    Jp = zeros(params.nPoints, params.nPoints/2);     % Initialize Jp for batch
    G = y(:,:,idx);                                   % Select current batch of images

    for ii=1:params.nPoints                           % Iterate the grid and if the wavevector
        for jj=1:params.nPoints/2                     % is within the limits, evaluate
            if Knorm(ii, jj)>cmin*FCutN && Knorm(ii, jj)<cmax*FCutN
                ktest = squeeze(K(ii, jj, :));        % Evaluate cost without filtering                 
                Jp(ii, jj) = EvalJ(ktest, wf, G, params, grids, 0, 0, 0);
            else
                Jp(ii,jj)=NaN;
            end
        end
    end
    
    disp(['   - Extracting the ',num2str(params.nMinima),' smallest local minima...']);
    k_est_landscape = zeros(params.nMinima, 2);       % Give user info and initialize...
    Jp_tmp = Jp;                                      % variables to get minima of J

    % Analyze Jp (through Jp_tmp) to extract the n local minima
    for nth = 1:params.nMinima
        [~,idxMin]=min(Jp_tmp(:));                     % Extract and store minima
        [ii,jj]=ind2sub([params.nPoints,params.nPoints/2],idxMin);
        k_est_tmp=[K1(ii,jj);K2(ii,jj)];               
        k_est_landscape(nth, :) = k_est_tmp;

        % Temporarily erase the vecinity of the found minima
        idx1 = ii - round(params.nPoints/50); idx1(idx1<1) = 1;
        idx2 = ii + round(params.nPoints/50); idx2(idx2<params.nPoints) = params.nPoints;
        idx3 = jj - round(params.nPoints/50); idx3(idx3<1) = 1;
        idx4 = jj + round(params.nPoints/50); idx4(idx4>params.nPoints/2) = params.nPoints/2;
        Jp_tmp(idx1:idx2, idx3:idx4)  = max(Jp(:));
    end
    
    if params.displ > 1                           % If requested, display J grid
        mJp=max(Jp(:));
        fg=figure; subplot(1,2,1); axis xy;hold on; view(45,45);
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14); axis([0 maxp -maxp maxp]);grid;
        for nth = 1:params.nMinima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20); %#ok<AGROW>            
        end        
        sgtitle(sprintf('J landscape for orientation #%d', OrientCount))
        title('Landscape and initial local minima');
        subplot(1,2,2); axis xy;hold on; title('Zoom');
        surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
        colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14)
        for nth = 1:params.nMinima
            h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20);           
        end  
        drawnow;
    end

%% Refinement of  each of the found minima with  an iterative approach    
    disp('   - Refine position of extracted local min...');  % Display user info
    fprintf('%s','     - local min #');
    for ithk = 1:params.nMinima                    % Iterate minima
        if mod(ithk,ceil(params.nMinima/10))==0
            if ithk==params.nMinima
                fprintf('%i\n',ithk);              % More user info
            else
                fprintf('%i, ',ithk);
            end
        end

        k=k_est_landscape(ithk, :); k = k(:);      % Select wavevector
        if params.displ>1                          % More display (optim curve)
            fig=figure;
        end
        count_it = 1;                              % Dummy variable for plotting and storing cost
        clear cf                                   % Clear cost array for the current iteration

        for nit_filt = 1:params.FilterRefinement   % Iterate filter refinement
            [att_filt, filt] = BuildFilter(k, sz, OTF, params, grids);
    
            tau=1e-2;                              % Gradient descent parameters
            tau_fact=1.5;
            tau_min=1e-10;            
            nit_tot=100;
            cf_tolerance = 1e-6;
            
            % Calculate current cost
            cf(count_it) = EvalJ(k, wf, G, params, grids, filt, att_filt, 0); %#ok<AGROW> 
            
            for jj=1:nit_tot                       % Update k with gradient descent (max nit_tot steps)
                [~, g] = EvalJ(k, wf, G, params, grids, filt, att_filt, 1);    % Get current gradient
                ktmp = k - tau * g';                                             % Update k
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

                count_it=count_it+1;                  % Calculate and store new cost 
                cf(count_it) = EvalJ(ktmp, wf, G, params, grids, filt, att_filt, 0);        
                
                % Stop grad descent if step size or cost improvement become negligible
                if (tau<tau_min) || abs(cf(count_it)-cf(count_it-1))/abs(cf(count_it-1))<cf_tolerance 
                    break;
                end        
            end
                        
            if params.displ>1                               % Draw curve before erasing cf
                figure(fig); plot(cf,'linewidth',2); 
                xlim([0 nit_tot*params.FilterRefinement]);grid;
                xlabel('Iteration'); ylabel('Cost');
                title(sprintf('Phase Iter Opti Curve (Orient #%d, Minima #%d)', OrientCount, ithk)); set(gca,'FontSIze',14); drawnow;
            end
        end

        % Store final wavevector and phase (with updated filter)
        k_alt_est(ithk, :) = k; %#ok<AGROW> 
        [att_filt, filt] = BuildFilter(k, sz, OTF, params, grids);
        A=BuildA(k, wf, filt, params, grids); AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));       
        if params.method == 2
            s = AA\A'*G_filt(:);
        else
            G_filt_tmp= G_filt(:,:,1); 
            s = AA\A'*G_filt_tmp(:);
        end           
        ac(ithk)=s(1); as(ithk)=s(2);        %#ok<AGROW>         
    end

    disp('   - Choosing the best wavevector...');    % Iterate the found minimas
    lowestJ = realmax;                               % And choose the one with the 
    for iii = 1:params.nMinima                       % lowest cost
        Jp = EvalJ(k_alt_est(iii,:), wf, G, params, grids, 0, 0, 0);
        if Jp < lowestJ  
            lowestJ = Jp;
            optIdx = iii;
        end
    end
    k_final(OrientCount, :) = k_alt_est(optIdx, :); %#ok<AGROW> Store in final array
    

    if params.displ > 1                                     
        figure(fg);subplot(1,2,1);                   % Paint minima in red
        plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        subplot(1,2,2);
        h = plot3(k_alt_est(optIdx, 1),k_alt_est(optIdx, 2),mJp,'Color', 'r','Marker', '.', 'markersize',25);
        axis([k_alt_est(optIdx, 1) - maxp*0.1 k_alt_est(optIdx, 1) + maxp*0.1 k_alt_est(optIdx, 2) - maxp*0.1 k_alt_est(optIdx, 2) + maxp*0.1 ]);grid;
        legend(h, 'Optimized wavevector')
        pause(0.1);

        kPix = k_alt_est(optIdx, :)' .* sz(1:2)' * params.res / pi;

        % Show corresponding filter
        [att_filt, filt] = BuildFilter(k_alt_est(optIdx, :)', sz, OTF, params, grids);        
        figure(); subplot(1,2,1); imshow(fftshift(att_filt), []); colormap(viridis); impixelinfo; 
        title('Filter $G_A$ for estimated $k$', 'interpreter', 'latex'); hold on
        plot(kPix(1) + sz(1)/2 + 1, kPix(2)+ sz(2)/2 + 1, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
        subplot(1,2,2); imshow(fftshift(filt), []); colormap(viridis); impixelinfo; hold on
        title('Filter $G_b$', 'interpreter', 'latex')
        plot(kPix(1) + sz(1)/2 + 1, kPix(2)+ sz(2)/2 + 1, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
        sgtitle(sprintf('Filters used for the detected wavevector (Orient #%d)',OrientCount))
        drawnow
    end

    tmp = atan(as(optIdx)/ac(optIdx));           % Calculate phase, ac, and as
    if ac(optIdx) <0, tmp=pi+tmp; end
    ph_init=mod(tmp,2*pi)/2;                     
    phase(OrientCount, 1) = ph_init;             % #ok<AGROW> % Store 
    A=BuildA(k_final(OrientCount, :), wf, 0, params, grids);
    AA = A'*A;                                   % Recalculate ac and as without...
    if params.method == 2                        % filters and extract amplitude
        s = AA\A'*G(:);
    else
        G_tmp= G(:,:,1); s = AA\A'*G_tmp(:);
    end           
    ac(ithk)=s(1); as(ithk)=s(2); a(OrientCount, 1) = as(ithk)./sin(2*ph_init); 
    
    % Calculate every phase if there is no assumption equidistant phase
    if params.method < 2
        k = k_alt_est(optIdx, :); 
        [att_filt, filt] = BuildFilter(k, sz, OTF, params, grids);
        A=BuildA(k, wf, filt, params, grids); AA = A'*A;
        G_filt=real(ifft2(fft2(G).*att_filt));

        % For amplitude estimation
        ANoFilt = BuildA(k_final(OrientCount, :), wf, filt, params, grids); 
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
    OrientCount = OrientCount + 1;
end
end