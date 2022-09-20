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
% Crop to the ROI, keeping the original size before
if isfield(params,'roi') && ~isempty(params.roi)   % Crop to ROI
    y = y(params.roi(1):params.roi(1)+params.roi(3)-1, ...
        params.roi(2):params.roi(2)+params.roi(3)-1,:);
end
sz = size(y);                                      % Calculate size of ROI
y = (y - min(y(:))) / (max(y(:)) - min(y(:)));     % Normalize stack images

% TODO : fix this convention once for all in FlexSIM.m (by reordering if
% necessary the stack y according to our code convention 'pa'.
imgIdxs = 1:params.nbPh*params.nbOr;               % Select the indexes to use (through imgIdxs)...
imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]); % according to the selected phase,z angle convention
if strcmp(params.AcqConv, "zap")                   % . If Zeiss convention (zap), flip
    imgIdxs = imgIdxs'; 
end 

% TODO: remove the adding of nbImgs in the user struct params
% TODO: improve regarding the choice of method
if params.method > 0                  % If there are assumptions on multiple...
    params.nbImgs = params.nbPh;      % images, we will use all the phases
else     
    imgIdxs = imgIdxs(1,:);           % Else we will only use the first... 
    params.nbImgs = 1;                % image of each orientation
end

% Common quantities of interest
[grids.I, grids.J] = meshgrid(0:sz(2)-1,0:sz(1)-1);        % Numerical mesh - multipurpose
grids.X = grids.I*params.res; grids.Y=grids.J*params.res;  % Scaled (by resolution) mesh 
OTF = GenerateOTF(params.Na, params.lamb, sz, params.res, 1); % Computation of the OTF

%% Loop over batch of images (1 batch = 1 orr + x phases)
% TODO: add back the displays (need to think how to manage the patch-based case)
OrientCount = 1; 
for idx = imgIdxs 
    disp([' Batch of images: ', num2str(idx')]);      % Display info to the user 
    
    disp('   - Remove WF and mask...');
    wf=mean(y(:,:,idx),3);                            % TODO : this is valid only when we have method 2 and we are sure that we can compute the wf like that
    [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);
    
    disp('   - Grid-based evaluation of J landscape...');
    [Jp,K1,K2] = GridEvalJ(params,wf,G,grids);
    
    disp(['   - Extracting the ',num2str(params.nMinima),' smallest local minima...']);
    k_init= ExtractLocMin(params,Jp,K1,K2);
    
    disp('   - Refine position of extracted local min...');
    fprintf('%s','     - local min #');
    for ithk = 1:params.nMinima
        if mod(ithk,ceil(params.nMinima/10))==0
            if ithk==params.nMinima, fprintf('%i\n',ithk);  else, fprintf('%i, ',ithk); end
        end
        k(ithk,:) = IterRefinementWavevec(k_init(ithk, :)',wf,G,grids,OTF,sz,params);
    end
    
    disp('   - Choosing the best wavevector...');    % Choose the best wavevector in terms of value of J   
    Jp=zeros(1,params.nMinima);
    for iii = 1:params.nMinima                       
        Jp(iii) = EvalJ(k(iii,:), wf, G, params, grids, 0, 0, 0);
    end
    [~,optIdx] = min(Jp);
    k_final(OrientCount, :) = k(optIdx, :); 
    
    disp('   - Computing phases and amplitutes...'); 
    % TODO: compute the phases accordingly to the choice of method
    [phase(OrientCount, 1),a(OrientCount, 1)]=GetPhaseAndAmp(k_final(OrientCount, :)',wf,G,grids,OTF,sz,params);
    
    OrientCount=OrientCount+1;
end

end

%% OLD CODE
%{
% Initializations
wfs =zeros(sz(1), sz(2),params.nbOr);              % Initialize variable for widefield

idwf=0;
for idx = imgIdxs 
    idwf=idwf+1;
    wf=mean(y(:,:,idx),3);     % TODO : this is valid only when we have method 2 and we are sure that we can compute the wf like that
    [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);
    y(:,:,idx')=G;
    wfs(:,:,idwf)=wf;          % TODO : same here, something to fix when not method==2
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

% -- Precomputations and evaluation of the cost function 
OrientCount = 1;                                      % Dummy var to count orientations
for idx = imgIdxs                                     % Loop through images to calculate J
    wf =  wfs(:,:,OrientCount);                       % Get corresponding wf
    disp([' Batch of images: ', num2str(idx')]);      % Display info to the user
    disp('   - Grid-based evaluation of J landscape...');
    G = y(:,:,idx);                                   % Select current batch of images
    [Jp,K1,K2] = GridEvalJ(params,wf,G,grids);
    
    disp(['   - Extracting the ',num2str(params.nMinima),' smallest local minima...']);   
    k_est_landscape= ExtractLocMin(params,Jp,K1,K2);
    
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

% -- Refinement of  each of the found minima with  an iterative approach    
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
        
        k = IterRefinementWavevec(k_est_landscape(ithk, :)',wf,G,grids,OTF,sz,params);     
        if params.displ>1                               % Draw curve before erasing cf
            fig=figure;
            figure(fig); plot(cf,'linewidth',2);
            xlim([0 nit_tot*params.FilterRefinement]);grid;
            xlabel('Iteration'); ylabel('Cost');
            title(sprintf('Phase Iter Opti Curve (Orient #%d, Minima #%d)', OrientCount, ithk)); set(gca,'FontSIze',14); drawnow;
        end

        % Store final wavevector and phase (with updated filter)
        k_alt_est(ithk, :) = k; %#ok<AGROW>     
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
    
    [phase(OrientCount, 1),a(OrientCount, 1)]=GetPhaseAndAmp(k_final(OrientCount, :)',wf,G,grids,OTF,sz,params);
    
    
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
            % Calculate amplitude with the unfiltered as, but the accurate phase
            a(OrientCount, 1) = asNoFilt/sin(2*tmp);
        end
    end
    OrientCount = OrientCount + 1;
end
%}