function [k_final, phase, a] = EstimatePatterns(params,y,k_init, wfAcq)
%--------------------------------------------------------------------------
% Function [k_final, phase, a] = EstimatePatterns(params,y,k_init)
%
% Estimate the parameters of 2D sinusoidal patterns of the form 
%    w(x) = 1 + a cos(<k,x> + phase),
% as described in [1].
%
% Inputs  : params  -> Structures with fields:
%                         - StackOrder: order of the SIM stack. Phase (p), angle (a) and time (z) convention. Choose one of ('paz', 'pza' or 'zap')
%                         - roi: region on interest on which the parameters will be estimated
%                         - lamb: Emission wavelength
%                         - Na: Objective numerical aperture
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
%           k_init  -> Initial wavevector. If provided, skips the grid search and runs only the
%                      iterative optimization. Set to 0 to perform grid search 
%           wfAcq   -> Widefield image. Provide if acquired. Otherwise set to 0                      
%
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
    if any(wfAcq, 'all')
        wfAcq = wfAcq(params.roi(1):params.roi(1)+params.roi(3)-1, ...
        params.roi(2):params.roi(2)+params.roi(3)-1,:);
    end
end
sz = size(y);                                      % Calculate size of ROI
y = (y - min(y(:))) / (max(y(:)) - min(y(:)));     % Normalize stack images

if ~k_init
    compute_k_init=true;
else
    compute_k_init=false;
end
 
% Check the acquisition convention of the user and convert to ap(w)
imgIdxs = 1:params.nbPh*params.nbOr;    % Select the indexes to use (through imgIdxs)...
if params.method                        % according to the selected phase,z angle convention
    imgIdxs = reshape(imgIdxs, [params.nbPh, params.nbOr]);  
end 
if params.method > 0                  % If there are assumptions on multiple...
    params.nbImgs = params.nbPh;      % images, we will use all the phases
else     
    imgIdxs = imgIdxs(1,:);           % Else we will only use the first... 
    params.nbImgs = 1;                % image of each orientation
end

% Common quantities of interest
[grids.I, grids.J] = meshgrid(0:sz(2)-1,0:sz(1)-1);           % Numerical mesh - multipurpose
grids.X = grids.I*params.res; grids.Y=grids.J*params.res;     % Scaled (by resolution) mesh 
OTF = GenerateOTF(params.Na, params.lamb, sz, params.res, 1); % Computation of the OTF
if ~any(wfAcq, 'all') && params.method < 2                   % If the WF has not been given and can be estimated 
    wfAcq  = mean(y,3);                                           % Calculate WF as the mean of all images
end

%% Loop over batch of images (1 batch = 1 orr + x phases)
% TODO: add back the displays (need to think how to manage the patch-based case)
OrientCount = 1; 

for idx = imgIdxs 
    disp([' Batch of images: ', num2str(idx')]);      % Display info to the user     
    disp('   - Remove WF and mask...');
%     if ~ismember('w', char(params.AcqConv)) && params.method == 2
    if ~any(wfAcq, 'all') && params.method == 2
        wf=mean(y(:,:,idx),3);                            % TODO : this is valid only when we have method 2 and we are sure that we can compute the wf like that
        [G,wf] = RemoveWFandMask(y(:,:,idx),wf,params);
    else
        [G,wf] = RemoveWFandMask(y(:,:,idx),wfAcq,params);
    end
    
%     if params.displ > 1 && OrientCount == 1           % Debug mode - display conditioned images
    if params.displ > 1                                 % Debug mode - display conditioned images
        if numel(size(G)) > 2
            figure; sliceViewer(G, "Colormap",viridis);   % Display masked `b` in [1]
            title('Conditioned Images');                  
            figure; sliceViewer(log10(abs(fftshift(fft2(G))+1)), "Colormap",viridis);   % Display FTs of both
            title('Conditioned Images FT');
        else
            figure; imdisp(G, 'Conditioned Images',0);   % Display masked `b` in [1]
            figure; imdisp(log10(abs(fftshift(fft2(G))+1)), 'Conditioned Images FT', 0);   % Display FTs of both
        end
        imdisp(wf, 'Conditioned Widefield', 1);  % Display widefields
        title('Conditioned widefields (masked on FT domain)');                       
        imdisp(log10(abs(fftshift(fft2(wf))+1)), "FT of Widefield", 1); % Display widefields
        title('Conditioned widefield FT');
    end
    
    if compute_k_init
        disp('   - Grid-based evaluation of J landscape...');
        [Jp,K1,K2] = GridEvalJ(params,wf,G,grids);

        

        if params.displ > 1                           % If requested, display J grid
            mJp=max(Jp(:));
            fg=figure; subplot(1,2,1); axis xy;hold on; view(45,45);
            surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
            colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14); %axis([0 maxp -maxp maxp]);grid;
%             for nth = 1:params.nMinima
%                 h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20); %#ok<AGROW>            
%             end        
            sgtitle(sprintf('J landscape for orientation #%d', OrientCount))
            title('Landscape and initial local minima');
            subplot(1,2,2); axis xy;hold on; title('Zoom');
            surf(K1,K2,Jp,'FaceColor','interp','EdgeColor','interp');
            colorbar; xlabel('k_1');ylabel('k_2'); set(gca,'fontsize',14)
%             for nth = 1:params.nMinima
%                 h(nth)=plot3(k_est_landscape(nth, 1),k_est_landscape(nth, 2),mJp,'Color', 'k','Marker', '.', 'markersize',20);           
%             end  
            drawnow;
        end
 
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
    else
        disp('   - Refine position of wavevector...');
        k_final(OrientCount, :) = IterRefinementWavevec(k_init(OrientCount, :)',wf,G,grids,OTF,sz,params);
    end
    
    disp('   - Computing phases and amplitutes...'); 
    % TODO: compute the phases accordingly to the choice of method
    [phase(OrientCount, :),a(OrientCount, :)]=GetPhaseAndAmp(k_final(OrientCount, :)',wf,G,grids,OTF,sz,params);
    
    OrientCount=OrientCount+1;
end
end