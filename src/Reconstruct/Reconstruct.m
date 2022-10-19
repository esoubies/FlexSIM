function rec = Reconstruct(y,patt,params)
%--------------------------------------------------------------------------
% function rec = Reconstruct(y,patt,params)
%
% SIM Reconstruction method described in [1].
%
% Inputs : y      - SIM data (stack)
%          patt   - Illumination patterns (stack, same order as y)
%          params - Structure with fields
%                     
% Output: rec - Super-resolved reconstructed image 
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also FlexSIM.m and EstimatePatterns.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Initial computations
sz=size(y);
szUp=sz*2;szUp=szUp(1:2); 
downFact=[2,2];
% -- Treat orientations separately or not
if params.sepOrr
    n1=params.nbOr;
    n2=params.nbPh;
else
    n1=1;
    n2=params.nbOr*params.nbPh;
end
% -- Generate OTF/PSF
otf = GenerateOTF(params.Na,params.lamb,min([256,256],sz(1:2)),params.res/2,params.damp);
if params.GPU
    otf = gpuArray(otf);
end
psf = fftshift(real(ifft2(otf)));
% -- Normalization
maxy=max(y(:));y=y/maxy;
sig=max(max(y(:,:,1)))/10;

%% Operators and Regul
% -- Data term
H=LinOpConv('PSF', psf,1,[1,2],'Centered','Pad',szUp+params.padSz,0);
P=LinOpSelectorPatch(szUp+params.padSz,[1 1],szUp);
S=LinOpDownsample(szUp,downFact);

% -- Regularization
G=LinOpGrad(P.sizeout,[1 2]); 
if params.regType==1      % Thikonov
    R=CostL2(G.sizeout)*G;
elseif params.regType==2  % TV
    R=CostHyperBolic(G.sizeout,1e-4,length(G.sizeout))*G;
elseif params.regType==3   % GR
    R=CostGoodRoughness(G,1e-1);
end

%% Cost and Optim
rec=zeros(szUp);
for id1=0:n1-1
    if params.sepOrr
        yy=y(:,:,id1*params.nbPh+1:(id1+1)*params.nbPh);
        pp=patt(:,:,id1*params.nbPh+1:(id1+1)*params.nbPh);
    else
        yy=y;
        pp=patt;
    end
    % -- Patterns normalization
    pp=pp/(mean(pp(:))*size(pp,3));
    % -- Data term
    wght=LinOpDiag([],1./(yy(:,:,1)+sig));
    if params.padSz==0
        F=CostL2([],yy(:,:,1),wght)*(S*H*LinOpDiag(P.sizeout,pp(:,:,1)));
    else
        F=CostL2([],yy(:,:,1),wght)*(S*P*H*P'*LinOpDiag(P.sizeout,pp(:,:,1)));
    end
    for id2=2:n2
        sig=max(max(yy(:,:,id2)))/10;
        wght=LinOpDiag([],1./(yy(:,:,id2)+sig));
        if params.padSz==0
            F=F+CostL2([],yy(:,:,id2),wght)*(S*H*LinOpDiag(P.sizeout,pp(:,:,id2)));
        else
            F=F+CostL2([],yy(:,:,id2),wght)*(S*P*H*P'*LinOpDiag(P.sizeout,pp(:,:,id2)));
        end
    end
    
    % -- Build cost and optimize. Try the faster VMLMB for Linux devices, or FBS for Windows devices 
    try        
        Opt=OptiVMLMB((1/numel(yy))*F+(params.mu/prod(P.sizeout))*R,0.,[]);            % Algorithm instanciation
    catch
        disp(['Reconstructing with Projected Gradient Descent...']);
        Opt=OptiFBS((1/numel(yy))*F+(params.mu/prod(P.sizeout))*R, CostNonNeg(R.sizein));
        Opt.fista=1;
        Opt.gam=numel(yy);
        Opt.updateGam='backtracking';
        Opt.eta=1.5;
    end
    Opt.OutOp=OutputOpti(1,round(params.maxIt/10));  % Verbose monitoring
    Opt.CvOp=TestCvgStepRelative(params.stepTol);    % Test convergence criterion
    Opt.ItUpOut=round(params.maxIt/10);              % Call OutputOpti update every ItUpOut iterations
    Opt.maxiter=params.maxIt;                        % Max number of iterations
    Opt.run(zeros(P.sizeout));                       % Run the algorithm zeros(H.sizein)
    rec=rec+Opt.xopt*maxy;
end
rec=rec/n1;
end