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
% if params.GPU
%     useGPU(1)
% else
%     useGPU(0)
% end

sz=size(y);
szUp=sz*2;szUp=szUp(1:2); 
downFact=[2,2];
% -- Treat orientations separately or not
if params.sepOrr
    n1=params.nbOr;
else
    n1=1;
end
% -- Generate OTF/PSF
otf = GenerateOTF(params.Na,params.lamb,min([256,256],sz(1:2)),params.res/2,params.damp);
% if params.GPU
%     otf = gpuArray(otf);
% end
psf = fftshift(real(ifft2(otf)));
% -- Normalization
maxy=max(y(:));y=y/maxy;
% -- OTF Attenuation
if params.OTFAttStr && params.OTFAttwdth
    OTFatt = GenerateOTFAttMask(params.Na,params.lamb,sz(1:2),params.res,params.OTFAttStr,params.OTFAttwdth);
    Hatt=LinOpConv('MTF', OTFatt,1,[1,2]);
else
    Hatt=LinOpIdentity(sz(1:2));
end
% -- Apodisation function
if params.apodize
    [X,Y]=meshgrid(linspace(-1,1,sz(2)),linspace(-1,1,sz(1)));
    Apo=1./(1+exp(100*(abs(X)-0.97)))./(1+exp(100*(abs(Y)-0.97)));
else
    Apo=1;
end

%% Operators and Regul
% -- Data term
H=LinOpConv('PSF', psf,1,[1,2],'Centered','Pad',szUp+params.padSz,0);
if params.padSz==0
    P=LinOpIdentity(szUp);
else
    P=LinOpSelectorPatch(szUp+params.padSz,[1 1],szUp);
end
S=LinOpDownsample(szUp,downFact);

% -- Regularization
G=LinOpGrad(P.sizein,[1 2]);
if params.regType==1      % Thikonov
    R=1/prod(G.sizein)*CostL2(G.sizeout)*G;
elseif params.regType==2  % TV
    R=1/prod(G.sizein)*CostHyperBolic(G.sizeout,1e-4,length(G.sizeout))*G;
elseif params.regType==3   % GR
    R=1/prod(G.sizein)*CostGoodRoughness(G,1e-1);
end

%% Cost and Optim
rec=zeros([szUp(1:2),n1]);
x=zeros(P.sizein);
if params.parallelProcess && (params.szPatch==0) && (n1>1)
    nbcores=feature('numcores');
    parfevalOnAll(@warning,0,'off','all');
else
    nbcores=0;
end

if params.sepOrr
    y=reshape(y,[size(y,[1,2]),params.nbPh,params.nbOr]);
    patt=reshape(patt,[size(patt,[1,2]),params.nbPh,params.nbOr]);
end
parfor (id1 = 1:n1,nbcores)
    if nbcores>0
        t = getCurrentTask();
        DispMsg(params.verbose,['-- [Worker #',num2str(t.ID),'] Reconstruct orientation  #',num2str(id1),'/',num2str(params.nbOr),' ...']);
    end

    % -- Patterns normalization
    pp=patt(:,:,:,id1);yy=y(:,:,:,id1);
    pp=pp/(mean(pp(:))*size(pp,3));
    % -- Data term
    sig=max(max(yy,[],1),[],2)/10;
    wght=LinOpDiag([],Apo./(yy(:,:,1)+sig(1)));
    F=CostL2([],yy(:,:,1),wght*Hatt)*(S*P*H*LinOpDiag(P.sizein,pp(:,:,1)));
    for id2=2:size(yy,3)
        wght=LinOpDiag([],Apo./(yy(:,:,id2)+sig(id2)));
        F=F+CostL2([],yy(:,:,id2),wght*Hatt)*(S*P*H*LinOpDiag(P.sizein,pp(:,:,id2)));
    end

    % -- Build cost and optimize.
    Opt=OptiVMLMB((1/numel(yy))*F+params.mu*R,0.,[]);            % Algorithm instanciation
    Opt.OutOp=OutputOpti(1,round(params.maxIt/10));  % Verbose monitoring
    Opt.verbose=(params.verbose==2)*(~params.parallelProcess);
    Opt.CvOp=TestCvgStepRelative(params.stepTol);    % Test convergence criterion
    Opt.ItUpOut=round(params.maxIt/10);              % Call OutputOpti update every ItUpOut iterations
    Opt.maxiter=params.maxIt;                        % Max number of iterations
    Opt.run(x);                        % Run the algorithm zeros(H.sizein)
    rec(:,:,id1)=P*Opt.xopt*maxy;
    if nbcores>0
        DispMsg(params.verbose,['-- [Worker #',num2str(t.ID),'] Reconstruct orientation  #',num2str(id1+1),'/',num2str(params.nbOr),' done.']);
    end
end
rec=mean(rec./mean(rec,[1,2]),3);

% Apodize result (as apotization is used in the cost)
if params.apodize
    [X,Y]=meshgrid(linspace(-1,1,szUp(2)),linspace(-1,1,szUp(1)));
    Apo=1./(1+exp(100*(abs(X)-0.97)))./(1+exp(100*(abs(Y)-0.97)));
    rec=rec.*Apo;
end
end