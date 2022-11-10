function [prm] = EvalRun(prm)
%--------------------------------------------------------------------------
% Function prm = EvalRun(prm)
% 
% Returns a structure that contains the results of the simulation. This
% inlcudes the accuracy of the parameter estimation, and the SNR of the 
% reconstruction as measured against the original image.
%
% Inputs : prm -> Structure with fields  
%                         - or: Orientations of the SIM images [rad]
%                         - ph: Phases of the SIM
%                         - DataPath: Path of the original image
%                         - k: Estimated wavevectors 
%                         - phase: Estimated phases
%                         - rec: reconstruction 
%           plot  -> Boolean to choose whether to plot or not the existing
%                    results. 
%
% Outputs: prm    -> Input structure with the additional fields
%                         - kGt: Groundtruth wavevector of the acquisition
%                         - kErr: Wavevector error expressed as norm
%                         - phGt: Groundtruth phase of the acqu
%                         - phErr: Phase error expressed as norm
%                         - MEP: Maximum expected number of photons during acq
%                         - snr: SNR of the reconstruction
%                         - ssi: SSI of the reconstruction
%
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

% Few initial tests and useful definitions

FCut = 2*prm.Na/prm.lamb*prm.res;    % Cut-off frequency

for i = 2:prm.nbPh                % Complete the phase values in the CCEqPh, J and Filt case
    prm.phCCEqPh(:,i) = prm.phCCEqPh(:,i-1) + pi/prm.nbPh;
    prm.phCCEqPhRef(:,i) = prm.phCCEqPhRef(:,i-1) + pi/prm.nbPh;
    prm.phCCEqPhFilt(:,i) = prm.phCCEqPhFilt(:,i-1) + pi/prm.nbPh;
%     prm.phJ(:,i) = prm.phJ(:,i-1) + pi/prm.nbPh;
%     prm.phFilt(:,i) = prm.phFilt(:,i-1) + pi/prm.nbPh;
end
prm.ph = repmat(prm.ph, [prm.nbOr, 1]);                 % Complete the GTarray for easy access
% prm.ph = padarray(prm.ph, [5-prm.nbPh 0], NaN, 'post'); % Pad with NaN if less than 5 phases (for table storage)
% prm.phCCEqPh = padarray(prm.phCCEqPh, [5-prm.nbPh 0], NaN, 'post');
% prm.phCC = padarray(prm.phCC, [5-prm.nbPh 0], NaN, 'post');
for or = 1:prm.nbOr                           % Iterate orientations (each will be processed as one)
    % GT wavevectorCalculate
    k = 2*pi*prm.ns/prm.lamb*[cos(prm.or(or)), sin(prm.or(or))]*prm.Na/prm.nl;
    prm.k = k;

    % Generate patterns for comparison
    [X,Y]=meshgrid(0:prm.sz(2)-1,0:prm.sz(1)-1); X=X*prm.res; Y=Y*prm.res;
    patt = ones([prm.sz(1:2) prm.nbPh]);    
    gen_patt = @(a, k,p) 1 + a*cos(2*p)*cos(2*(k(1)*X+k(2)*Y)) - a*sin(2*p)*sin(2*(k(1)*X+k(2)*Y));
    for jj = 1:prm.nbPh                                      % Iterate orientation and phase
        patt(:,:,jj) =   gen_patt(prm.a, k, prm.ph(jj, or));
    end

    % Estimate parameters from pattern (to test resolution limit)
    fac = 4;                                                              % 
    for jj=1:prm.nbPh
        fftpatt=fft2(padarray(patt(:,:,jj),size(patt(:,:,jj))*fac,'post'));   % pad and compute fft
        fftpattMask=MaskFT(fftshift(fftpatt),FCut,[0.3 1]);                   % mask the central peak
        [~,id]=max(abs(fftpattMask(:)));                                      % detect one of the two remaining peaks
        [i,j]=ind2sub(size(fftpattMask),id);     
        if all(sign([j,i]-size(fftpattMask)/2)==sign(k)), sg=1; else sg=-1; end  % to know if we detected the one with same sign as simulated k (if not need to change the sign of the arg in the next line)
        prm.kPatt=sg*([j,i]-floor(size(fftpattMask)/2)-1)*pi/prm.res./size(fftpattMask);        
        prm.phPatt(jj)=mod(sg*(angle(fftpattMask(id))),2*pi)/2;               % Simple arg{wavevecto}
        [~, idx_tmp] = min([abs(prm.phPatt(jj)-prm.ph(jj, or)),abs(prm.phPatt(jj)+pi-prm.ph(jj, or)),abs(prm.phPatt(jj)-prm.ph(jj, or)-pi)]);        
        switch idx_tmp
            case 1
            case 2
                prm.phPatt(jj) = prm.phPatt(jj) + pi;               
            case 3
                prm.phPatt(jj) = prm.phPatt(jj) - pi;
        end
    end       

    % Pad with NaNs the missing phase values (we're storing up to 5)
%     prm.phPatt = padarray(prm.phPatt, [5-prm.nbPh 0], NaN, 'post'); 

    % Generate Patterns and calculate Frob Norm
    froPatt = zeros(3, 1);
    froCCEqPh = zeros(3, 1);
    froCCEqPhRef = zeros(3, 1);
    froCCEqPhFilt = zeros(3, 1);
    froCC = zeros(3, 1);
    froCCRef = zeros(3, 1);
    froCCFilt = zeros(3, 1);
    fro= zeros(3, 1);

    for jj = 1:prm.nbPh                  % Iterate orientation and phase
        % Build the corresponding patterns
        pattPatt = gen_patt(prm.a, prm.kPatt, prm.phPatt(jj));       % From Pattern
        pattCCEqPh = gen_patt(prm.a, prm.kCCEqPh(or, :), prm.phCCEqPh(or, jj)); % CC EqPh
        pattCCEqPhRef = gen_patt(prm.a, prm.kCCEqPhRef(or, :), prm.phCCEqPhRef(or, jj)); % CC EqPh Ref
        pattCCEqPhFilt = gen_patt(prm.a, prm.kCCEqPhFilt(or, :), prm.phCCEqPhFilt(or, jj)); % CC EqPh Filt
        pattCC = gen_patt(prm.a, prm.kCC(or, :), prm.phCC(or, jj));             % CC
        pattCCRef = gen_patt(prm.a, prm.kCCRef(or, :), prm.phCCRef(or, jj));             % CC Ref
        pattCCFilt = gen_patt(prm.a, prm.kCCFilt(or, :), prm.phCCFilt(or, jj));             % CC Filt       
        % And the corresponsing Frobenius norms
        fro(jj) = 0;
        froPatt(jj) = norm(pattPatt - patt(:,:,jj),'fro');
        froCCEqPh(jj) = norm(pattCCEqPh - patt(:,:,jj),'fro');
        froCCEqPhRef(jj) = norm(pattCCEqPhRef - patt(:,:,jj),'fro');
        froCCEqPhFilt(jj) = norm(pattCCEqPhFilt - patt(:,:,jj),'fro');
        froCC(jj) = norm(pattCC - patt(:,:,jj),'fro');
        froCCRef(jj) = norm(pattCCRef- patt(:,:,jj),'fro');
        froCCFilt(jj) = norm(pattCCFilt - patt(:,:,jj),'fro');
    end

    % Table variable names
%     VarNames = ["Contrast", "MEP", "Kx (GT)", "Ky (GT)", ...
%                                      "Kx (From Patt)", "Ky (From Patt)", ... %"Kx (Peak D)", "Ky (Peak D)", ...
%                                      "Kx (CC eq-ph)", "Ky (CC eq-ph)",  ...
%                                      "Kx (CC)", "Ky (CC)", ...
%                                      "Ph #1 (GT)", "Ph #1 (from Patt)", "Ph #1 (CC eq-ph)", "Ph #1 (CC)", ...
%                                      "Ph #2 (GT)", "Ph #2 (from Patt)", "Ph #2 (CC eq-ph)", "Ph #2 (CC)", ...
%                                      "Ph #3 (GT)", "Ph #3 (from Patt)", "Ph #3 (CC eq-ph)", "Ph #3 (CC)", ...
%                                      "Ph #4 (GT)", "Ph #4 (from Patt)", "Ph #4 (CC eq-ph)", "Ph #4 (CC)", ...
%                                      "Ph #5 (GT)", "Ph #5 (from Patt)", "Ph #5 (CC eq-ph)", "Ph #5 (CC)",...
%                                      "FrobNorm #1 (fromPatt)", "FrobNorm #1 (CC eq-ph)", "FrobNorm #1 (CC)",...
%                                      "FrobNorm #2 (fromPatt)", "FrobNorm #2 (CC eq-ph)", "FrobNorm #2 (CC)",...
%                                      "FrobNorm #3 (fromPatt)", "FrobNorm #3 (CC eq-ph)", "FrobNorm #3 (CC)",...
%                                      "FrobNorm #4 (fromPatt)", "FrobNorm #4 (CC eq-ph)", "FrobNorm #4 (CC)",...
%                                      "FrobNorm #5 (fromPatt)", "FrobNorm #5 (CC eq-ph)", "FrobNorm #5 (CC)"];
% 
%     % Load table to store the results
%     T = table(prm.a, prm.MEP, prm.k(1), prm.k(2), prm.kPatt(1), prm.kPatt(2), ...
%     prm.kCCEqPh(or, 1), prm.kCCEqPh(or, 2), prm.kCC(or, 1), prm.kCC(or, 2), ...
    VarNames = ["Contrast", "MEP", "K (GT)", ...
                 "K (From Patt)", "K (CC eq-ph)", "K (CC eq-ph Ref)", "K (CC eq-ph Filt)", "K (CC)", "K (CC Ref)", "K (CC Filt)",...
                 "Ph #1 (GT)", "Ph #1 (from Patt)", "Ph #1 (CC eq-ph)", "Ph #1 (CC eq-ph Ref)", "Ph #1 (CC eq-ph Filt)",...
                 "Ph #1 (CC)", "Ph #1 (CC Ref)", "Ph #1 (CC Filt)", ...
                 "Ph #2 (GT)", "Ph #2 (from Patt)", "Ph #2 (CC eq-ph)", "Ph #2 (CC eq-ph Ref)", "Ph #2 (CC eq-ph Filt)",... 
                 "Ph #2 (CC)", "Ph #2 (CC Ref)", "Ph #2 (CC Filt)", ...
                 "Ph #3 (GT)", "Ph #3 (from Patt)", "Ph #3 (CC eq-ph)", "Ph #3 (CC eq-ph Ref)", "Ph #3 (CC eq-ph Filt)",... 
                 "Ph #3 (CC)", "Ph #3 (CC Ref)", "Ph #3 (CC Filt)", ...                 
                 "FrobNorm #1 (GT)", "FrobNorm #1 (from Patt)", "FrobNorm #1 (CC eq-ph)", "FrobNorm #1 (CC eq-ph Ref)", "FrobNorm #1 (CC eq-ph Filt)",...
                 "FrobNorm #1 (CC)", "FrobNorm #1 (CC Ref)", "FrobNorm #1 (CC Filt)", ...
                 "FrobNorm #2 (GT)", "FrobNorm #2 (from Patt)", "FrobNorm #2 (CC eq-ph)", "FrobNorm #2 (CC eq-ph Ref)", "FrobNorm #2 (CC eq-ph Filt)",...
                 "FrobNorm #2 (CC)", "FrobNorm #2 (CC Ref)", "FrobNorm #2 (CC Filt)", ...
                 "FrobNorm #3 (GT)", "FrobNorm #3 (from Patt)", "FrobNorm #3 (CC eq-ph)", "FrobNorm #3 (CC eq-ph Ref)", "FrobNorm #3 (CC eq-ph Filt)",...
                 "FrobNorm #3 (CC)", "FrobNorm #3 (CC Ref)", "FrobNorm #3 (CC Filt)"];

    % Create table with results. Wavevector is converted to positive x convention
    T = table(prm.a, prm.MEP, ...
        sign(prm.k(1))*prm.k, ...
        sign(prm.kPatt(1))*prm.kPatt, ...
        sign(prm.kCCEqPh(or, 1))*prm.kCCEqPh(or, :), ...
        sign(prm.kCCEqPhRef(or, 1))*prm.kCCEqPhRef(or, :), ...
        sign(prm.kCCEqPhFilt(or, 1))*prm.kCCEqPhFilt(or, :), ...
        sign(prm.kCC(or, 1))*prm.kCC(or, :), ...
        sign(prm.kCCRef(or, 1))*prm.kCCRef(or, :), ...
        sign(prm.kCCFilt(or, 1))*prm.kCCFilt(or, :), ...
        prm.ph(1, or), prm.phPatt(1),  prm.phCCEqPh(1, or), prm.phCCEqPhRef(1, or), prm.phCCEqPhFilt(1, or), prm.phCC(1, or), prm.phCCRef(1, or), prm.phCCFilt(1, or),...
        prm.ph(2, or), prm.phPatt(2),  prm.phCCEqPh(2, or), prm.phCCEqPhRef(2, or), prm.phCCEqPhFilt(2, or), prm.phCC(2, or), prm.phCCRef(2, or), prm.phCCFilt(2, or),...
        prm.ph(3, or), prm.phPatt(3),  prm.phCCEqPh(3, or), prm.phCCEqPhRef(3, or), prm.phCCEqPhFilt(3, or), prm.phCC(3, or), prm.phCCRef(3, or), prm.phCCFilt(3, or),...
        fro(1), froPatt(3),  froCCEqPh(1), froCCEqPhRef(1), froCCEqPhFilt(1), froCC(1), froCCRef(1), froCCFilt(1),...
        fro(2), froPatt(2),  froCCEqPh(2), froCCEqPhRef(2), froCCEqPhFilt(2), froCC(2), froCCRef(2), froCCFilt(2),...
        fro(3), froPatt(3),  froCCEqPh(3), froCCEqPhRef(3), froCCEqPhFilt(3), froCC(3), froCCRef(3), froCCFilt(3),...
        'VariableNames', VarNames);
    if exist('Results.mat', 'file')
        T = [load('Results.mat').T; T];
    end
    save('Results.mat','T');
    
end
end
