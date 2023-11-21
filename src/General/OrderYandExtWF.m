function [y, wf]  = OrderYandExtWF(y, params)
%--------------------------------------------------------------------------
% Function y = OrderYandExtWF(y, params)
% 
% Checks the acquisition convention of the SIM raw images, and reorders
% stack to follow the convention `ap`. Extract or compute the WF.
%
% Inputs : y       -> Raw SIM data
%
%          params  -> Structures with acquisition fields (see EstimatePatterns 
%                     for details). 
%
% Outputs: y       -> Reordered ("ap") SIM data
%          wf      -> Widefield image if provided, otherwise compute it as
%                     the mean over phases for each orientations
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

sz=size(y);
if isempty(strfind(params.StackOrder,"w"))
    nimgs=params.nbOr*params.nbPh;    
    nt=sz(3)/nimgs; 
    assert(mod(nt,1)==0,'Number of raw images is not a multiple of nbOr*nbPh');
else
    nimgs=params.nbOr*params.nbPh+1;
    nt=sz(3)/nimgs; 
    assert(mod(nt,1)==0,'Number of raw images is not a multiple of nbOr*nbPh + 1 (the +1 comes from the fact that params.StackOrder contains w)');
end
y=reshape(y,[sz(1:2),nimgs,nt]);
% Check the acquisition convention of the user and convert to ap(w)
switch string(params.StackOrder)
    case "pa" 
        wf=[];
        % Do nothing
    case "ap" 
        wf=[];
        % Reorder stack in angle-phase mode - reshape(reshape(1:9, [3, 3])', 1, [])
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []);
        y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
    case "apw"
        wf = y(:,:,end,:); 
        y = y(:,:,1:end-1,:); 
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []); 
        y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
    case "paw"
        wf = y(:,:,end,:);
        y = y(:,:,1:end-1,:);
    case "wpa"
        wf = y(:,:,1,:);
        y = y(:,:,2:end,:);
    case "wap"
        wf = y(:,:,1,:);
        y = y(:,:,2:end,:);
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []);
        y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
end
end