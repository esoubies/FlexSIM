function [y, wf]  = OrderY(y, params)
%--------------------------------------------------------------------------
% Function y = OrderY(y, params)
% 
% Checks the acquisition convention of the SIM raw images, and reorders
% stack to follow the convention `ap`. 
%
% Inputs : y       -> Raw SIM data
%
%          params  -> Structures with acquisition fields (see EstimatePatterns 
%                     for details). 
%
% Outputs: y       -> Reordered ("ap") SIM data
%          wf      -> Widefield image if provided, otherwise 0
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------
% Check the acquisition convention of the user and convert to ap(w)
switch string(params.StackOrder)
    case "pa" 
        wf = 0; 
        % Do nothing
    case "ap" 
        % Reorder stack in angle-phase mode - reshape(reshape(1:9, [3, 3])', 1, [])
        wf = 0;
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []); 
        y(:,:,1:params.nbOr*params.nbPh) = y(:,:,newOrder); 
    case "apw"
        wf = y(:,:,end); 
        y = y(:,:,1:end-1); 
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []); 
        y(:,:,1:params.nbOr*params.nbPh) = y(:,:,newOrder);
    case "paw"
        wf = y(:,:,end); 
        y = y(:,:,1:end-1); 
    case "wpa"
        wf = y(:,:,1); 
        y = y(:,:,2:end); 
    case "wap"
        wf = y(:,:,1); 
        y = y(:,:,2:end); 
        newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []); 
        y(:,:,1:params.nbOr*params.nbPh) = y(:,:,newOrder); 
end
end