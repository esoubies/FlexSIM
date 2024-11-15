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
if strcmp(params.StackOrder,'axp')
    assert(mod(sz(1)/params.nbOr,1)==0 && mod(sz(2)/params.nbPh,1)==0,'(x,y) size of raw data should be multiple of nbOr and nbPh  when StackOrder is axp or axp');
    nimgs=params.nbOr*params.nbPh;
    if length(sz)==2, nt=1; else nt=sz(3); end
    sz1=sz(1)/params.nbOr;sz2=sz(2)/params.nbPh;
    tmp=zeros([sz1,sz2,nimgs,nt]);
    for it=1:nt
        idx=1;
        for ior=1:params.nbOr
            for iph=1:params.nbPh
                tmp(:,:,idx,it)=y((ior-1)*sz1+1:ior*sz1,(iph-1)*sz2+1:iph*sz2,it);
                idx=idx+1;
            end
        end
    end
    wf=[];y=tmp;
elseif strcmp(params.StackOrder,'pxa')
    assert(mod(sz(1)/params.nbPh,1)==0 && mod(sz(2)/params.nbOr,1)==0,'(x,y) size of raw data should be multiple of nbPh and nbOr when StackOrder is axp or pxa');
    nimgs=params.nbOr*params.nbPh;
    if length(sz)==2, nt=1; else nt=sz(3); end
    sz1=sz(1)/params.nbPh;sz2=sz(2)/params.nbOr;
    tmp=zeros([sz1,sz2,nimgs,nt]);
    for it=1:nt
        idx=1;
        for ior=1:params.nbOr
            for iph=1:params.nbPh
                tmp(:,:,idx,it)=y((iph-1)*sz1+1:iph*sz1,(ior-1)*sz2+1:ior*sz2,it);
                idx=idx+1;
            end
        end
    end
    wf=[];y=tmp;
else
    if isempty(strfind(params.StackOrder,'w'))
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
    switch params.StackOrder
        case 'pa'
            wf=[];
            % Do nothing
        case 'ap'
            wf=[];
            % Reorder stack in angle-phase mode - reshape(reshape(1:9, [3, 3])', 1, [])
            newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []);
            y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
        case 'apw'
            wf = y(:,:,end,:);
            y = y(:,:,1:end-1,:);
            newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []);
            y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
        case 'paw'
            wf = y(:,:,end,:);
            y = y(:,:,1:end-1,:);
        case 'wpa'
            wf = y(:,:,1,:);
            y = y(:,:,2:end,:);
        case 'wap'
            wf = y(:,:,1,:);
            y = y(:,:,2:end,:);
            newOrder = reshape(reshape(1:params.nbOr*params.nbPh, [params.nbOr, params.nbPh])', 1, []);
            y(:,:,1:params.nbOr*params.nbPh,:) = y(:,:,newOrder,:);
    end
end
end