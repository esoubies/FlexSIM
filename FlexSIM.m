function res = FlexSIM(params)
%--------------------------------------------------------------------------
% Function res = FlexSIM(params)
% 
% FlexSIM [1] should be called from a script defining the structure `params`,
% which  includes all the paths and optical parameters necessary for a 
% reconstruction. (See folder Examples)
%
% Inputs : params -> Structure containing all the necessary data (optical
%                    and reconstruction parameters, paths, etc.)  
%
% Outputs: res    -> Structure with the final image (in the field `res`)
%                    and other intermediate results, like patterns.  
%
% [1] FlexSIM: ADD REF TO PAPER
%
% See also EstimatePatterns.m and Reconstruct.m
%
% Copyright (2022) A. Nogueron (anogueron.1996@gmail.com)
%                  E. Soubies (emmanuel.soubies@irit.fr) 
%--------------------------------------------------------------------------

%% Data loading + routinary checks
y = double(loadtiff(params.DataPath));    % Read data 
sz_y=size(y);                             % Raw SIM data size
CheckParams(params);                      % Check conformity of parameters

% - Displays
if params.displ > 0
    figure;sliceViewer(y);title('SIM Raw data');
    drawnow;set(gcf,'Visible','on');
end

%% Pattern Estimation
disp('=== Patterns parameter estimation START ...');
% - Estimate pattern parameters
[k, phase, a] = EstimatePatterns(params, y);
disp('=== Patterns parameter estimation END.');
% - Generate Patterns for reconstruction
a=a./a; % Hardcode to 1 for now (to be as in previous version) 
[patterns] = GenerateReconstructionPatterns(params,k,phase,a,sz_y); 

% - Displays
if params.displ > 0          
    % - Display of estimated patterns
    figure;sliceViewer(patterns);  title('Estimated Patterns');
    % - Displays related to estimated parameters
    % Image with wavevectors super-imposed on the FT
    fig=figure('Position',[500 500 600 800]);
    ax = axes(fig,'Units','pixels','Position',[60 250 500 500]);   
    imagesc(log10(sum(abs(fftshift(fft2(y))),3)+1), 'Parent', ax);colormap(ax,viridis);
    axis(ax,'equal','off');hold(ax,'on');
    fc = 2*params.Na/params.lamb*params.res;
    drawellipse(ax,'Center',sz_y(1:2)/2,'SemiAxes',fc.*sz_y(2:-1:1),'StripeColor','w','InteractionsAllowed','none');
    for i = 1:params.nbOr
        tmp = k(i, :) * sz_y(1) * params.res / pi + sz_y(1)/2+1;
        plot(ax,tmp(1), tmp(2), 'ro', 'MarkerSize', 8, 'LineWidth', 3);
    end  
    text(sz_y(2)/2,-15,'\bf OTF support + Detected wavevectors','HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',14);
    
    % Table with estimated parameters
    if params.method==0
        % TODO
    elseif params.method==1
        % TODO
    elseif params.method==2        
        patternParams = horzcat(k, phase,phase + pi/3,phase + 2*pi/3, a);       % Initialize the data
        for ii=0:params.nbOr-1
            text(16,300+ii*9,['\bf Orr #',num2str(ii+1)] ,'HorizontalAlignment' ,'left','VerticalAlignment', 'top');
        end
    end
    text(sz_y(2)/2,275,'\bf Estimated parameters','HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',14);
    format short;
    TString = evalc('disp(patternParams)');
    text(sz_y(2)/2+10,300,TString,'HorizontalAlignment' ,'center','VerticalAlignment', 'top');
    text(55,290,'\bf Kx         Ky       Ph #1    Ph #2    Ph #3     Amp' ,'HorizontalAlignment' ,'left','VerticalAlignment', 'top');
end
drawnow;

%% Reconstruction
rec = Reconstruct(y,patterns,params);

% - Displays
if params.displ > 0
    figure;
    s1=subplot(1,2,1);imdisp(rec,'SIM Reconstruction',0);
    s2=subplot(1,2,2);imdisp(imresize(mean(y,3),size(rec)),'Widefield',0);
    Link = linkprop([s1, s2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
end

%% Save
% - Save reconstruction / patterns / reconst. parameters / pattern parameters
if params.sav
    prefix=params.DataPath(1:end-4);
    saveastiff(single(rec),strcat(prefix,'_Rec.tif'));
    saveastiff(single(patterns),strcat(prefix,'_Patt.tif'));
    save(strcat(prefix,'_Params'),'params');
end

% - Fill output variable
res.rec=rec;
res.patt=patterns;

end

