function fig_id=DisplayPattParams(y,params,k,phase,a,fig_id,id_patch)
%% Pre-computations
sz_y=size(y); first=0;
if ~isgraphics(fig_id)
    fig_id=figure; first=1;
    fig_id.Position=[500 500 600 800];
end
set(0,'CurrentFigure',fig_id);
%% Image with wavevectors super-imposed on the FT
if id_patch>0
    if first || (id_patch > length(fig_id.Children.Children))
        uitab( 'Title', ['Patch #',num2str(id_patch)]);
    else
        delete(fig_id.Children.Children(id_patch));
        uitab( 'Title', ['Patch #',num2str(id_patch)]);
        perm=1:length(fig_id.Children.Children);perm=[perm(1:max(id_patch-1,1)),perm(end),perm(id_patch:end-1)];
        fig_id.Children.Children=fig_id.Children.Children(perm);
    end
    ax = axes(fig_id.Children.Children(id_patch),'Units','pixels','Position',[60 250 500 500]);
    cla(ax, 'reset');
else
    ax = axes(fig_id,'Units','pixels','Position',[60 250 500 500]);
end
y = (y - min(y(:))) / (max(y(:)) - min(y(:))); 
G = RemoveWFandMask(y,mean(y,3),params);
imagesc(log10(sum(abs(fftshift(fft2(G))),3)+1), 'Parent', ax);colormap(ax,viridis);
axis(ax,'equal','off');hold(ax,'on');
fc = 2*params.Na/params.lamb*params.res;
drawellipse(ax,'Center',sz_y(1:2)/2,'SemiAxes',fc.*sz_y(2:-1:1),'StripeColor','w','InteractionsAllowed','none');
for i = 1:params.nbOr
    tmp = k(i, :) * sz_y(1) * params.res / pi + sz_y(1)/2+1;
    plot(ax,tmp(1), tmp(2), 'ro', 'MarkerSize', 8, 'LineWidth', 3);
end
text(sz_y(2)/2,-sz_y(1)*0.04,'\bf Pre-processed Data (FT) + OTF support + Detected wavevectors','HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',12);

%% Table with estimated parameters
if params.method==0
    % TODO
elseif params.method==1
    % TODO
elseif params.method==2
    patternParams = horzcat(k, phase,phase + pi/3,phase + 2*pi/3, a);       % Initialize the data
    for ii=0:params.nbOr-1
        text(sz_y(2)*0.06,sz_y(1)*(1+0.2+ii*0.035),['\bf Orr #',num2str(ii+1)] ,'HorizontalAlignment' ,'left','VerticalAlignment', 'top');
    end
end
text(sz_y(2)/2,sz_y(1)*1.1,'\bf Estimated parameters','HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',14);
format short;
TString = evalc('disp(patternParams)');
text(sz_y(2)*0.54,sz_y(1)*1.2,TString,'HorizontalAlignment' ,'center','VerticalAlignment', 'top');
text(sz_y(2)*0.2,sz_y(1)*1.16,'\bf Kx         Ky        Ph #1    Ph #2    Ph #3     Amp' ,'HorizontalAlignment' ,'left','VerticalAlignment', 'top');

drawnow;
end