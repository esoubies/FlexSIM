function [fig_id,tab]=DisplayPattParams(y,params,k,phase,a,fig_id,id_patch,tab)
%% Pre-computations
sz_y=size(y);
if ~isgraphics(fig_id), fig_id=figure; end
figure(fig_id);
%% Image with wavevectors super-imposed on the FT
fig_id.Position=[500 500 600 800];
if id_patch>0
    if id_patch > length(tab)
        tab{id_patch}=uitab( 'Title', ['Patch #',num2str(id_patch)]);
    end
    ax = axes(tab{id_patch},'Units','pixels','Position',[60 250 500 500]);
    cla(ax, 'reset');
else
    ax = axes(fig_id,'Units','pixels','Position',[60 250 500 500]);
end
imagesc(log10(sum(abs(fftshift(fft2(y))),3)+1), 'Parent', ax);colormap(ax,viridis);
axis(ax,'equal','off');hold(ax,'on');
fc = 2*params.Na/params.lamb*params.res;
drawellipse(ax,'Center',sz_y(1:2)/2,'SemiAxes',fc.*sz_y(2:-1:1),'StripeColor','w','InteractionsAllowed','none');
for i = 1:params.nbOr
    tmp = k(i, :) * sz_y(1) * params.res / pi + sz_y(1)/2+1;
    plot(ax,tmp(1), tmp(2), 'ro', 'MarkerSize', 8, 'LineWidth', 3);
end
text(sz_y(2)/2,-sz_y(1)*0.04,'\bf OTF support + Detected wavevectors','HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',14);

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