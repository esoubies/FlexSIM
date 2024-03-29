function fig_id=DisplayPattParams(y,params,k,phase,fig_id,id_patch,titl)
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
axis(ax,'square','on');hold(ax,'all');axis off;
fc = 2*params.Na/params.lamb*params.res;

h=[];leg={};
cen=floor(sz_y(2:-1:1)/2)+1;
[xp,yp]=circle(fc.*sz_y(2:-1:1));h(end+1)=plot(cen(1)+xp,cen(2)+yp,'Linewidth',1.5,'color','k');  % fc
[xpin,ypin]=circle(params.ringRegionSearch(1)*fc.*sz_y(2:-1:1));plot(cen(1)+xpin,cen(2)+ypin,'Linewidth',2,'color',[0 0.4470 0.7410],'LineStyle','--');  % ring in
[xpout,ypout]=circle(params.ringRegionSearch(2)*fc.*sz_y(2:-1:1));h(end+1)=plot(cen(1)+xpout,cen(2)+ypout,'Linewidth',2,'color',[0 0.4470 0.7410],'LineStyle','--');  % ring out
fill([xpin flip(xpout)]+cen(1),[ypin flip(ypout)]+cen(2),[0 0.4470 0.7410],'FaceAlpha',0.3,'EdgeColor','none');
[xp,yp]=circle(params.maskWF*fc.*sz_y(2:-1:1));h(end+1)=plot(cen(1)+xp,cen(2)+yp,'Linewidth',2,'color','r','LineStyle','--');  % ring mask WF
fill([zeros(size(xp)) flip(xp)]+cen(1),[zeros(size(yp))  flip(yp)]+cen(2),'r','FaceAlpha',0.1,'EdgeColor','none');
%drawellipse(ax,'Center',floor(sz_y(2:-1:1)/2)+1,'SemiAxes',fc.*sz_y(2:-1:1),'StripeColor','w','InteractionsAllowed','none','FaceAlpha',0,'Color','k');
leg{end+1}='OTF cutoff';
leg{end+1}='Ring Region Search';
leg{end+1}='Mask WF';
col=[0 0.4470 0.7410
0.8500 0.3250 0.0980
0.9290 0.6940 0.1250
0.4940    0.1840    0.5560
0.4660    0.6740    0.1880
0.3010    0.7450    0.9330
0.6350    0.0780    0.1840];
for i = 1:size(k, 1)
    tmp = k(i, :) .* sz_y(2:-1:1) * params.res / pi + sz_y(2:-1:1)/2+1;
    h(end+1)=plot(ax,tmp(1), tmp(2), 'o', 'MarkerSize', 8, 'LineWidth', 3,'Color',col(i,:));
    leg{end+1}=['Or #',num2str(i)];
end

legend(h,leg,'NumColumns',2);
text(sz_y(2)/2,-sz_y(1)*0.04,['\bf',titl],'HorizontalAlignment' ,'center','VerticalAlignment', 'top','FontSize',12);

%% Table with estimated parameters
% - Find monospaced font
mono_fonts={ ...
'Liberation Mono',...
'DejaVu Sans Mono',...
'Monospaced',...
'Courier',...
'FreeMono',...
'Go Mono',...
'Noto Mono',...
'Mitra Mono'};
id_font=find(cell2mat(cellfun(@(x) sum(strcmp(x,listfonts)),mono_fonts,'UniformOutput',false)),1,'first');

% - Display Table
if params.eqPh
    patternParams=k.*sz_y(2:-1:1) * params.res / pi;
    Col_name={' Kx',' Ky'};
    if ~isempty(phase)
        for ii=1:params.nbPh
            patternParams = horzcat(patternParams,mod(phase + (ii-1)*pi/params.nbPh,pi));
            Col_name{end+1}=['Ph',num2str(ii)];
        end
    end
    for ii=1:params.nbOr
        Row_name{ii}=['Or',num2str(ii)];
    end
else
    Col_name={' Kx',' Ky'};
    if ~isempty(phase)
        patternParams = horzcat(k.*sz_y(2:-1:1) * params.res / pi, phase);       % Initialize the data
        for ii=1:params.nbPh
            Col_name{end+1}=['Ph',num2str(ii)];
        end
    else
        patternParams = horzcat(k.*sz_y(2:-1:1) * params.res / pi);       % Initialize the data
    end    
    for ii=1:params.nbOr
        Row_name{ii}=['Or',num2str(ii)];
    end
end
T=array2table(patternParams,'VariableName',Col_name,'Rowname',Row_name);
tableCell = [' ',T.Properties.VariableNames;T.Properties.RowNames, table2cell(T)]; 
tableCell(cellfun(@isnumeric,tableCell)) = cellfun(@(x) sprintf('% 1.2f',x), tableCell(cellfun(@isnumeric,tableCell)),'UniformOutput',false); 
tableChar = splitapply(@strjoin,pad(tableCell),(1:length(Row_name)+1)');
if id_patch>0
    axes(fig_id.Children.Children(id_patch),'Position',[.1,.1,.85,.18], 'Visible','off');
else
    axes(fig_id,'Position',[.1,.1,.85,.18], 'Visible','off');
end
text(.5,.95,tableChar,'VerticalAlignment','Top','HorizontalAlignment','Center','FontName',mono_fonts{id_font},'Interpreter','Tex');
text(.5,1.1,'\bf Estimated parameters','VerticalAlignment','Top','HorizontalAlignment','Center');

drawnow;
end
function [xp,yp]=circle(d)
ang=linspace(0,2*pi,100); 
xp=d(1)*cos(ang);yp=d(2)*sin(ang);
end
