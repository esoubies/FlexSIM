function fig_id=DisplayReconstruction(rec,wf,fig_id)

if ~isgraphics(fig_id), fig_id=figure; end
set(0,'CurrentFigure',fig_id);
s1=subplot(1,2,1);imdisp(rec,'SIM Reconstruction',0);
s2=subplot(1,2,2);imdisp(wf,'Widefield',0);
Link = linkprop([s1, s2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);drawnow;
end