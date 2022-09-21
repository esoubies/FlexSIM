function fig_id=DisplayStack(stack,titl,fig_id)
if ~isgraphics(fig_id), fig_id=figure; end
set(0,'CurrentFigure',fig_id);
sliceViewer(stack);
title(titl);
drawnow;set(gcf,'Visible','on');
end