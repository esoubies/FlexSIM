function fig_id=DisplayStack(stack,titl,fig_id)
if ~isgraphics(fig_id), fig_id=figure; end
sliceViewer(stack,'Parent',fig_id);
title(titl);
drawnow;set(gcf,'Visible','on');
end