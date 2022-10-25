function fig_id=DisplayStack(stack,titl,fig_id)
if ~isgraphics(fig_id), fig_id=figure; end
set(0,'CurrentFigure',fig_id);
if exist('sliceViewer')
    sliceViewer(stack);
else
    imagesc(stack(:,:,1));axis image; axis off; colormap gray;
end
title(titl);
drawnow;set(gcf,'Visible','on');
end