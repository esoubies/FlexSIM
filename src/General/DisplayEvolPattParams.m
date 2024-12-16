function fig_id=DisplayEvolPattParams(params,k_est,phase_est,k,phase,fig_id)
%% Pre-computations
if ~isgraphics(fig_id)
    fig_id=figure; 
    fig_id.Position=[500 500 1200 900];
end
set(0,'CurrentFigure',fig_id);

% Plot evolution of wavevector norm
if ~isempty(params.framePattEsti), fr = params.framePattEsti; else fr = 1:params.nframes; end
for ii=1:params.nbOr
    subplot(3,3,ii);hold on;
    plot(fr,squeeze(sqrt((k_est(ii,1,:)/(2*params.Na)*params.lamb/pi).^2 + (k_est(ii,2,:)/(2*params.Na)*params.lamb/pi).^2)),'k.','markersize',15);
    if params.cstTimePatt || params.rollMed>0, plot(squeeze(sqrt((k(ii,1,:)/(2*params.Na)*params.lamb/pi).^2 + (k(ii,2,:)/(2*params.Na)*params.lamb/pi).^2)),'rx','markersize',8,'linewidth',2);end
    title(['Orientation #',num2str(ii)]);grid; set(gca,'fontsize',14);
    xticks(1:max(round(params.nframes/5),1):params.nframes);axis([1 params.nframes params.ringRegionSearch(1) params.ringRegionSearch(2)]);
    if ii==1, ylabel('||k||/f_c'); end
end
% Plot wavevector angle
agl_min=Inf;agl_max=-Inf;
for ii=1:params.nbOr
    subplot(3,3,3+ii); hold on;
    agl = squeeze(atan(k_est(ii,2,:)./k_est(ii,1,:)));
    agl_min = min(agl_min,min(agl));agl_max=max(agl_max,max(agl));
    plot(fr,agl,'k.','markersize',15);
    if params.cstTimePatt || params.rollMed>0, plot(squeeze(atan(k(ii,2,:)./k(ii,1,:))),'rx','markersize',8,'linewidth',2); end
    grid; set(gca,'fontsize',14); %axis([1 params.nframes -pi/2 pi/2]);
    xticks(1:max(round(params.nframes/5),1):params.nframes);
    if ii==1, ylabel('angle(k)'); end
end
for ii=1:params.nbOr
    subplot(3,3,3+ii);
    axis([1 params.nframes agl_min-0.1 agl_max+0.1]);
end
% Plot phase offset
for ii=1:params.nbOr
    subplot(3,3,6+ii); hold on;
    plot(fr,squeeze(phase_est(ii,1,:)),'k.','markersize',15);
    if params.cstTimePatt || params.rollMed>0, plot(squeeze(phase(ii,1,:)),'rx','markersize',8,'linewidth',2); end
    xlabel('Frames');grid; set(gca,'fontsize',14);
    xticks(1:max(round(params.nframes/5),1):params.nframes); % axis([1 params.nframes 0 pi]);
    if ii==1, ylabel('Phase offset'); end
end
if params.cstTimePatt, Lgnd = legend('Estimated','Median (all frames)');
elseif params.rollMed>0, Lgnd = legend('Estimated',['Rolling Median (size ',num2str(params.rollMed),')']);
else  Lgnd = legend('Estimated','Mean (all frames)'); end
legend('NumColumns',2);
Lgnd.Position(1) = 0.35;
Lgnd.Position(2) = 0.01;
end