%% compute cross-validated decoding error for all sessions
% uses code from van der Meer et al. (2017) Hippocampus
% select cfg_master.TCsmooth = 1; and cfg_master.QsmoothSD = 0.002;
% choose cfg_master.dt = 0.025;
GENERATE_decErr_paramSweep; % by default, doesn't equalize number of laps for L and R

%% compute variables to plot

fw = [2 1 2 1 2 1 2 1 2 1 2 1 1 2 1 2 1 2 1 2 1 2 1 2]; % 1 food, 2 water
lerr = ALL_decErr.left.meanErr;
rerr = ALL_decErr.right.meanErr;
skip_idx = [1 7 8 9 12];

fw(skip_idx) = [];
lerr(skip_idx) = [];
rerr(skip_idx) = [];

flr = lerr(fw == 1); frr = rerr(fw == 1); % food, raw errors
wlr = lerr(fw == 2); wrr = rerr(fw == 2);

flp = lerr(fw == 1)./(lerr(fw == 1) + rerr(fw == 1));
frp = rerr(fw == 1)./(lerr(fw == 1) + rerr(fw == 1));
wlp = lerr(fw == 2)./(lerr(fw == 2) + rerr(fw == 2));
wrp = rerr(fw == 2)./(lerr(fw == 2) + rerr(fw == 2));

%% plot
subplot(241); % raw
x = [1 2 4 5]; xl = [0 6]; fs = 10;

bar(x(1:2),[nanmean(flr) nanmean(frr)],'FaceColor',[1 0 0],'EdgeColor','none'); hold on;
bar(x(3:4),[nanmean(wlr) nanmean(wrr)],'FaceColor',[0 0 1],'EdgeColor','none');

plot(x(1:2),[flr frr],'.','Color',[0.7 0 0]);
plot(x(3:4),[wlr wrr],'.','Color',[0 0 0.7]);

set(gca,'XLim',xl,'XTick',x,'XTickLabel',{'L','R','L','R'},'LineWidth',1,'FontSize',fs,'YTick',0:20:60); box off;
%set(gca,'XLim',xl,'XTick',x,'XTickLabel',{'L','R','L','R'},'LineWidth',1,'FontSize',fs,'YTick',0:10:30); box off;
ylabel('cross-validated decoding error (cm)');

subplot(242); % proportional
bar(x(1:2),[nanmean(flp) nanmean(frp)],'FaceColor',[1 0 0],'EdgeColor','none'); hold on;
bar(x(3:4),[nanmean(wlp) nanmean(wrp)],'FaceColor',[0 0 1],'EdgeColor','none');

plot(x(1:2),[flp frp],'.','Color',[0.7 0 0]);
plot(x(3:4),[wlp wrp],'.','Color',[0 0 0.7]);

set(gca,'XLim',xl,'XTick',x,'XTickLabel',{'L','R','L','R'},'LineWidth',1,'FontSize',fs,'YTick',0:0.5:1,'YLim',[0 1]); box off;
ylabel('cross-validated decoding error (p)');

subplot(243); % null
bar(x(1:2),[1-nanmean(flp) 1-nanmean(frp)],'FaceColor',[1 0 0],'EdgeColor','none'); hold on;
bar(x(3:4),[1-nanmean(wlp) 1-nanmean(wrp)],'FaceColor',[0 0 1],'EdgeColor','none');
set(gca,'XLim',xl,'XTick',x,'XTickLabel',{'L','R','L','R'},'LineWidth',1,'FontSize',fs,'YTick',0:0.5:1,'YLim',[0 1]); box off;
ylabel('SWR content null (p)');


%% xport
cfg.output_fn = 'xval_decAcc';
maximize;
print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-painters','-dpdf','-r300',cfg.output_fn);
print(gcf,'-painters','-depsc','-r300',cfg.output_fn);

