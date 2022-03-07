function FUNC_plotPSI(h,this_session_data)

axes(h);

fvec = this_session_data.psi.freq;
psl = this_session_data.psi.psispctrm;
psl_sem = this_session_data.psi.psispctrmsem;
icoh = this_session_data.coh.cohspctrm;

lgx = [45 45 65 65]; lgc = [0.8 1 0.8];
hgx = [70 70 90 90]; hgc = [1 0.8 0.8];
fa = 1;

[ax,h1,h2] = plotyy(fvec,icoh,fvec,psl);
set(ax(1),'YLim',[-20 20])
set(ax(1),'YTick',-20:10:20)
set(ax(1),'YTickLabel',{'-20','','0','','20'})
set(ax(2),'YLim',[-1 1])
set(ax(2),'YTick',[-1 -0.5 0 0.5 1])
set(ax(2),'YTickLabel',{'-1','','0','','1'})
set(ax(1),'LineWidth',1,'XLim',[40 100],'XTick',40:20:100,'FontSize',18);
set(ax(2),'LineWidth',1,'XLim',[40 100],'XTick',40:20:100,'FontSize',18);
set(h1,'LineWidth',2);
set(h2,'LineWidth',2);
axes(ax(1)); box off; ylabel('phase lag');
axes(ax(2)); box off; hold on; plot([35 105],[0 0],'k--','LineWidth',1); ylabel('phase slope');

yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))