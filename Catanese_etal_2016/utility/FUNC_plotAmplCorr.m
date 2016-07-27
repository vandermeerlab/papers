function stats = FUNC_plotAmplCorr(cfg_in,shift_values,lg_mat,hg_mat)
% function stats = FUNC_plotAmplCorr(cfg_in,shift_values,lg_mat,hg_mat)
% 
% plotting function for amplitude correlations

cfg_def = [];
cfg_def.function_handle = gcf;
cfg_def.title_string = [];
cfg_def.fontsize = 14;
cfg_def.linewidth = 1;
cfg_def.caxis = [0.5 1];
cfg_def.cmap = 'hot';

cfg = ProcessConfig(cfg_def,cfg_in);

shift_values = shift_values * 1000; % convert to ms

%%
figure(cfg.function_handle)
subplot(121)
[m_val,m_idx_lg] = max(lg_mat,[],2);
[~,s_idx] = sort(m_idx_lg,'ascend');

nEvents = length(m_val);
evt_idx = 1:nEvents;

imagesc(shift_values,evt_idx,lg_mat(s_idx,:));
hold on;
plot([0 0],ylim,'k:','LineWidth',1);
plot(shift_values(m_idx_lg(s_idx)),evt_idx,'.w');

stats.p_lg = signrank(shift_values(m_idx_lg));
stats.n_lg = nEvents;

set(gca,'FontSize',cfg.fontsize,'LineWidth',cfg.linewidth);
title(sprintf('%s mean lag %1.2e (p = %1.2e)',cfg.title_string,nanmean(shift_values(m_idx_lg)),stats.p_lg),'Interpreter','none');
xlabel('time (ms)');
caxis(cfg.caxis); colormap(cfg.cmap);

%%
subplot(122)
[m_val,m_idx_hg] = max(hg_mat,[],2);
[~,s_idx] = sort(m_idx_hg,'ascend');

nEvents = length(m_val);
evt_idx = 1:nEvents;

imagesc(shift_values,evt_idx,hg_mat(s_idx,:));
hold on;
plot([0 0],ylim,'k:','LineWidth',1);
plot(shift_values(m_idx_hg(s_idx)),evt_idx,'.w');

stats.p_hg = signrank(shift_values(m_idx_hg));
stats.n_hg = nEvents;

set(gca,'FontSize',cfg.fontsize,'LineWidth',cfg.linewidth);
title(sprintf('%s mean lag %1.2e (p = %1.2e)',cfg.title_string,nanmean(shift_values(m_idx_hg)),stats.p_hg),'Interpreter','none');
xlabel('time (ms)');
caxis(cfg.caxis); colormap(cfg.cmap);

%% difference
stats.p_lgVShg = ranksum(shift_values(m_idx_lg),shift_values(m_idx_hg));

%%
stats.lgp = shift_values(m_idx_lg);
stats.hgp = shift_values(m_idx_hg);