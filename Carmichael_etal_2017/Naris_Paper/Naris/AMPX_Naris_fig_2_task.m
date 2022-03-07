function AMPX_Naris_fig_2_task(all_data_pre, all_data_task, varargin)
%% AMPX_Naris_fig_2_task: generates figure 2 which contains the heatmaps
% for the average rewarded/approach gamma 50 and gamma 80 events for the
% two rats that reached task criterion
%
%
% EC - 2016-05-31

%% extract inputs
extract_varargin

global PARAMS;

%% load the data for the task events
% load('C:\temp\Naris_all_data_task')

%% R049
h1 = figure(21);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);

Session_list = fieldnames(all_data_task);
bands = {'lg_reward', 'hg_reward', 'lg_approach', 'hg_approach'};
for iSess = 1:length(Session_list)
    for  iBand = 1:length(bands)
        if isfield(all_data_task.(Session_list{iSess}).(bands{iBand}).power, 'power_avg') || isempty(all_data_task.(Session_list{iSess}).(bands{iBand}).power.power_distrib_avg)
            all_data_task.(Session_list{iSess}).(bands{iBand}).power.power_distrib_avg = NaN*zeros(8,8);
            warning(['Session ' Session_list{iSess} '  ' bands{iBand} '  did not contain any events'])
        end
    end
end
subplot(2,4,1)
R049.low.approach = nanmean(cat(3,all_data_task.R049_2014_02_07.lg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_08.lg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_10.lg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.low.approach/min(min(R049.low.approach)))

subplot(2,4,2)
R049.low.reward = nanmean(cat(3,all_data_task.R049_2014_02_07.lg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_08.lg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_10.lg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.low.reward/min(min(R049.low.reward)))

subplot(2,4,5)
R049.high.approach = nanmean(cat(3,all_data_task.R049_2014_02_07.hg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_08.hg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_10.hg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.high.approach/min(min(R049.high.approach)))

subplot(2,4,6)
R049.high.reward = nanmean(cat(3,all_data_task.R049_2014_02_07.hg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_08.hg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_10.hg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.high.reward/min(min(R049.high.reward)))

% R045  (R4)
figure(21)
subplot(2,4,3)
R045.low.approach = nanmean(cat(3,all_data_task.R045_2014_04_16.lg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_17.lg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_15.lg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.low.approach/min(min(R045.low.approach)))

subplot(2,4,4)
R045.low.reward = nanmean(cat(3,all_data_task.R045_2014_04_16.lg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_17.lg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_15.lg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.low.reward/min(min(R045.low.reward)))

subplot(2,4,7)
R045.high.approach = nanmean(cat(3,all_data_task.R045_2014_04_16.hg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_17.hg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_15.hg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.high.approach/min(min(R045.high.approach)))

subplot(2,4,8)
R045.high.reward = nanmean(cat(3,all_data_task.R045_2014_04_16.hg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_17.hg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_15.hg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.high.reward/min(min(R045.high.reward)))

axesHandles = findobj(get(h1,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
for isub = 1:8
    subplot(2,4,isub)
    set(gca, 'xtick', [], 'ytick', []);
    og_size = get(gca, 'position');
    colorbar('location', 'eastoutside');
    set(gca, 'position', og_size);
end

%% export part 1
if exist('save_fig')
    saveas(h1, cat(2, PARAMS.figure_dir,'\Fig2_A'), 'epsc');
    saveas(h1, cat(2, PARAMS.figure_dir,'\Fig2_A'), 'fig')
end
%% get the counts across all pre and post recording sessions
clear all_data_task
%
% load('C:\temp\Naris_all_data_pre.mat')
% load('C:\temp\Naris_all_data_post.mat')

%% get R^2 across sessions
sess_list = fieldnames(all_data_pre);
bands = {'lg','lg_ran', 'hg', 'hg_ran'};% 'hg'};
all.lg.rsq = []; all.hg.rsq = [];all.lg_ran.rsq = []; all.hg_ran.rsq = [];
all.lg.rsq2 = []; all.hg.rsq2 = [];all.lg_ran.rsq2 = []; all.hg_ran.rsq2 = [];

for iSess = 1:length(sess_list)
    for iband = 1:length(bands)
        if strcmp(sess_list{iSess}(1:4), 'R045') %|| strcmp(sess_list{iSess}(1:4), 'R049')

        else
            all.(bands{iband}).rsq = [all.(bands{iband}).rsq, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq];
            all.(bands{iband}).rsq2 = [all.(bands{iband}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq2];
        end
        
    end
end
%%
[~, plane_stats.low] = AMPX_plane_count_hist(all.lg, all.lg_ran, 'low');
hax1 = gca;
SetFigure([], gcf)
if exist('save_fig')
    saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_B'), 'epsc')
    saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_B'), 'fig')
end
hf1 = figure(22);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);
s1 = subplot(5,2,[1 3]);
pos=get(s1,'Position');
delete(s1);
hax2=copyobj(hax1,hf1);
set(hax2, 'Position', pos);
% ylim([0 100])


[~, plane_stats.high] = AMPX_plane_count_hist(all.hg, all.hg_ran, 'high');
hax11 = gca;
SetFigure([], gcf)
if exist('save_fig')
    saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_C'), 'epsc')
    saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_C'), 'fig')
end
hf1 = figure(22);
s2 = subplot(5,2,[2 4]);
pos=get(s2,'Position');
delete(s2);
hax2=copyobj(hax11,hf1);
set(hax2, 'Position', pos);

axesHandles = findobj(get(hf1,'Children'), 'flat','Type','axes');
% set(axesHandles, 'xtick', 0:20:100, 'ytick', 0:25:100, 'xlim', [0 100]);
og_size = get(axesHandles, 'position');
% set(axesHandles, 'position', og_size);
set(axesHandles, 'fontsize', 18, 'fontname', 'helvetica')

SetFigure([], gcf)


maximize
%% Same plot but with each subject as a line
for iBand = 1:length(bands)  % initialize the R^2 for each rat
    R061.(bands{iBand}).rsq = [];
    R061.(bands{iBand}).rsq2 = [];
    R054.(bands{iBand}).rsq = [];
    R054.(bands{iBand}).rsq2 = [];
    R045.(bands{iBand}).rsq = [];
    R045.(bands{iBand}).rsq2 = [];
    R049.(bands{iBand}).rsq = [];
    R049.(bands{iBand}).rsq2 = [];
end
% collect all the R^2 values for each rat.
for iSess = 1:length(sess_list)
    for iBand = 1:length(bands)
        if strcmp(sess_list{iSess}(1:4), 'R061')
            R061.(bands{iBand}).rsq = [R061.(bands{iBand}).rsq , all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq];
            R061.(bands{iBand}).rsq2 = [R061.(bands{iBand}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq2];
        end
        if strcmp(sess_list{iSess}(1:4), 'R054')
            R054.(bands{iBand}).rsq = [R054.(bands{iBand}).rsq , all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq];
            R054.(bands{iBand}).rsq2 = [R054.(bands{iBand}).rsq2 , all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq2];
        end
        if strcmp(sess_list{iSess}(1:4), 'R049')
            R049.(bands{iBand}).rsq = [R049.(bands{iBand}).rsq, all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq];
            R049.(bands{iBand}).rsq2 = [R049.(bands{iBand}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq2];
        end
        if strcmp(sess_list{iSess}(1:4), 'R045')
            R045.(bands{iBand}).rsq = [R045.(bands{iBand}).rsq, all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq];
            R045.(bands{iBand}).rsq2 = [R045.(bands{iBand}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iBand}).power.stats.rsq2];
        end
    end
end

all_rats.R061 = R061;
all_rats.R054 = R054;
all_rats.R045 = R045;
all_rats.R049 = R049;

%% Plot the R^2 histograms on a rat by rat bassis
cfg = [];
cfg.bins = 5;
cfg.c_ord = linspecer(4);
Rats = {'R061', 'R054', 'R045', 'R049'};

figure(22)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);
for ii =  1:2
    subplot(2,2,ii+2) % low gamma
    hold on
    for iRat = 1:length(Rats)
        if ii ==1
            a = all_rats.(Rats{iRat}).lg.rsq2;
            b = all_rats.(Rats{iRat}).lg_ran.rsq2;
        else
            a = all_rats.(Rats{iRat}).hg.rsq2;
            b = all_rats.(Rats{iRat}).hg_ran.rsq2;
        end
        
        [h, p, ci, stats] = ttest2(a, b);
        [n1, p1] = hist(b, 50);
        [n2, p2] = hist(a, 50);
        N1 = conv(n1,ones(1,cfg.bins),'same')/cfg.bins; %smooth over bins
        H = line(p1, N1);
        set(H, 'color', cfg.c_ord(iRat,:), 'linewidth', 2, 'LineStyle','--');
        
        N2 = conv(n2,ones(1,cfg.bins),'same')/cfg.bins; %smooth over bins
        H2 = line(p2, N2);
        set(H2, 'color', cfg.c_ord(iRat,:), 'linewidth', 2);
        %         legend([h h2], {'Control', 'Gamma'})
        xlim([0 100])
    end
    SetFigure([], gcf)

end
maximize

%% export part 1
if exist('save_fig')
    saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_B_C_new'), 'epsc')
    % saveas(gcf, cat(2, PARAMS.figure_dir,'\Fig2_B_C_new'), 'fig') % fails
    % with obscure error -- MvdM
end