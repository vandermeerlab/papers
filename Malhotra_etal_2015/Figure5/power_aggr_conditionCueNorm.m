clear all; close all;
rng('shuffle');
%addpath('C:\Users\mvdm\Dropbox\projects\Sushant\shared');

%cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
%load fd_isidro
cd('D:\My_Documents\My Dropbox\projects\Sushant\2014-09-19') % athena
load fd; % list of files to process, obtained from remove_cueviolation.m

% only certain rat
%fd = fd(end-2:end); % R014
%files = fd(1); % R020
%fd = fd(2:8); % R018
%fd = fd(9:13); % R016

%% build data set -- all trials together!
twin = [-2 4];
f_lg = [50 65];
f_hg = [70 85];

baseline = [-0.5 0];

data_all = []; data_all_cue = [];
for iFD = 1:length(fd)
    
    cd(fd{iFD});
    fprintf('\n***MAIN: Entering folder %s (%d/%d)\n\n',pwd,iFD,length(fd));
       
    run(FindFile('*keys.m'));
    fc = ExpKeys.goodGamma(1);
    data = ft_read_neuralynx_interp(fc);
    
    fprintf('\n***MAIN: Finished loading data\n\n');
    
    data.label{1} = 'vStr';
    
    % somehow restrict to TimeOnTrack for normalization
    cfg = []; cfg.t = [ExpKeys.TimeOnTrack(1)-5 ExpKeys.TimeOffTrack(2)+5];
    data = ft_restrict_data(cfg,data);
    
    % add hilberted power -- filters could be improved!
    data2 = data;
    
    % low-gamma
    data2f = ft_filterLFP(data2,f_lg,'ford',4);
    data.trial{1}(2,:) = data2f.trial{1};
    data.label{2} = 'lg_pow';
    data.trial{1}(2,:) = abs(hilbert(data.trial{1}(2,:))).^2;
    
    data.label{4} = 'lg_powz';
    data.trial{1}(4,:) = zscore(data.trial{1}(2,:));
    
    % high-gamma
    data2f = ft_filterLFP(data2,f_hg,'ford',4);
    data.trial{1}(3,:) = data2f.trial{1};
    data.label{3} = 'hg_pow';
    data.trial{1}(3,:) = abs(hilbert(data.trial{1}(3,:))).^2;
    
    data.label{5} = 'hg_powz';
    data.trial{1}(5,:) = zscore(data.trial{1}(3,:));
    
    clear data2;
    
    %% get all nosepokes
    
    all_cues = {'c1','c3','c5','lo','hi'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
    
    for iCue = 1:length(all_cues)
        
        fprintf('\n***MAIN: Obtaining cue set (%d/%d)\n\n',iCue,length(all_cues));
        
        % nosepoke data 
        cfg = [];
        cfg.trialdef.cue = all_cues(iCue); % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
        cfg.trialdef.block = 'both';
        cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
        cfg.trialdef.location = 'both';
        cfg.trialfun = 'ft_trialfun_lineartracktone2';
        cfg.trialdef.hdr = data.hdr;
        cfg.trialdef.pre = 2; cfg.trialdef.post = 4;
        
        [trl, event] = ft_trialfun_lineartracktone2(cfg); cfg.trl = trl;
        this_data = ft_redefinetrial(cfg,data);
        
        event = repmat(iCue,[length(this_data.trial) 1]);
        
        if isempty(data_all)
            data_all = this_data;
            event_all = event;
        else
            fprintf('\n***MAIN: Appending data\n\n',iCue,length(all_cues));
            data_all = ft_appenddata([],data_all,this_data);
            event_all = cat(1,event_all,event);
        end
        
        % cue data (for baseline)
        cfg.trialdef.eventtype = 'cue'; % could be 'nosepoke', 'reward', 'cue'
        [trl, event] = ft_trialfun_lineartracktone2(cfg); cfg.trl = trl;
        this_data = ft_redefinetrial(cfg,data);
        
        event = repmat(iCue,[length(this_data.trial) 1]);
        
        if isempty(data_all_cue)
            data_all_cue = this_data;
            event_all_baseline = event;
        else
            fprintf('\n***MAIN: Appending data\n\n',iCue,length(all_cues));
            data_all_cue = ft_appenddata([],data_all_cue,this_data);
            event_all_baseline = cat(1,event_all_baseline,event);
        end
        
    end
    
end

%% get various hilbert powers
% restrict data
cfg = []; cfg.latency = twin;
data_all = ft_selectdata(cfg,data_all);

trl_len = length(data_all.trial{1}(1,:));

tvec = data_all.time{1}; tvec_cue = data_all_cue.time{1};

% get data
lg_idx = find(strcmp(data_all.label,'lg_pow'));
hg_idx = find(strcmp(data_all.label,'hg_pow'));
lgz_idx = find(strcmp(data_all.label,'lg_powz'));
hgz_idx = find(strcmp(data_all.label,'hg_powz'));

% do analysis for each trial type separately
trial_types = [1 3];
for iTT = 1:length(trial_types)
    
    cur_tt_idx = find(event_all == trial_types(iTT));
    nTrials = length(cur_tt_idx);

    d = cell2mat(data_all.trial(cur_tt_idx)); % trials concatenated

    lg_pow = reshape(d(lg_idx,:),[trl_len nTrials])';
    hg_pow = reshape(d(hg_idx,:),[trl_len nTrials])';
    lg_powz = reshape(d(lgz_idx,:),[trl_len nTrials])';
    hg_powz = reshape(d(hgz_idx,:),[trl_len nTrials])';
    
    % obtain baseline from cue data
    d_cue = cell2mat(data_all_cue.trial(cur_tt_idx)); % trials concatenated

    lg_pow_cue = reshape(d_cue(lg_idx,:),[trl_len nTrials])';
    hg_pow_cue = reshape(d_cue(hg_idx,:),[trl_len nTrials])';
    lg_powz_cue = reshape(d_cue(lgz_idx,:),[trl_len nTrials])';
    hg_powz_cue = reshape(d_cue(hgz_idx,:),[trl_len nTrials])';

    % normalize
    t_idxb = nearest(tvec_cue,baseline);
    lg_baseline = nanmean(lg_pow_cue(:,t_idxb(1):t_idxb(2)),2); lg_baseline = repmat(lg_baseline,[1 trl_len]);
    hg_baseline = nanmean(hg_pow_cue(:,t_idxb(1):t_idxb(2)),2); hg_baseline = repmat(hg_baseline,[1 trl_len]);
    lgz_baseline = nanmean(lg_powz_cue(:,t_idxb(1):t_idxb(2)),2); lgz_baseline = repmat(lgz_baseline,[1 trl_len]);
    hgz_baseline = nanmean(hg_powz_cue(:,t_idxb(1):t_idxb(2)),2); hgz_baseline = repmat(hgz_baseline,[1 trl_len]);

    % average etc.
    lg_pow_m{iTT} = nanmean(lg_pow); lg_pow_sd{iTT} = nanstd(lg_pow);
    hg_pow_m{iTT} = nanmean(hg_pow); hg_pow_sd{iTT} = nanstd(hg_pow);
    lg_pow_mb{iTT} = nanmean(lg_pow./lg_baseline); lg_pow_sdb{iTT} = nanstd(lg_pow./lg_baseline);
    hg_pow_mb{iTT} = nanmean(hg_pow./hg_baseline); hg_pow_sdb{iTT} = nanstd(hg_pow./hg_baseline);
    
    lg_powz_m{iTT} = nanmean(lg_powz); lg_powz_sd{iTT} = nanstd(lg_powz);
    hg_powz_m{iTT} = nanmean(hg_powz); hg_powz_sd{iTT} = nanstd(hg_powz);
    lg_powz_mb{iTT} = nanmean(lg_powz-lgz_baseline); lg_powz_sdb{iTT} = nanstd(lg_powz-lgz_baseline);
    hg_powz_mb{iTT} = nanmean(hg_powz-hgz_baseline); hg_powz_sdb{iTT} = nanstd(hg_powz-hgz_baseline);

end

%% plot hilbert versions
stat_twin1 = [-1.25 0];
stat_twin2 = [0 1.25];
%stat_twin1 = [-0.5 0];
%stat_twin2 = [0 0.5];

for iTT = 1:length(trial_types)
    
    figure;
    
    subplot(221);
    h(1) = plot(tvec,lg_pow_m{iTT},'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(tvec,lg_pow_m{iTT}+lg_pow_sd{iTT},':','Color',[1 0 0]);
    plot(tvec,lg_pow_m{iTT}-lg_pow_sd{iTT},':','Color',[1 0 0]);
    
    h(2) = plot(tvec,hg_pow_m{iTT},'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(tvec,hg_pow_m{iTT}+hg_pow_sd{iTT},':','Color',[0 1 0]);
    plot(tvec,hg_pow_m{iTT}-hg_pow_sd{iTT},':','Color',[0 1 0]);
    title('raw');
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    
    %
    subplot(222);
    h(1) = plot(tvec,lg_pow_mb{iTT},'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(tvec,lg_pow_mb{iTT}+lg_pow_sdb{iTT},':','Color',[1 0 0]);
    plot(tvec,lg_pow_mb{iTT}-lg_pow_sdb{iTT},':','Color',[1 0 0]);
    
    h(2) = plot(tvec,hg_pow_mb{iTT},'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(tvec,hg_pow_mb{iTT}+hg_pow_sdb{iTT},':','Color',[0 1 0]);
    plot(tvec,hg_pow_mb{iTT}-hg_pow_sdb{iTT},':','Color',[0 1 0]);
    
    title('raw baseline-divided');
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    set(gca,'YLim',[0.5 2]);
    
    % stats interlude
    % becomes a pain because we don't have the full matrix
    
    %
    subplot(223);
    h(1) = plot(tvec,lg_powz_m{iTT},'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(tvec,lg_powz_m{iTT}+lg_powz_sd{iTT},':','Color',[1 0 0]);
    plot(tvec,lg_powz_m{iTT}-lg_powz_sd{iTT},':','Color',[1 0 0]);
    
    h(2) = plot(tvec,hg_powz_m{iTT},'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(tvec,hg_powz_m{iTT}+hg_powz_sd{iTT},':','Color',[0 1 0]);
    plot(tvec,hg_powz_m{iTT}-hg_powz_sd{iTT},':','Color',[0 1 0]);
    
    title('z-scored');
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    set(gca,'YLim',[-1 1]);
    
    %
    subplot(224);
    h(1) = plot(tvec,lg_powz_mb{iTT},'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(tvec,lg_powz_mb{iTT}+lg_powz_sdb{iTT},':','Color',[1 0 0]);
    plot(tvec,lg_powz_mb{iTT}-lg_powz_sdb{iTT},':','Color',[1 0 0]);
    
    h(2) = plot(tvec,hg_powz_mb{iTT},'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(tvec,hg_powz_mb{iTT}+hg_powz_sdb{iTT},':','Color',[0 1 0]);
    plot(tvec,hg_powz_mb{iTT}-hg_powz_sdb{iTT},':','Color',[0 1 0]);
    
    title('z-scored baseline-subtracted');
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    set(gca,'YLim',[-0.5 0.5]);
    
end

%% remove channels to speed up what follows (channel selection doesn't work in ft?)
cfg = [];
cfg.channel = 'vStr';
data_all = ft_selectdata(cfg,data_all);

%% compute spectral version (1)
for iTT = 1:length(trial_types)
    
    cfg              = []; % start with empty cfg
    cfg.output       = 'pow';
    cfg.keeptrials   = 'yes'; % need this for stats later
    cfg.trials       = find(event_all == trial_types(iTT));
    cfg.channel      = {'vStr'};
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = 30:1:100; % frequencies of interest
    cfg.t_ftimwin    = 19*(1./cfg.foi);
    cfg.toi          = twin(1):0.01:twin(2);
    
    TFR{iTT} = ft_freqanalysis(cfg, data_all);
    TFR_cue{iTT} = ft_freqanalysis(cfg, data_all_cue);
    
end

%% plot/stats
vStr_idx = find(strcmp(data_all.label,'vStr'));

F_idx = find(TFR{1}.freq > f_hg(1) & TFR{1}.freq <= f_hg(2)); % high-gamma power
F_idx2 = find(TFR{1}.freq > f_lg(1) & TFR{1}.freq <= f_lg(2)); % low-gamma power
T_idx = find(TFR{1}.time > stat_twin1(1) & TFR{1}.time <= stat_twin1(2)); % pre
T_idx2 = find(TFR{1}.time > stat_twin2(1) & TFR{1}.time <= stat_twin2(2)); % post

clear h nh lg_pwr_pre lg_pwr_post hg_pwr_pre hg_pwr_post lg_pwr_preb lg_pwr_postb hg_pwr_preb hg_pwr_postb;

for iTT = 1:length(trial_types)
    
    hg_pwr = sq(nanmean(sq(TFR{iTT}.powspctrm(:,vStr_idx,F_idx,:)),2)); % average across frequencies
    lg_pwr = sq(nanmean(sq(TFR{iTT}.powspctrm(:,vStr_idx,F_idx2,:)),2));
    
    % stats interlude
    hg_pwr_pre{iTT} = nanmean(hg_pwr(:,T_idx),2); hg_pwr_post{iTT} = nanmean(hg_pwr(:,T_idx2),2);
    lg_pwr_pre{iTT} = nanmean(lg_pwr(:,T_idx),2); lg_pwr_post{iTT} = nanmean(lg_pwr(:,T_idx2),2);
    p_hg = signrank(hg_pwr_pre{iTT},hg_pwr_post{iTT});
    p_lg = signrank(lg_pwr_pre{iTT},lg_pwr_post{iTT});
    fprintf('\n%d raw Hanning lg: pre %.2f, post %.2f, p %3.2e\n',trial_types(iTT),nanmean(lg_pwr_pre{iTT}),nanmean(lg_pwr_post{iTT}),p_lg);
    fprintf('%d raw Hanning hg: pre %.2f, post %.2f, p %3.2e\n',trial_types(iTT),nanmean(hg_pwr_pre{iTT}),nanmean(hg_pwr_post{iTT}),p_hg);
    
    figure(iTT);
    subplot(221);
    h(1) = plot(TFR{iTT}.time,nanmean(lg_pwr),'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(TFR{iTT}.time,nanmean(lg_pwr)+nanstd(lg_pwr),':','Color',[1 0 0]);
    plot(TFR{iTT}.time,nanmean(lg_pwr)-nanstd(lg_pwr),':','Color',[1 0 0]);
    
    h(2) = plot(TFR{iTT}.time,nanmean(hg_pwr),'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(TFR{iTT}.time,nanmean(hg_pwr)+nanstd(hg_pwr),':','Color',[0 1 0]);
    plot(TFR{iTT}.time,nanmean(hg_pwr)-nanstd(hg_pwr),':','Color',[0 1 0]);
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    
    title('raw hanning 15-cycle');
    
    % baseline correct (manual, because relative to cue data)
    baseline = [-1 0];
    t_idxb = nearest(TFR_cue{iTT}.time,baseline);
    temp_baseline = nanmean(sq(TFR_cue{iTT}.powspctrm(:,vStr_idx,:,t_idxb(1):t_idxb(2))),3); % average over time
    temp_baseline = repmat(temp_baseline,[1 1 length(TFR{iTT}.time)]);
    
    all_pwr = sq(TFR{iTT}.powspctrm)./temp_baseline;
    
    lg_pwr = sq(nanmean(all_pwr(:,F_idx2,:),2)); % average across frequencies
    hg_pwr = sq(nanmean(all_pwr(:,F_idx,:),2));
    
    % stats interlude
    hg_pwr_preb{iTT} = nanmean(hg_pwr(:,T_idx),2); hg_pwr_postb{iTT} = nanmean(hg_pwr(:,T_idx2),2);
    lg_pwr_preb{iTT} = nanmean(lg_pwr(:,T_idx),2); lg_pwr_postb{iTT} = nanmean(lg_pwr(:,T_idx2),2);
    p_hg = signrank(hg_pwr_preb{iTT},hg_pwr_postb{iTT});
    p_lg = signrank(lg_pwr_preb{iTT},lg_pwr_postb{iTT});
    fprintf('%d raw-bd Hanning lg: pre %.2f, post %.2f, p %3.2e\n',trial_types(iTT),nanmean(lg_pwr_preb{iTT}),nanmean(lg_pwr_postb{iTT}),p_lg);
    fprintf('%d raw-bd Hanning hg: pre %.2f, post %.2f, p %3.2e\n',trial_types(iTT),nanmean(hg_pwr_preb{iTT}),nanmean(hg_pwr_postb{iTT}),p_hg);
    
    subplot(222);
    h(1) = plot(TFR{iTT}.time,nanmean(lg_pwr),'Color',[1 0 0],'LineWidth',2);
    hold on;
    plot(TFR{iTT}.time,nanmean(lg_pwr)+nanstd(lg_pwr),':','Color',[1 0 0]);
    plot(TFR{iTT}.time,nanmean(lg_pwr)-nanstd(lg_pwr),':','Color',[1 0 0]);
    
    h(2) = plot(TFR{iTT}.time,nanmean(hg_pwr),'Color',[0 1 0],'LineWidth',2);
    hold on;
    plot(TFR{iTT}.time,nanmean(hg_pwr)+nanstd(hg_pwr),':','Color',[0 1 0]);
    plot(TFR{iTT}.time,nanmean(hg_pwr)-nanstd(hg_pwr),':','Color',[0 1 0]);
    legend(h,{'low-gamma','high-gamma'}); legend boxoff;
    set(gca,'YLim',[0.5 2]);
    
    title('baseline-corrected hanning 11-cycle');
    
    figure(3);
    subplot(222);
    nh(iTT) = plot(TFR{iTT}.time,nanmean(lg_pwr),'Color',[1 0 0],'LineWidth',2);
    hold on;
    nh(iTT+2) = plot(TFR{iTT}.time,nanmean(hg_pwr),'Color',[0 1 0],'LineWidth',2);
    
    % put low and high gamma in different subplots
    subplot(223);
    nh2(iTT) = plot(TFR{iTT}.time,nanmean(lg_pwr),'Color',[1 0 0],'LineWidth',2);
    hold on;
    set(gca,'FontSize',14); xlabel('time (s)');
    ylabel('power (baseline-normalized)');
    %legend(nh,{'1p lg','5p lg','1p hg','5p hg'});
    xlim([-1 3]); set(gca,'YLim',[0.5 2]); box off;
        
    subplot(224);
    
    nh2(iTT+2) = plot(TFR{iTT}.time,nanmean(hg_pwr),'Color',[0 1 0],'LineWidth',2);
    hold on;
    set(gca,'FontSize',14); xlabel('time (s)');
    ylabel('power (baseline-normalized)');
    %legend(nh,{'1p lg','5p lg','1p hg','5p hg'});
    xlim([-1 3]); set(gca,'YLim',[0.5 2]); box off;
    
    
end
subplot(222)
box off
set(gca,'FontSize',14); xlabel('time (s)');
ylabel('power (baseline-normalized)');
set(nh(1),'Color',[0.6 0.3 0.3],'LineStyle',':');
set(nh(3),'Color',[0.3 0.6 0.3],'LineStyle',':');
%legend(nh,{'1p lg','5p lg','1p hg','5p hg'});
xlim([-1 3]); set(gca,'YLim',[0.5 2]);

subplot(223);
set(nh2(1),'Color',[0.6 0.3 0.3],'LineStyle',':');
set(nh2(3),'Color',[0.3 0.6 0.3],'LineStyle',':');
subplot(224);
set(nh2(1),'Color',[0.6 0.3 0.3],'LineStyle',':');
set(nh2(3),'Color',[0.3 0.6 0.3],'LineStyle',':');

%% stats
p_lg_pre = ranksum(lg_pwr_preb{1},lg_pwr_preb{2});
p_hg_pre = ranksum(hg_pwr_preb{1},hg_pwr_preb{2});
fprintf('\nraw-bd Hanning lg pre 1p %.2f, 5p %.2f, ranksum p %3.2e\n',nanmean(lg_pwr_preb{1}),nanmean(lg_pwr_preb{2}),p_lg_pre);
fprintf('raw-bd Hanning hg pre 1p %.2f, 5p %.2f, ranksum p %3.2e\n',nanmean(hg_pwr_preb{1}),nanmean(hg_pwr_preb{2}),p_hg_pre);
    
p_lg_post = ranksum(lg_pwr_postb{1},lg_pwr_postb{2});
p_hg_post = ranksum(hg_pwr_postb{1},hg_pwr_postb{2});
fprintf('raw-bd Hanning lg post 1p %.2f, 5p %.2f, ranksum p %3.2e\n',nanmean(lg_pwr_postb{1}),nanmean(lg_pwr_postb{2}),p_lg_post);
fprintf('raw-bd Hanning hg post 1p %.2f, 5p %.2f, ranksum p %3.2e\n',nanmean(hg_pwr_postb{1}),nanmean(hg_pwr_postb{2}),p_hg_post);
    
[~,p_lg_pre] = ttest2(lg_pwr_preb{1},lg_pwr_preb{2});
[~,p_hg_pre] = ttest2(hg_pwr_preb{1},hg_pwr_preb{2});
fprintf('raw-bd Hanning lg pre 1p %.2f, 5p %.2f, ttest p %3.2e\n',nanmean(lg_pwr_preb{1}),nanmean(lg_pwr_preb{2}),p_lg_pre);
fprintf('raw-bd Hanning hg pre 1p %.2f, 5p %.2f, ttest p %3.2e\n',nanmean(hg_pwr_preb{1}),nanmean(hg_pwr_preb{2}),p_hg_pre);
    
[~,p_lg_post] = ttest2(lg_pwr_postb{1},lg_pwr_postb{2});
[~,p_hg_post] = ttest2(hg_pwr_postb{1},hg_pwr_postb{2});
fprintf('raw-bd Hanning lg post 1p %.2f, 5p %.2f, ttest p %3.2e\n',nanmean(lg_pwr_postb{1}),nanmean(lg_pwr_postb{2}),p_lg_post);
fprintf('raw-bd Hanning hg post 1p %.2f, 5p %.2f, ttest p %3.2e\n',nanmean(hg_pwr_postb{1}),nanmean(hg_pwr_postb{2}),p_hg_post);
    
figure(3);
subplot(221);
bv = [nanmean(lg_pwr_postb{1}) nanmean(lg_pwr_postb{2}) nanmean(hg_pwr_postb{1}) nanmean(hg_pwr_postb{2})];
bstd = [nanstd(lg_pwr_postb{1}) nanstd(lg_pwr_postb{2}) nanstd(hg_pwr_postb{1}) nanstd(hg_pwr_postb{2})];
bstd = bstd./sqrt(length(fd)); % SEM
cols = 'rrgg';
for iB = 1:length(bv)
   h = bar(iB,bv(iB));
   set(h,'FaceColor',cols(iB));
   hold on;
end
h = errorbar(1:4,bv,bstd,'k'); set(h,'LineStyle','none');
set(gca,'XTick',1:4,'FontSize',14,'XTickLabel',{'1p','5p','1p','5p'});
ylabel('power (baseline-normalized)'); box off