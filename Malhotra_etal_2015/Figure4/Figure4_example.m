%% example plot for Figure 3
clear all; pack

%% load data
%fd = {'D:\vandermeerlab\R002\R002-2012-04-30_promoted'};
cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
load fd_isidro

for iFD = 12% previously 12

cd(fd{iFD});
run(FindFile('*keys.m'));
fc = ExpKeys.goodGamma(1);
data = ft_read_neuralynx_interp(fc);
data = ft_filterLFP(data,200,'fmode','lowpass');
data.label{1} = 'vStr';

% restrict to TimeOnTrack for normalization
cfg = []; cfg.t = [ExpKeys.TimeOnTrack(1)-5 ExpKeys.TimeOffTrack(end)+5];
data = ft_restrict_data(cfg,data);

%% add hg and lg hilbert for gamma event detection
f_lg = [50 65];
f_hg = [70 85];

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

%% detect gamma events
%rmpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab-old\util\MClust-3.5'));
lgp_z = temp_tsd(data.time{1},data.trial{1}(4,:));
hgp_z = temp_tsd(data.time{1},data.trial{1}(5,:));

cfg = [];
cfg.method = 'raw';
cfg.threshold = 2;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.05; % minimum interval length
 
lg_iv = TSDtoIV(cfg,lgp_z); % detect intervals where z-scored power is > 3 SD above mean
hg_iv = TSDtoIV(cfg,hgp_z); % detect intervals where z-scored power is > 3 SD above mean

%% create copy of signal with non-gamma NaNd out
lg_ivx = InvertIV(lg_iv,lgp_z);
data.trial{1}(6,:) = data.trial{1}(1,:);
data.label{6} = 'lg_evt';
for iEvt = 1:length(lg_ivx.tstart)
    t_idx = find(data.time{1} > lg_ivx.tstart(iEvt) & data.time{1} <= lg_ivx.tend(iEvt));
    data.trial{1}(6,t_idx) = NaN;
end

hg_ivx = InvertIV(hg_iv,hgp_z);
data.trial{1}(7,:) = data.trial{1}(1,:);
data.label{7} = 'hg_evt';
for iEvt = 1:length(hg_ivx.tstart)
    t_idx = find(data.time{1} > hg_ivx.tstart(iEvt) & data.time{1} <= hg_ivx.tend(iEvt));
    data.trial{1}(7,t_idx) = NaN;
end

%% trialification
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2; cfg.trialdef.post = 4;

cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'value'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)

cfg.ExpKeys = ExpKeys;

% cue presentation
[trl, event] = ft_trialfun_lineartracktone2(cfg); cfg.trl = trl;
data_np_c1 = ft_redefinetrial(cfg,data);

%% plot
figure;

cfg = []; cfg.channel = 'vStr';
cfg.xlim = [-1 2];
cfg.yscale = 200; % rescaling of each trial
cfg.ystep = 2; % spacing between trials
cfg.lw = 1; % line width
cfg.plotaverage = 0;
cfg.plotcolor = [0.7 0.7 0.7]; 
%cfg.trials = 9:29;

ft_mvdmlab_trialplot(cfg,data_np_c1); % should return scaling factors for overplotting
hold on;

cfg.channel = 'lg_evt';
cfg.plotcolor = [1 0 0];
ft_mvdmlab_trialplot(cfg,data_np_c1); % should return scaling factors for overplotting

cfg.channel = 'hg_evt';
cfg.plotcolor = [0 1 0];
ft_mvdmlab_trialplot(cfg,data_np_c1); % should return scaling factors for overplotting

title(fd{iFD});
end

%% spectogram
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'vStr';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:0.5:100;
cfg.keeptrials   = 'yes'; % need this for stats later
cfg.t_ftimwin    = 15./cfg.foi;  % 20 cycles per time window
 
cfg.toi          = -2:0.01:4; % pre-nosepoke baseline
TFR = ft_freqanalysis(cfg, data_np_c1);

%%
figure
cfg = []; cfg.channel = 'vStr';
cfg.xlim = [-0.5 2];
ft_singleplotTFR(cfg, TFR);

%%
figure;
cfg.baseline     = [-2 0];
cfg.xlim = [-0.5 2];
cfg.baselinetype = 'relative';
cfg.channel = 'vStr';
ft_singleplotTFR(cfg, TFR);