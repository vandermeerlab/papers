%% PLOT_AmplCorr.m
% plots amplitude cross-correlations on gamma events
% generated by MASTER_AmplCorr
%
% Julien Catanese & Matthijs van der Meer
%
% TODO make function out of actual plotter

%% set paths
restoredefaultpath;
cd('D:\My_Documents\GitHub\fieldtrip');
ft_defaults;

rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal\');

addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\tasks\Julien_linear_track')); % Detect events, CountCycles live here

%% load gamma events (use MASTER_CollectGammaEvents.m to obtain) -- puts ALL_evt variable in workspace
cd('D:\My_Documents\Dropbox\projects\Julien_multiLFP\2016-01-07');
%load(FindFile('gamma*.mat'));

%% define what to run
rats = {'R026','R032','R033','R039'};

%%
PARAM_Fs = 2000; % sampling rate; used to convert xcorr lags (in samples) into seconds

PARAM_shuffleTest = 1; % use shuffle test to determine event inclusion
PARAM_shuffle_alpha = 0.05; % significance threshold for vStr-mPFC independence test
PARAM_cdfx = -5:0.01:5;
PARAM_cdfy = normcdf(PARAM_cdfx,0,1);

%%
available_rats = fieldnames(ALL_evt);
nRats = 0; nSessions = 0;

all_lg = [];
all_hg = [];

for iRat = 1:length(rats)
    
    this_rat = rats{iRat};
    
    if ~strmatch(this_rat,available_rats)
       warning('Rat %s not available -- skipping...',rats{iRat});
       continue;
    end
    
    rat_lg = []; rat_hg = [];
    
    available_sessions = fieldnames(ALL_evt.(this_rat));
    
    for iSession = 1:length(available_sessions)
    
        nSessions = nSessions + 1;
        
        this_session = available_sessions{iSession};
        this_session_data = ALL_evt.(this_rat).(this_session);
        
        this_lg = this_session_data.xcorr.lg;
        this_hg = this_session_data.xcorr.hg;
        shift_values = this_session_data.xcorr.shift_values .* (1./PARAM_Fs);
        
        if PARAM_shuffleTest % only keep those events that pass non-independence test vStr-mPFC
            % otherwise xcorr doesn't have clear interpretation
            nEvents_lg = size(this_lg,1); nEvents_hg = size(this_hg,1);
            
            % find threshold in number of SDs
            temp_idx_lg = nearest_idx3(1-(PARAM_shuffle_alpha./nEvents_lg),PARAM_cdfy);
            SD_lg = PARAM_cdfx(temp_idx_lg);
            
            temp_idx_hg = nearest_idx3(1-(PARAM_shuffle_alpha./nEvents_hg),PARAM_cdfy);
            SD_hg = PARAM_cdfx(temp_idx_hg);
            
            % convert to correlation threshold
            corr_thr_lg = this_session_data.xcorr.shuf_lg_maxcorr_sd.*SD_lg;
            corr_thr_hg = this_session_data.xcorr.shuf_hg_maxcorr_sd.*SD_hg;
            
            keep_idx_lg = max(this_session_data.xcorr.lg,[],2) > corr_thr_lg';
            fprintf('\nkept %d/%d lg events.\n',sum(keep_idx_lg),nEvents_lg);
            this_lg = this_lg(keep_idx_lg,:);
            
            keep_idx_hg = max(this_session_data.xcorr.hg,[],2) > corr_thr_hg';
            fprintf('\nkept %d/%d hg events.\n',sum(keep_idx_hg),nEvents_hg);
            this_hg = this_hg(keep_idx_hg,:);
        end
        %%
        figure(nSessions)
        cfg_plot = [];
        cfg_plot.title_string = this_session;
        FUNC_plotAmplCorr(cfg_plot,shift_values,this_lg,this_hg)
        
        %% add to rat
        rat_lg = cat(1,rat_lg,this_lg);
        rat_hg = cat(1,rat_hg,this_hg);
        
    end % over sessions
    
    %% plot rat
    figure(10+iRat)
    
    cfg_plot = [];
    cfg_plot.title_string = this_rat;
    FUNC_plotAmplCorr(cfg_plot,shift_values,rat_lg,rat_hg)

    %% add to all
    all_lg = cat(1,all_lg,rat_lg);
    all_hg = cat(1,all_hg,rat_hg);
    
end % over rats


%% plot all
figure(100)
cfg_plot = [];
cfg_plot.title_string = 'ALL';
out = FUNC_plotAmplCorr(cfg_plot,shift_values,all_lg,all_hg);

%% av
c1 = [1 0 0]; c2 = [0 0.2 0.5];
figure(101);
subplot(421)
lgm = nanmean(all_lg);
h1 = plot(shift_values*1000,lgm,'LineWidth',2,'Color',c1);
hold on;
hgm = nanmean(all_hg);
h2 = plot(shift_values*1000,hgm,'LineWidth',2,'Color',c2);
set(gca,'FontSize',14,'LineWidth',1);
axis([-10 10 0.5 0.75]); grid on; xlabel('time (ms)');
lh = legend([h1 h2],{'lg','hg'}); legend boxoff

%%
c1 = [1 0 0]; c2 = [0 0.2 0.5];
figure(101);
subplot(423)
lgm = histc(out.lgp,shift_values(1:2:end)*1000)./length(out.lgp);
hgm = histc(out.hgp,shift_values(1:2:end)*1000)./length(out.hgp);
h = bar(shift_values(1:2:end)*1000,[lgm; hgm]');

set(gca,'FontSize',14,'LineWidth',1);
axis([-10 10 0 0.25]); grid on; xlabel('time (ms)');
lh = legend(h,{'lg','hg'}); legend boxoff
set(h(1),'FaceColor',c1,'EdgeColor','none'); set(h(2),'FaceColor',c2,'EdgeColor','none');

%% compare directly
fprintf('lg mean %1.3e, median %1.3e\n',nanmean(out.lgp),median(out.lgp));
fprintf('hg mean %1.3e, median %1.3e\n',nanmean(out.hgp),median(out.hgp));
[p,h] = ranksum(out.lgp,out.hgp)

