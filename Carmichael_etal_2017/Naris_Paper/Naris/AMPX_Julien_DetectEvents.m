function [evts, detect_csc, detect_chan] = AMPX_Julien_DetectEvents(data, ExpKeys, varargin)

%% MASTER_CollectGammaEvents.m
% detects and organizes gamma events from raw LFP data
%
% Julien Catanese & Matthijs van der Meer

Remove_ft()
%% Convert data to TSD and extract the channel for detection from the ExpKeys

if exist('detect_chan') ~=1
    detect_chan = Naris_BestChan_remap(ExpKeys, 'location', 'vl');
end

data_tsd = AMPX_to_tsd(data);
csc = data_tsd;
csc.data = csc.data(detect_chan,:);
csc.detect_chan = detect_chan;
%% set params
% gamma event detection
PARAM_f_label = {'low','high'};%, 'low_98_tr', 'high_98_tr'};
PARAM_f_bandpass = {[45 65],[70 90],[45 65], [70 90]}; % frequency bands for event detection
PARAM_detect_thr = [0.95 .95 0.98 .98]; % [2.5 2.5 1 1]; %threshold for event detection: 95th percentile of (amplitude) envelope
PARAM_detect_method = 'percentile'; % 'raw', 'zscore', 'percentile'
PARAM_detect_nCycles = 4; % require minimum number of gamma cycles
PARAM_var_thr = 1.5; % def: variance/mean of cycle peaks and throughs must be smaller than this
PARAM_detect_epoch = 'all'; % 'all', 'post', 'task'; % set threshold based on what data (for events)
if strcmp(ExpKeys.Subject, 'R054')
    PARAM_ampl_min = [40 30 ]; % all peaks of event must be larger than this (in V)
    PARAM_mean_filt = [100 40 ]; % all peaks of event must be larger than this (in V)
elseif strcmp(ExpKeys.Subject, 'R061')
    PARAM_ampl_min = [40 30 ]; % all peaks of event must be larger than this (in V)
    PARAM_mean_filt = [100 40 ]; % all peaks of event must be larger than this (in V)
elseif strcmp(ExpKeys.Subject, 'R045')
    PARAM_detect_thr = [0.95 .98 0.98 .98]; % [2.5 2.5 1 1]; %threshold for event detection: 95th percentile of (amplitude) envelope
    PARAM_ampl_min = [100 35]; % all peaks of event must be larger than this (in V)
    PARAM_mean_filt = [100 40]; % all peaks of event must be larger than this (in V)
else
    PARAM_ampl_min = [100 35]; % all peaks of event must be larger than this (in V)
    PARAM_mean_filt = [100 40]; % all peaks of event must be larger than this (in V)
end
PARAM_var_raw_max = [1, 1, 1, 1];
% artifact, chewing, and spindle detection
PARAM_artif_thr =  std(csc.data)*4;   %0.75 * 10^-3; % raw amplitude must be smaller than this (in V) to pass artifact detection
PARAM_chew_thr = 3; % z-score of chew band envelope must be smaller than this, default 0.25
PARAM_spindle_thr = 4; % z-score of spindle band envelope must be smaller than this

% some flags to enable visualization
debug = 0; debug2 = 0;

extract_varargin


%% detect major transient artifacts: this is on UNFILTERED data because of sustained artifacts in R026 (see project log)
csc_artif = csc;
csc_artif.data = abs(csc_artif.data); % detect artifacts both ways

cfg_artif_det = [];
cfg_artif_det.method = 'raw';
cfg_artif_det.threshold = PARAM_artif_thr;
cfg_artif_det.minlen = 0;
evt_artif = TSDtoIV(cfg_artif_det,csc_artif);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
evt_artif = ResizeIV(cfg_temp,evt_artif);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_artif,csc);
    pause; close all;
end

%% NaN out artifacts to improve reliability of subsequent z-scoring
artif_idx = TSD_getidx2(csc,evt_artif); % if error, try TSD_getidx (slower)
csc.data(artif_idx) = NaN;

%% find chewing artifacts
cfg_chew = [];
cfg_chew.epoch = 'all'; % chewing occurs during task mostly
cfg_chew.minlen = 0.02;
cfg_chew.filter_cfg.f = [200 500]; % default [200 300]
cfg_chew.threshold = PARAM_chew_thr; % 0.25 for session 2, 0.5 for session 1?
cfg_chew.smooth = 0.05; % convolve with Gaussian of this SD
evt_chew = Julien_DetectEvents(cfg_chew,csc,ExpKeys);

cfg_temp = []; cfg_temp.d = [-0.1 0.1];
evt_chew = ResizeIV(cfg_temp,evt_chew);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_chew,csc);
    pause; close all;
end

%% NaN out artifacts to improve reliability of subsequent z-scoring
chew_idx = TSD_getidx2(csc,evt_chew); % if error, try TSD_getidx (slower)
csc.data(chew_idx) = NaN;
%% find High voltage spindles to exclude from the data. 
cfg_spindl = [];
cfg_spindl.epoch = 'all'; % chewing occurs during task mostly
cfg_spindl.minlen = 0.005;
cfg_spindl.filter_cfg.f = [7 11]; % default [25 35]
cfg_spindl.threshold = PARAM_spindle_thr; % 0.25 for session 2, 0.5 for session 1?
cfg_spindl.smooth = 0.05; % convolve with Gaussian of this SD
evt_spindl = Julien_DetectEvents(cfg_spindl,csc,ExpKeys);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
evt_spindl = ResizeIV(cfg_temp,evt_spindl);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_spindl,csc);
    pause; close all;
end
%% now, loop over frequency bands to process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iFreq = 1:length(PARAM_f_label)
    
    %% set up filter
    cfg_filter = [];
    cfg_filter.f = PARAM_f_bandpass{iFreq};
    cfg_filter.type = 'cheby1';
    cfg_filter.order = 5;
    
    % basic gamma detection
    cfg_evt = [];
    cfg_evt.epoch = PARAM_detect_epoch;
    cfg_evt.epochLength = 10*60; % was 5 * 60
    cfg_evt.filter_cfg = cfg_filter;
    cfg_evt.minlen = PARAM_detect_nCycles./mean(PARAM_f_bandpass{iFreq}); % or, 0.05
    cfg.merge_thr = 0.025; % merge events closer than this
    %cfg_evt.minlen = 0.05;
    cfg_evt.epochLength = 0;
    cfg_evt.smooth = 0.05; % convolve with Gaussian of this SD
    cfg_evt.threshold = PARAM_detect_thr(iFreq);
    cfg_evt.method = PARAM_detect_method;
    
    [evt,evt_thr] = Julien_DetectEvents(cfg_evt,csc,ExpKeys);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events detected initially.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    if debug
        cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center'; cfg_plot.width = 0.2;
        PlotTSDfromIV(cfg_plot,evt,csc);
        pause; close all;
    end
    
    % remove artifacts
    evt = DifferenceIV([],evt,evt_artif);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after artifact removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    % remove chewing
    evt = DifferenceIV([],evt,evt_chew);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after chewing removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    % remove spindles
    evt = DifferenceIV([],evt,evt_spindl);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after spindle removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    
    %% exclude events with insufficient gamma cycles - count how many exist above same threshold as used for detection
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filter;
    evt = CountCycles(cfg_cc,csc,evt);
    
    cfg_cc = [];
    cfg_cc.operation = '>=';
    cfg_cc.threshold = PARAM_detect_nCycles-1;
    evt = SelectIV(cfg_cc,evt,'nCycles');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after cycle count removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    % exclude events with excessive variance in amplitude (only those in within evt detect boundaries)
    %         iv_temp = IVcenters(evt);
    %         iv_temp = iv(iv_temp-cfg_evt.minlen/2,iv_temp+cfg_evt.minlen/2);
    %
    %         cfg_cc = [];
    %         cfg_cc.threshold_type = 'raw';
    %         cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    %         cfg_cc.filter_cfg = cfg_filter;
    %         iv_temp = CountCycles(cfg_cc,csc,iv_temp);
    %
    %         cfg_cc = [];
    %         cfg_cc.operation = '<';
    %         cfg_cc.threshold = PARAM_var_thr;
    %         [~,keep_idx] = SelectIV(cfg_cc,iv_temp,'var_raw');
    %
    %         evt = SelectIV([],evt,keep_idx); iv_temp = SelectIV([],iv_temp,keep_idx);
    %
    %         cfg_cc = [];
    %         cfg_cc.operation = '<';
    %         cfg_cc.threshold = PARAM_var_thr;
    %         [~,keep_idx] = SelectIV(cfg_cc,iv_temp,'var');
    %
    %         evt = SelectIV([],evt,keep_idx);
    
    
    % exclude events with excessive variance in amplitude (all peaks and troughs)
    cfg_cc = [];
    cfg_cc.operation = '<';
    cfg_cc.threshold = PARAM_var_thr;
    evt = SelectIV(cfg_cc,evt,'var_raw');
    
    cfg_cc = [];
    cfg_cc.operation = '<';
    cfg_cc.threshold = PARAM_var_thr;
    evt = SelectIV(cfg_cc,evt,'var');
    
    if strcmp(ExpKeys.Subject, 'R049') % add artifact where some false positives are picked up and have very low variance for some reason
        cfg_cc = [];
        cfg_cc.operation = '>';
        cfg_cc.threshold = 0.1;
        evt = SelectIV(cfg_cc,evt,'var');
    end
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after amplitude variance removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    %% exclude events with insufficient mean (or min) -- could try doing this on minlen part only
    cfg_cc = [];
    cfg_cc.operation = '<';
    cfg_cc.threshold = PARAM_var_raw_max(iFreq);
    evt = SelectIV(cfg_cc,evt,'var_raw');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after mean peak removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    % minlen only version
    %     iv_temp = IVcenters(evt);
    %     iv_temp = iv(iv_temp-cfg_evt.minlen/2,iv_temp+cfg_evt.minlen/2);
    
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filter;
    iv_temp = CountCycles(cfg_cc,csc,evt);
    
    %     %%
    %     cfg_cc = [];
    %     cfg_cc.operation = '>';
    %     cfg_cc.threshold = PARAM_ampl_min(iFreq);
    %     [evt,keep_idx] = SelectIV(cfg_cc,iv_temp,'min_filt');
    %
    if PARAM_f_bandpass{iFreq}(1) >= 70
        fprintf('\n This is a high gamma detection and will use the mean_freq threshold to remove false positives\n')
        cfg_cc = [];
        cfg_cc.operation = '>';
        cfg_cc.threshold = PARAM_mean_filt(iFreq);
        evt = SelectIV(cfg_cc,iv_temp,'mean_filt');
    end
    %
    %     evt = SelectIV([],evt,keep_idx);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after mean peak removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    
    %% visualize
    if debug2
        % raw LFP only
        close all
        cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center'; cfg_plot.width = 0.25; cfg_plot.title = 'var';
        PlotTSDfromIV(cfg_plot,evt,csc);
        cfg_plot.fname = cd;
        cfg_plot.fname = cfg_plot.fname(end-14:end)
        if exist(['C:\temp\Gamma_detection\' num2str(cfg_plot.fname)]) ~=1
            mkdir(['C:\temp\Gamma_detection\' num2str(cfg_plot.fname)]);
        end
        h = get(0, 'children');
        for ii = 1:length(h)
            saveas(h(ii), ['C:\temp\Gamma_detection\' num2str(cfg_plot.fname) '\' PARAM_f_label{iFreq} '_' num2str(ii), 'fig'])
            saveas(h(ii), ['C:\temp\Gamma_detection\' num2str(cfg_plot.fname) '\' PARAM_f_label{iFreq} '_' num2str(ii), '.png'])
        end
        close all;
        
        
        %%
%         h3 = figure(2000)
%         set(gcf,'PaperPositionMode','auto')
%         set(gcf, 'position', [48 53 1824 936])
%         csc_example = restrict(data_tsd, evt.tstart(1), evt.tend(1))
%         offsets = -2000:400:2000;
%         c_ord = linspecer(8);
%         for isub = 1:8
%             off = 1;
%             subtightplot(1,8,isub);
%             chans = fliplr((1:8)+((isub-1)*8));
%             for iChan = chans
%                 if isub == 4 && iChan == 25;
%                     title([ExpKeys.HSIds{1} ' ' ExpKeys.HSIds{2} ' ' ExpKeys.HSConnection{1} ' ' ExpKeys.HSConnection{2}])
%                 end
%                 hold on
%                 if ismember(iChan, ExpKeys.BadChannels)
%                     plot(csc_example.tvec, csc_example.data(iChan,:)+offsets(off), 'color', [1 1 1])
%                 else
%                     plot(csc_example.tvec, csc_example.data(iChan,:)+offsets(off), 'color', c_ord(off,:))
%                 end
%                 text(csc_example.tvec(1), csc_example.data(iChan,1)+offsets(off),['\bf ' num2str(iChan)])
%                 off = off+1;
%             end
%             axis off
%         end
%         saveas(h3, ['C:\temp\Gamma_detection\' num2str(cfg_plot.fname) '\' PARAM_f_label{iFreq} '_trace', 'fig'])
%         print(h3, ['C:\temp\Gamma_detection\' num2str(cfg_plot.fname) '\' PARAM_f_label{iFreq} '_trace' ], '-dpng', '-r300')
        %%
        close all;
        
        % TFR version
        %             cfg_convert = []; cfg_convert.mode = 'resample';
        %             csc_ft = TSDtoFT(cfg_convert,csc);
        %
        %             evt_temp = evt;
        %             evt_temp.tstart = evt_temp.tstart-csc.tvec(1);
        %             evt_temp.tend = evt_temp.tend-csc.tvec(1);
        %
        %             cfg_temp = []; cfg_temp.foi = 1:5:300; cfg_temp.clim = [0 10^-9]; cfg_temp.twin = [-0.1 0.1];
        %             PlotTSDfromIV_TFR(cfg_temp,evt_temp,csc_ft);
        %             pause; close all;
    end
    
    %% store data into collector variable
    evts.(PARAM_f_label{iFreq}) = evt;
    evts.(PARAM_f_label{iFreq}).firstTimestamp = csc.tvec(1);
    evt_thr.(PARAM_f_label{iFreq}) = evt_thr; 
    %         sess_id_field = regexprep(SessionLIST{iFD},'-','_');
    %         ALL_evt.(this_RatID).(sess_id_field).(PARAM_f_label{iFreq}) = evt;
    %         ALL_evt.(this_RatID).(sess_id_field).fd = fd{iFD};
    %         ALL_evt.(this_RatID).(sess_id_field).firstTimestamp = csc.tvec(1); % need this for trialification later
    
end % of loop over frequencies

%% remove overlapping gamma events.
temp_lg = []; temp_hg = []; %temp_lg_98 = []; temp_hg_98 = [];

temp_lg = DifferenceIV([],evts.low,evts.high);
temp_hg = DifferenceIV([],evts.high,evts.low);

% temp_lg_98 = DifferenceIV([],evts.low_98_tr,evts.high_98_tr);
% temp_hg_98 = DifferenceIV([],evts.high_98_tr,evts.low_98_tr);


evts.low = temp_lg;
evts.high = temp_hg;
% evts.high_98 = temp_hg_98;
% evts.low_98 = temp_lg_98;

 %% find high voltage spindles for PCA
%  cfg_filter = [];
%     cfg_filter.f = [7 12];
%     cfg_filter.type = 'butter';
%     cfg_filter.order = 4;
% cfg_spindl = [];
% cfg_spindl.epoch = 'all'; % chewing occurs during task mostly
% cfg_spindl.minlen = 0.5;
% cfg_spindl.filter_cfg = cfg_filter; % default [25 35]
% cfg_spindl.threshold = 2; % 0.25 for session 2, 0.5 for session 1?
% % cfg_spindl.smooth = 0.05; % convolve with Gaussian of this SD
% [evt_spindl, evt_spin_thr] = Julien_DetectEvents(cfg_spindl,csc,ExpKeys);
% 
%     cfg_cc = [];
%     cfg_cc.threshold_type = 'raw';
%     cfg_cc.threshold = evt_spin_thr; % use same threshold as for orignal event detection
%     cfg_cc.filter_cfg = cfg_filter;
%     evts.spindles = CountCycles(cfg_cc,csc,evt_spindl);
%     
%     cfg_cc = [];
%     cfg_cc.operation = '>=';
%     cfg_cc.threshold = 8;
%     evts.spindles = SelectIV(cfg_cc,evts.spindles,'nCycles');
% 
% evts.spindles.firstTimestamp = csc.tvec(1);  
% cfg_plot.display = 'iv'; cfg_plot.title = 'nCycles'; % tsd, iv
%     PlotTSDfromIV(cfg_plot,evts.spindles,csc);
%     close all;
%%


% for iFreq = 1:length(PARAM_f_label)
%     fprintf(['\n ' PARAM_f_label{iFreq} ' found: ' num2str(length(evts.(PARAM_f_label{iFreq}).tstart)) ' events'])
% end

% fprintf(['\nSpindles found:'  num2str(length(evts.spindles.tstart)) ' events'])
% fprintf('\n')

% output the same csc used to detect the events (used for pseudo random
% non-gamma epochs selection
detect_csc = csc;
