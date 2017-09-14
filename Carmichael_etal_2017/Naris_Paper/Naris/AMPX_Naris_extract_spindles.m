% Sandbox to extract spindle power gradients
session_list = fieldnames(all_data)
%%
for iSess = 1:length(session_list)
    Remove_ft(); % make sure fieldtrip is not in the path.  this can cause issues with filters
    %% prep data in same manner as AMPX_Naris_pipeline
    fname = strrep(session_list{iSess},'_','-');
    cd(['D:\DATA\' fname(1:4) '\' fname(1:15) ])
    [data, data_ft] = AMPX_Naris_preprocess([],strrep(session_list{iSess}, '_', '-'),'pre');
    
    LoadExpKeys
    data_remap_AMPX = AMPX_remap_sites(data, ExpKeys);
    
    ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys, data_remap_AMPX);
    detect_chan = Naris_BestChan_remap(ExpKeys, 'location', 'dm');
    csc = AMPX_to_tsd(data_remap_AMPX);
    csc.data = csc.data(detect_chan,:);
    csc.detect_chan = detect_chan;
    % removed artifacts
    csc_artif = csc;
    csc_artif.data = abs(csc_artif.data); % detect artifacts both ways
    
    cfg_artif_det = [];
    cfg_artif_det.method = 'raw';
    cfg_artif_det.threshold = std(csc.data)*4;
    cfg_artif_det.minlen = 0;
    evt_artif = TSDtoIV(cfg_artif_det,csc_artif);
    
    cfg_temp = []; cfg_temp.d = [-0.5 0.5];
    evt_artif = ResizeIV(cfg_temp,evt_artif);
    
    %%
    cfg_chew = [];
    cfg_chew.epoch = 'all'; % chewing occurs during task mostly
    cfg_chew.minlen = 0.02;
    cfg_chew.filter_cfg.f = [200 500]; % default [200 300]
    cfg_chew.threshold = 2; % 0.25 for session 2, 0.5 for session 1?
    cfg_chew.smooth = 0.05; % convolve with Gaussian of this SD
    evt_chew = Julien_DetectEvents(cfg_chew,csc,ExpKeys);
    
    cfg_temp = []; cfg_temp.d = [-0.1 0.1];
    evt_chew = ResizeIV(cfg_temp,evt_chew);
    %% NaN out artifacts to improve reliability of subsequent z-scoring
    
    artif_idx = TSD_getidx2(csc,evt_artif); % if error, try TSD_getidx (slower)
    csc.data(artif_idx) = NaN;
    chew_idx = TSD_getidx2(csc,evt_chew); % if error, try TSD_getidx (slower)
    csc.data(chew_idx) = NaN;
    %%
    % spindles
    cfg_spindl = [];
    cfg_spindl.minlen = 0.1;
    cfg_spindl.filter_cfg.f = [8 15]; % based on Berke 2004
    cfg_spindl.filter_cfg.order = 4;  % methods used for theta detection
    cfg_spindl.filter_cfg.display_filter = 0;% methods used for theta detection
    cfg_spindl.filter_cfg.type = 'fdesign';% methods used for theta detection
    cfg_spindl.threshold = 2;
    cfg_spindl.smooth = 0.05; % convolve with Gaussian of this SD
    [evt_spindl,evt_thr] = Julien_DetectEvents(cfg_spindl,csc,ExpKeys);
    evt_spindl.firstTimestamp = csc.tvec(1);
    % check cycle count and get other metrics
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_spindl.filter_cfg;
    evt_spindl = CountCycles(cfg_cc,csc,evt_spindl);
    % select events with some criteria
    %     cfg_sel = [];
    %     cfg_sel.operation = '<';
    %     cfg_sel.threshold = 1.5;
    %     evt_spindl = SelectIV(cfg_sel,evt_spindl,'var_raw');
    % exclude events with less than 4 cycles
    cfg_cc = [];
    cfg_cc.operation = '>=';
    cfg_cc.threshold = 4;
    evt_spindl = SelectIV(cfg_cc,evt_spindl,'nCycles');
    %use the mean filt .usr parameter to remove false positives.
    cfg_mf = [];
    cfg_mf.operation = '>=';
    cfg_mf.threshold = 75;
    evt_spindl = SelectIV(cfg_mf,evt_spindl,'mean_filt');
    %use the mean filt .usr parameter to remove false positives.
    cfg_mf = [];
    cfg_mf.operation = '<=';
    cfg_mf.threshold = 400;
    evt_spindl = SelectIV(cfg_mf,evt_spindl,'mean_filt');
    % check them
    cfg_plot = [];
    cfg_plot.display = 'iv';
    cfg_plot.mode = 'edges';
    cfg_plot.width = 0.5;
    cfg_plot.title = 'mean_filt';
    PlotTSDfromIV(cfg_plot, evt_spindl, csc)
    
    %% convert spindles that are longer than a gamma event into multiple evets
    Split_IV = evt_spindl;
    Split_IV.tstart = [];
    Split_IV.tend = [];
    Split_IV.usr = [];
    for iEvt = 1:length(evt_spindl.tstart)
        if  (evt_spindl.tend(iEvt) - evt_spindl.tstart(iEvt)) > .4 % is it longer than a typical gamma event
            segments = floor((evt_spindl.tend(iEvt) - evt_spindl.tstart(iEvt))/0.4);
            temp_IV = [];
            for iSeg = 1:segments
                if iSeg == 1
                    temp_IV.tstart(iSeg) = evt_spindl.tstart(iEvt);
                    temp_IV.tend(iSeg) = evt_spindl.tstart(iEvt) +0.4;
                else
                    temp_IV.tstart(iSeg) = temp_IV.tend(iSeg-1);
                    temp_IV.tend(iSeg) = temp_IV.tstart(iSeg) + 0.4;
                end
            end
            Split_IV.tstart = [Split_IV.tstart, temp_IV.tstart];
            Split_IV.tend = [Split_IV.tend, temp_IV.tend];
        else
            Split_IV.tstart = [Split_IV.tstart, evt_spindl.tstart(iEvt)];
            Split_IV.tend = [Split_IV.tend, evt_spindl.tend(iEvt)];
        end
    end
    fprintf(['\n' num2str(evt_spindl.tstart') '\n' num2str(Split_IV.tstart) '\n'])
    Split_IV.tstart = Split_IV.tstart';
    Split_IV.tend = Split_IV.tend';
    
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_spindl.filter_cfg;
    Split_IV = CountCycles(cfg_cc,csc,Split_IV);
    
    cfg_cc = [];
    cfg_cc.operation = '>=';
    cfg_cc.threshold = 2;
    Split_IV = SelectIV(cfg_cc,Split_IV,'nCycles');
    % check again
    cfg_plot = [];
    cfg_plot.display = 'iv';
    cfg_plot.mode = 'edges';
    cfg_plot.width = 0.5;
%     cfg_plot.title = 'nCycles';
    PlotTSDfromIV(cfg_plot, Split_IV, csc)
    
    close
    
    evt_spindl = Split_IV;
    %%
    % Convert Spindle IVs to the FT format for power extract
    if isempty(evt_spindl.tstart) ~=1
        clear data_remap_AMPX data
        data_remap_FT = AMPX_remap_sites(data_ft, ExpKeys);
        clear data_ft
        
        cfg = [];
        cfg.twin = 0.25;
        cfg.check = 1; % make a plot of 9 random trials to see the quality
        Add_ft()
        [data_spindles_ft] = AMPX_to_FT_trial_split(cfg, evt_spindl, data_remap_FT);
        
        % extract power from the spindle evets
        cfg = [];
        cfg.freq = cfg_spindl.filter_cfg.f;
        all_data.(session_list{iSess}).spindles.power = AMPX_get_power(cfg, data_spindles_ft, ExpKeys_remap);
        all_data.(session_list{iSess}).spindles.evts = evt_spindl;
    end
    figure(100+iSess)
    if length(evt_spindl.tstart) > 16
        plots = 16;
    else
        plots = length(evt_spindl.tstart);
    end
    pow = 1:2:32;
    evt = 2:2:32;
    for iEvt = 1:plots
        subplot(8,8,pow(iEvt))
        %         all_data.(session_list{iSess}).spindles.power.power_distrib{iEvt}(8,1) =NaN;
        nan_imagesc_ec(all_data.(session_list{iSess}).spindles.power.power_distrib{iEvt})
        subplot(8,8,evt(iEvt))
        plot(data_spindles_ft.time{iEvt}, data_spindles_ft.trial{iEvt}(detect_chan,:), 'linewidth', 1, 'color', 'k')
        xlim([data_spindles_ft.time{iEvt}(1) data_spindles_ft.time{iEvt}(end)])
    end
    
    clear data_remap_ft
    
end
all_data_spin = all_data;

%% append to the all_data_pre
load('C:\temp\Naris_all_data_pre_Paper.mat')
for iSess = 1:length(session_list)
    if isfield(all_data_spin, session_list{iSess})
        all_data.(session_list{iSess}).spindles = all_data_spin.(session_list{iSess}).spindles;
    else
        all_data.(session_list{iSess}).spindles = [] ;
    end
end


%% collect all the events across all the detection channels by combining all the IVs and then comparing them to the VL
fname = 'R049-2014-02-07';
fname = strrep(fname, '_', '-');
cd(['D:\DATA\' fname(1:4) '\' fname(1:15) ])
[data, data_ft] = AMPX_Naris_preprocess([],fname,'pre');

LoadExpKeys
% Organize the data to fit the probe lay out in the 8x8 format.
% Data about the probemapping

data_remap_AMPX = AMPX_remap_sites(data, ExpKeys)

data_remap_FT = AMPX_remap_sites(data_ft, ExpKeys)
clear data_ft

ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys, data_remap_AMPX);



detect_chan = fieldnames(detect_all);
all_events = []; 
bands = {'low', 'high', 'ctrl_low', 'ctrl_high', 'spindles'};

for iChan = 1:length(detect_chan)
    if iChan ==1
        all_events = detect_all.(detect_chan{iChan});
        all_events.spindles = all_events.spindles.evts;
    end
    for ibands = 1:length(bands)
        if isfield(detect_all.(detect_chan{iChan}), bands{ibands})
            if ibands <5
                all_events.(bands{ibands}) = UnionIV([], all_events.(bands{ibands}), detect_all.(detect_chan{iChan}).(bands{ibands}));
            else
                all_events.(bands{ibands}) = UnionIV([], all_events.(bands{ibands}), detect_all.(detect_chan{iChan}).(bands{ibands}).evts);
            end
        end
    end
end

VL_events.low = all_data.R049_2014_02_07.lg.evts;
VL_events.high = all_data.R049_2014_02_07.hg.evts;
VL_events.ctrl_low = all_data.R049_2014_02_07.lg_ran.evts;
VL_events.ctrl_high = all_data.R049_2014_02_07.hg_ran.evts;
VL_events.spindles = all_data.R049_2014_02_07.spindles.evts;

    for ibands = 1:length(bands)
        no_VL_events.(bands{ibands}) = DifferenceIV([],all_events.(bands{ibands}), VL_events.(bands{ibands}));
    end

    detect_all.VL_events = VL_events;
    detect_all.no_VL_events = no_VL_events;

    temp = DifferenceIV([], no_VL_events.(bands{ibands}), no_VL_events.(bands{ibands}))

%% define FT trials based on evts
addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip')
ft_defaults
% low gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_lg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.low, data_remap_FT);


% high Gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_hg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.high, data_remap_FT);


% low control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_lg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.ctrl_low, data_remap_FT);


% high Control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_hg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.ctrl_high, data_remap_FT);

% Spindles
if isfield(no_VL_events, 'spindles') &&  isempty(no_VL_events.spindles.tstart) ~=1
    
    cfg = [];
    cfg.twin = 0.25;
    cfg.check = 1; % make a plot of 9 random trials to see the quality
    data_spindles_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.spindles, data_remap_FT);
    cfg = [];
    cfg.freq =[7 12];
    data_out.spindles.power = AMPX_get_power(cfg, data_spindles_ft, ExpKeys_remap);
end
close all

% low gamma
cfg = [];
cfg.freq = [40 55];
data_out.lg.power = AMPX_get_power(cfg, data_lg_ft, ExpKeys_remap);

% if max(max(data_out.lg.power.power_distrib_avg)) > 200
%     pause
% end
% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg.power = AMPX_get_power(cfg, data_hg_ft, ExpKeys_remap);

% low random
cfg = [];
cfg.freq = [40 55];
data_out.lg_ran.power = AMPX_get_power(cfg, data_ctrl_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg_ran.power = AMPX_get_power(cfg, data_ctrl_hg_ft, ExpKeys_remap);

data_out_no_VL = data_out;
clear data_out
 
%%addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip')
ft_defaults
% low gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_lg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.low, data_remap_FT);


% high Gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_hg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.high, data_remap_FT);


% low control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_lg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.ctrl_low, data_remap_FT);


% high Control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_hg_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.ctrl_high, data_remap_FT);

% Spindles
if isfield(no_VL_events, 'spindles') &&  isempty(no_VL_events.spindles.tstart) ~=1
    
    cfg = [];
    cfg.twin = 0.25;
    cfg.check = 1; % make a plot of 9 random trials to see the quality
    data_spindles_ft = AMPX_to_FT_trial_split(cfg, no_VL_events.spindles, data_remap_FT);
    cfg = [];
    cfg.freq =[7 12];
    data_out.spindles.power = AMPX_get_power(cfg, data_spindles_ft, ExpKeys_remap);
end
close all

% low gamma
cfg = [];
cfg.freq = [40 55];
data_out.lg.power = AMPX_get_power(cfg, data_lg_ft, ExpKeys_remap);

% if max(max(data_out.lg.power.power_distrib_avg)) > 200
%     pause
% end
% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg.power = AMPX_get_power(cfg, data_hg_ft, ExpKeys_remap);

% low random
cfg = [];
cfg.freq = [40 55];
data_out.lg_ran.power = AMPX_get_power(cfg, data_ctrl_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg_ran.power = AMPX_get_power(cfg, data_ctrl_hg_ft, ExpKeys_remap);

data_out_no_VL = data_out;
clear data_out