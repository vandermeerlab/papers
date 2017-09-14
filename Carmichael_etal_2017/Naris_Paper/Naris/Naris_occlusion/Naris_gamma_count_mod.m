function [cfg, Naris_stats_out, naris_gamma] = Naris_gamma_count_mod(cfg_in)
%%Naris_gamma_detect : detects all the gamma events within a file and
% returns a trialified structure of all the events required for
% Naris_gamma_stats which gives all the detected gamma events.
%
cfg_def.fname = mkfile;
cfg_def.df = 10;
cfg_def.low_gamma= [45 65];
cfg_def.high_gamma = [70 90];
cfg = ProcessConfig2(cfg_def, cfg_in);

%% set up defualts:
ExpKeys = cfg.ExpKeys;
if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R060')
    Naris_exp = {'control', 'pre', 'right', 'left', 'post'};
else
    Naris_exp = {'control', 'pre', 'ipsi', 'contra', 'post'};
end

cfg.cd = cd;
% this should be carried over from the "Naris_fast" function in the Naris
% oclussion sandbox


%% load data

for exp = 2:5
    cfg.naris_type = Naris_exp{exp};
    if isunix
        cd([cfg.data_path '/' cfg.fname(1:4) '/' cfg.fname '/' cfg.fname '_' cfg.naris_type])
    else
        cd([cfg.data_path '\' cfg.fname(1:4) '\' cfg.fname '\' cfg.fname '_' cfg.naris_type])
    end
    if length(dir([ cd '\' '*.dat'])) >0 % workaround for finding files using a wildcard
        Naris.(cfg.naris_type).data = AMPX_loadData([cfg.fname '-' cfg.naris_type '.dat'], cfg.chan_to_load,cfg.df);
    elseif length(dir([ cd '\' '*.ncs'])) >0% identify if the recordings are using Neuralynx
        cell_idx = strfind(ExpKeys.Chan_to_use_labels, 'NAc');
        chan_idx = find(not(cellfun('isempty', cell_idx)));
        cfg_data.fc = ExpKeys.Chan_to_use(chan_idx);
        cfg_data.resample = 2000; 
        data = LoadCSC(cfg_data);
        evt = LoadEvents([]);
        % this section is to restrict to only the pot sections.
        idx = strfind(evt.label, 'Starting Recording');
        start_idx = find(not(cellfun('isempty', idx)));
        idx = strfind(evt.label, 'Stopping Recording');
        stop_idx = find(not(cellfun('isempty', idx)));
        
        % check to make sure NLX didn't mess up the start stop events (this happened for R104-2016-09-26_ipsi. It has 19 start times for some reason)
        if length(evt.t{start_idx}) >= 3 % should only have 2
            [~, trk_idx]  = max(diff(evt.t{start_idx}(2:end))); % find the largest gap
            trk_idx = trk_idx +1; % offset by one to compensate for the "diff"
        else
            trk_idx = [];
        end
        data_r = restrict(data,evt.t{start_idx}(1),evt.t{stop_idx}(1));
        Naris.(cfg.naris_type).data = tsd2AMPX(data_r);
        disp(['Sampling Frequency: ' num2str(Naris.(cfg.naris_type).data.hdr.Fs)])
    end
    for iChan = size(Naris.(cfg.naris_type).data.channels,2):-1:1
        data_temp{1,iChan} = Naris.(cfg.naris_type).data.channels{1,iChan} - mean(Naris.(cfg.naris_type).data.channels{1,iChan});
    end
    Naris.(cfg.naris_type).data.channels = data_temp;
    Naris.(cfg.naris_type).data.dc_remove = 'yes';
    clear data_temp;
    Naris.(cfg.naris_type).data = AMPX_to_tsd(Naris.(cfg.naris_type).data);
end

%% merge the pre and post sessions to create a 'control' used for
% extracting the raw gamma thresholds used for each phase
% modify the
Naris.post.data.tvec =  Naris.post.data.tvec + (Naris.pre.data.tvec(end)+(mode(diff(Naris.post.data.tvec))));
Naris.control.data = UnionTSD([], Naris.pre.data, Naris.post.data);
Naris.control.data.cfg.hdr{1}.SamplingFrequency = Naris.pre.data.cfg.hdr{1}.SamplingFrequency;
%% detect gamma events across each data_type
fprintf('\n\n=====================\nGamma Detection\n=====================\n')
LoadExpKeys
cfg_thresh.f_bandpass = {cfg.low_gamma,cfg.high_gamma,cfg.low_gamma, cfg.high_gamma}; 
for iExp = 1:length(Naris_exp)
    if strcmp(Naris_exp{iExp}, 'control') ==1
        [evts.(Naris_exp{iExp}), ~, evt_thr.control] = JK_Julien_DetectEvents_thresholds(cfg_thresh, Naris.(Naris_exp{iExp}).data, ExpKeys);
    else
        [evts.(Naris_exp{iExp})] = AMPX_Julien_DetectEvents_naris(Naris.(Naris_exp{iExp}).data, ExpKeys, 'PARAM_detect_method', 'raw','PARAM_detect_thr', evt_thr.control, 'PARAM_f_bandpass', cfg_thresh.f_bandpass);
        fprintf(['low events:  ' num2str(length(evts.(Naris_exp{iExp}).low.tstart)) '\n']);
        fprintf(['high events: ' num2str(length(evts.(Naris_exp{iExp}).high.tstart)) '\n']);
    end
end
%% correct for the length of the recording
for iExp = 1:length(Naris_exp)
    if strcmp(Naris_exp{iExp}, 'control')
        pre = length(evts.(Naris_exp{iExp}).low.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
        post = length(evts.(Naris_exp{iExp}).low.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
        count_L(iExp) = median([pre, post]);
        pre = length(evts.(Naris_exp{iExp}).high.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
        post = length(evts.(Naris_exp{iExp}).high.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
        count_H(iExp) = median([pre, post]);
    else
        count_L(iExp) = length(evts.(Naris_exp{iExp}).low.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
        count_H(iExp) = length(evts.(Naris_exp{iExp}).high.tstart)/((length(Naris.(Naris_exp{iExp}).data.data)/Naris.(Naris_exp{iExp}).data.cfg.hdr{1}.SamplingFrequency)/60);
    end
end

for iExp = 2:5
    cfg.naris_type = Naris_exp{iExp};
    fprintf(['_________\n' cfg.naris_type '\n'])
    fprintf(['low events:  ' num2str(count_L(iExp)) '\n']);
    fprintf(['high events: ' num2str(count_H(iExp)) '\n']);
end
Naris_stats_out.labels = Naris_exp; 
Naris_stats_out.low.count =  [length(evts.(Naris_exp{1}).low.tstart), length(evts.(Naris_exp{2}).low.tstart), length(evts.(Naris_exp{3}).low.tstart), length(evts.(Naris_exp{4}).low.tstart), length(evts.(Naris_exp{5}).low.tstart)];
Naris_stats_out.high.count =  [length(evts.(Naris_exp{1}).high.tstart), length(evts.(Naris_exp{2}).high.tstart), length(evts.(Naris_exp{3}).high.tstart), length(evts.(Naris_exp{4}).high.tstart), length(evts.(Naris_exp{5}).high.tstart)];
Naris_stats_out.low.rate = count_L;
Naris_stats_out.high.rate = count_H;
%% add it to the master gamma events file:
gamma_current = evts;  clear Naris
load('G:\Naris\Paper_naris_gamma_dec.mat');
naris_gamma.(strrep(cfg.fname, '-', '_')) = gamma_current;
naris_gamma.(strrep(cfg.fname, '-', '_')).count_L = count_L;
naris_gamma.(strrep(cfg.fname, '-', '_')).count_H = count_H;
save('G:\Naris\Paper_naris_gamma_dec.mat', 'naris_gamma', '-v7.3');
