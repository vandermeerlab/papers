function [cfg, Naris] = Naris_gamma_count(cfg)
%%Naris_gamma_detect : detects all the gamma events within a file and
% returns a trialified structure of all the events required for
% Naris_gamma_stats which gives all the detected gamma events.
%


%% set up defualts:
cfg.fname = mkfile;
LoadExpKeys
[cfg] = Naris_cfgs(cfg);
if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R060')
    Naris_exp = {'control', 'pre', 'right', 'left', 'post'};
else
    Naris_exp = {'control', 'pre', 'ipsi', 'contra', 'post'};
end
cfg.df = 10;
cfg.low_gamma= [40 55];
cfg.high_gamma = [70 85];
if ~isfield(cfg, 'tetrodes')
    error('Channel to proces not assigned.  Should be determine during the Naris power processing')
end
% this should be carried over from the "Naris_fast" function in the Naris
% oclussion sandbox


cfg.chan = 1;
%% load data

for exp = 2:5
    cfg.naris_type = Naris_exp{exp};3
    if exist('*.dat')
        Naris.(cfg.naris_type).data = AMPX_loadData([cfg.fname '-' cfg.naris_type '.dat'], cfg.tetrodes,cfg.df);
    elseif exist('*.ncs') % identify if the recordings are using Neuralynx
        cell_idx = strfind(ExpKeys.Chan_to_use_labels, 'NAc');
        chan_idx = find(not(cellfun('isempty', cell_idx)));
        cfg_data.fc = ExpKeys.Chan_to_use(chan_idx);
        data = LoadCSC(cfg_data);
        Naris.(cfg.naris_type).data = tsd2AMPX(data);
        disp(['Sampling Frequency: ' num2str(data_ampx.hdr.Fs)])
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
for iExp = 1:length(Naris_exp)
    if strcmp(Naris_exp{iExp}, 'control') ==1
        [evts.(Naris_exp{iExp}), ~, evt_thr.control] = JK_Julien_DetectEvents_thresholds([], Naris.(Naris_exp{iExp}).data, ExpKeys);
    else
        [evts.(Naris_exp{iExp})] = AMPX_Julien_DetectEvents_naris(Naris.(Naris_exp{iExp}).data, ExpKeys, 'PARAM_detect_method', 'raw','PARAM_detect_thr', evt_thr.control);
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
%% add it to the master gamma events file:
gamma_current = evts;  clear Naris
load('G:\Naris\Paper_naris_gamma_dec.mat');
naris_gamma.(strrep(cfg.fname, '-', '_')) = gamma_current;
naris_gamma.(strrep(cfg.fname, '-', '_')).count_L = count_L;
naris_gamma.(strrep(cfg.fname, '-', '_')).count_H = count_H;
save('G:\Naris\Paper_naris_gamma_dec.mat', 'naris_gamma', '-v7.3');
