function [ data ] = Naris_psd( cfg )
%Naris_psd(cfg) loads the data, decimates it (optional) and the computes
%the psds based on some cfg info
% INPUTS
%     cfg.file_name
%     cfg.naris_type = 'pre', left, ipsi...whatever
%     cfg.hann_win  = this is the size of the hanning window.
%     cfg.tetrodes = tetrode channel array from tts_R0XX file
%     cfg.df = decimation factor (default is 10)
%% preprocess the data
data = AMPX_loadData([cfg.fname '-' cfg.naris_type '.dat'], cfg.tetrodes,cfg.df);
% if isfield(cfg, 'butt_ord') == 0
%     cfg.butt_ord = 4; %default butterworth filter order
% end
% if isfield(cfg,'notch') && cfg.notch == 1;
%     [z,p,k] = butter(cfg.butt_ord, [59 61] * 2 / data.hdr.Fs, 'stop'); % note, we ask for 3 outputs instead of 2
%     [sos,g] = zp2sos(z,p,k); % convert to SOS format
%     %     h = dfilt.df2sos(sos,g); % create filter object
%     for iCh = 1:length(data.channels)
%         data_val{iCh} = filtfilt(sos, g, data.channels{iCh});
%     end
%     data.hdr.filter = 'notch';
% else
%     for iCh = 1:length(data.channels)
%         data_val{iCh} = data.channels{iCh};
%     end
% end
%%
% if isfield(cfg, 'whitefilter')
%     white_data = cell(length(data.channels),1);
%     for iChan = 1:length(data.channels)
%         F = sin((data.channels{iChan}.*data.hdr.Fs));
%         white_data{iChan} = F*cos(F.* data.hdr.Fs);
%     end
% end
%% get some psds
% hann_win = round(data.hdr.Fs*cfg.hann_win);
% cfg.hann_win = 1024*8;

Hs = spectrum.welch('Hann',cfg.hann_win,50);
rmpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip\external\signal\');
rmpath('D:\Users\mvdmlab\My_Documents\GitHub\BIOL680\FieldTrip\external\signal\');
for iCh = 1:length(data.channels)
    
    fprintf('psd %d\n',iCh);
    if isfield(cfg, 'whitefilter')
        data.psd{iCh} = psd(Hs,diff(data.channels{iCh}),'Fs',data.hdr.Fs);
    else
        data.psd{iCh} = psd(Hs,data.channels{iCh},'Fs',data.hdr.Fs);
    end
    data.Data{iCh} = data.psd{iCh}.Data.*data.psd{iCh}.Frequencies;
end
end

