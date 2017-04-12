function [naris, cfg] = Naris_fast_mod(cfg_in, data)
%Naris_Fast  Loads the data for s single set of channels (normally a
%tetrode) and then pltos all four session PSDs in a 2x2 plot as well as a
%overlapping PSD plot.  Will save different file extensions for different
%hanning window sizes.
%
%Inputs:
%   cfg [struct]: can be blank if you want to use the defaults.
%% default parameters
cfg_def.chan = 1; % used to pick the first channel in the tetrode. left over from an older version
cfg_def.hann_win_fac = 4;
cfg_def.hann_win = 1024*cfg_def.hann_win_fac;
cfg_def.whitefilter = 'on';
cfg_def.gamma= [45 65; 70 90];
cfg  = ProcessConfig2(cfg_def, cfg_in);

%%
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

   naris.data = data.psd;
    for ich = 1:length(data.labels)
        naris.labels{ich} = data.labels;
        naris.cfg = cfg;
    end



