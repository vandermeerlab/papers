function [z2raw_l, z2raw_h] = Naris_gamma_thresholds(cfg, gamma)

%%

low_gamma_band = cfg.low_gamma;  % default 40-55
high_gamma_band = cfg.high_gamma; % default 70-85
noise_band = [110 225];   % should be in the chewing/scratching range (110-225)
noise_band2 = [225 400];  % should be high ~225-400Hz
channel_of_interest = cfg.chan; % which channel is used for detection
display = 'off';           % plots the events and lets your see each event

%% merge pre and post

gamma.control.data.channels{1} = [gamma.pre.data.channels{1}' , gamma.post.data.channels{1}'];
gamma.control.data.tvec = [gamma.pre.data.tvec', ((gamma.pre.data.tvec(end)+1/gamma.pre.data.hdr.Fs)+gamma.post.data.tvec)'];
gamma.control.data.hdr = gamma.pre.data.hdr;

%%
data = gamma.control.data;
data_std = std(data.channels{channel_of_interest});

data.channels{channel_of_interest}(data.channels{channel_of_interest} > 3*data_std) =0;
data.channels{channel_of_interest}(data.channels{channel_of_interest} < -3*data_std) =0;

%% prepare data for event detection
cfg = [];
data.channels{channel_of_interest} = data.channels{channel_of_interest}';
data.cgf = cfg;
data_tsd = tsd(data.tvec,data.channels{channel_of_interest});
if size(data_tsd.tvec,2) ==1
    data_tsd.tvec = data_tsd.tvec'; % transpose the data (not the tvec) or else you will get an erroe in the
    disp('data.tvec needs to be transposed')
end
if size(data_tsd.data,2) ==1
    data_tsd.data = data_tsd.data'; % transpose the data (not the tvec) or else you will get an erroe in the
    disp('data.tvec needs to be transposed')
end
data_tsd.cfg.hdr{1}.SamplingFrequency = data.hdr.Fs;
data_tsd.data(isnan(data_tsd.data)) = 0;
data_tsd.label = num2str(channel_of_interest);
cfg = []; cfg.decimateFactor = 1; % reduce sampling frequency to speed things up
csc = decimate_tsd(cfg,data_tsd);


%% low gamma detection
cfg = []; cfg.f = low_gamma_band;
cfg.threshold = 2;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
z2raw_l = (mean(cscF.data) +cfg.threshold*std(cscF.data));
fprintf(['Raw low gamma detection threshold is ' num2str(z2raw_l) '\n'])

plot(cscF.data); hold on; hline(z2raw_l, 'g', 'g')

%% high gamma detection
cfg = []; cfg.f = high_gamma_band;
cfg.threshold = 2;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
z2raw_h = (mean(cscF.data) +cfg.threshold*std(cscF.data));
fprintf(['Raw high gamma detection threshold is ' num2str(z2raw_h) '\n'])
clear gamma
%return these values
