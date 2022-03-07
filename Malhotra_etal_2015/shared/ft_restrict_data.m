function data = ft_restrict_data(cfg,data)

t0 = double(data.hdr.FirstTimeStamp*10^-6); % actual time of first sample
keep = data.time{1}+t0 > cfg.t(1) & data.time{1}+t0 < cfg.t(2);

timestamps_vec = data.hdr.FirstTimeStamp:data.hdr.TimeStampPerSample:data.hdr.LastTimeStamp;

timestamps_vec = timestamps_vec(keep);
data.hdr.FirstTimeStamp = timestamps_vec(1);
data.hdr.LastTimeStamp = timestamps_vec(end);

data.hdr.nSamples = sum(keep);

data.time{1} = data.time{1}(keep);
for iS = 1:length(data.trial)
data.trial{iS} = data.trial{iS}(:,keep);
end

data.sampleinfo = [1 length(data.time{1})];