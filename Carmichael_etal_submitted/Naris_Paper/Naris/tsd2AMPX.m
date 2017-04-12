function AMPX_out = tsd2AMPX(data_in)
%% tsd2AMPX: converts data in the 'time-stamped data' (tsd) format to the Amplipex format
%
% inputs: 
%     - data_in: [struct] data in the tsd format (must contain tvec, data,
%     hdr, labels)
%
% outputs:
%     - AMPX_out: [struct] data in the AMPX format contains:
%         - type: [str] the data type "AMPX"
%         - tvec: [N x 1] time stamps 
%         - channels {nChan x length data} data for each channel
%         - labels {nChan} string label for each channel
%         - hdr [struct] contains header information such as sampling
%         frequency (Fs)
%
% EC 2016-12-23

AMPX_out.type = 'AMPX';
AMPX_out.hdr = data_in.cfg.hdr{1};
AMPX_out.hdr.Fs = AMPX_out.hdr.SamplingFrequency;

for ichan = size(data_in.data,1):-1:1
    AMPX_out.channels{ichan} = data_in.data(ichan,:);
end
AMPX_out.labels = data_in.label; 
AMPX_out.tvec = data_in.tvec;
