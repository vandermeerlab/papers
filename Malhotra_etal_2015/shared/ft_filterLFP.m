function [data,nan_idx] = ft_filterLFP(data,f,varargin)
% function [data,nan_idx] = ft_filterLFP(data,f,varargin)
%
% f should be specified as passband, e.g [1 150]
%
% varargins with defaults:
% ford = 6;
% fmode = 'bandpass'; % 'high', 'low'

ford = 6;
fmode = 'bandpass';
nan_mode = 'zero';
Fs = data.hdr.Fs;
extract_varargin;

if ~exist('B','var')
    
    % build filter
    Wp = f * 2 / Fs; % pass band for filtering
    switch fmode
        case 'bandpass'
            [B,A] = butter(ford,Wp); % builds filter
        case 'highpass'
            [B,A] = butter(ford,Wp(1),'high'); % builds filter
        case 'lowpass'
            [B,A] = butter(ford,Wp(1),'low'); % builds filter
        otherwise
            error('Unknown fmode.');
    end
    
end

% filter
for iS = 1:length(data.trial)
    
    % first check for nans
    d = data.trial{iS};
    
    nan_idx = [];
    if any(isnan(d))
        fprintf('WARNING (ft_filterLFP.m): trial %d has NaNs.\n',iS);
        
        nan_idx = find(isnan(d));
        d(nan_idx) = 0;
        
    end
    
    d = filtfilt(B, A, d); % runs filter
    
    switch nan_mode
        case 'keep'
            d(nan_idx) = NaN;
        case 'interp'
            non_nan_idx = setxor(nan_idx,1:length(d));
            d = interp1(non_nan_idx,d(non_nan_idx),1:length(d));
        case 'zero'
            d(nan_idx) = 0;
    end
    
   data.trial{iS} = d;
   
end

