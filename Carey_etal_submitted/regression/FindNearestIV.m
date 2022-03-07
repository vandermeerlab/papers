function idx = FindNearestIV(cfg_in,iv_in,t)
% function idx = FindNearestIV(cfg_in,iv_in,t)
%
% find index of iv_in with start nearest to input time t
%
% INPUTS:
%
% iv_in: input IV
% t: input time (must be scalar)
%
% OUTPUTS:
%
% idx: idx

cfg_def = [];
cfg_def.mode = 0; % 0: nearest, -1: nearest preceding, 1: nearest following

cfg = ProcessConfig(cfg_def,cfg_in);

% input checks
if ~isscalar(t)
    error('Input must be scalar format');
end

if ~CheckIV(iv_in)
    error('Input must be IV.');
end

%for iT = length(t):-1:1 % vector form not implemented because then outputs
%cannot be empty
    
    % make some iv starts ineligible depending on mode requested
    iv_t = iv_in.tstart;
    switch cfg.mode
        case -1 % preceding only requested, so remove all following
            iv_t(iv_t > t) = NaN;
        case 1 % following only requested, so remove all preceding
            iv_t(iv_t < t) = NaN;
    end
    
    if all(isnan(iv_t))
        idx = [];
    else
        % obtain idx of nearest iv
        [~,idx] = min(abs(iv_t - t));
    end
%end