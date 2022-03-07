function p = ft_getAvgTFRPow(cfg_in,TFR)
% function p = ft_getAvgTFRPow(cfg_in,TFR)
%
% extract average power in a time-frequency tile from a TFR object
%
% may exist already as a ft function
%
% cfg.f = [f1 f2]; % frequencies >= f1 and < f2 are included
% cfg.t = [t1 t2]; % times >= t1 and < t2 are included
% cfg.label = 'label';
%
% MvdM 2015-01-19

cfg = [];
cfg.f = [45 55];
cfg.t = [-1 1];

ProcessConfig;

% do some checks
if cfg.t(1) < TFR.time(1) | cfg.t(2) > TFR.time(end)
   warning('cfg.t [%.2f %.2f] wider than TFR.time [%.2f %.2f]',cfg.t(1),cfg.t(2),TFR.time(1),TFR.time(end));
end

if cfg.f(1) < TFR.freq(1) | cfg.f(2) > TFR.freq(end)
   warning('cfg.t [%.2f %.2f] wider than TFR.time [%.2f %.2f]',cfg.t(1),cfg.t(2),TFR.time(1),TFR.time(end));
end

% set target label
if ~isfield(cfg,'label')
    if length(TFR.label) > 1
        error('No label specified.');
    else
        cfg.label = TFR.label{1};
    end
end

% get idxs to pick up
t_idx = find(TFR.time >= cfg.t(1) & TFR.time < cfg.t(2));
f_idx = find(TFR.freq >= cfg.f(1) & TFR.freq < cfg.f(2));
ch_idx = strmatch(cfg.label,TFR.label);

p = TFR.powspctrm(:,ch_idx,f_idx,t_idx);
p = p(:);

% warn if NaN
if any(isnan(p))
   warning('TFR.powspctrm contains NaNs'); 
end

p = nanmean(p);