function ft_mvdmlab_trialplot(cfg,data,varargin)
% function ft_mvdmlab_trialplot(cfg,data,varargin)
%
% plots multiple trials of raw data in a single plot
% (ft_singleplotER does not do this, taking an average instead)
%
% not yet fieldtripped
%
% MvdM 2013-07-05 initial version
%
% example cfg:
%
%cfg = [];
%cfg.channel = 'R016-2012-10-03-CSC04a';
%cfg.trials = 1:20; 
%cfg.xlim = [-1 3];
%cfg.yscale = [-1 1]; % rescaling of each trial
%cfg.ystep = 2; % spacing between trials
%cfg.lw = 1; % line width
%cfg.plotaverage = 1;
%cfg.plotcolor = [0 0 0];

chan_no = strmatch(cfg.channel,data.label);

nTrials = length(data.trial);
if isfield(cfg,'trials')
    cfg.trials = cfg.trials(cfg.trials < nTrials);
    nTrials = length(cfg.trials);
else
    cfg.trials = 1:nTrials;
end

hold on;

for iT = 1:nTrials
   
    cur_trial = data.trial{iT}(chan_no,:);
    cur_time = data.time{iT};
    
    % rescale
    if length(cfg.yscale) == 2
        cur_trial = rescale(cur_trial,cfg.yscale(1),cfg.yscale(2));
    elseif length(cfg.yscale) == 1
        %cur_trial = (cur_trial-nanmean(cur_trial))./cfg.yscale;
        cur_trial = cur_trial./cfg.yscale;
    else
        error('Unknown yscale.');
    end

    % add correct spacing
    cur_trial = cur_trial + iT*cfg.ystep;
    
    plot(cur_time,cur_trial,'Color',cfg.plotcolor,'LineWidth',cfg.lw);
    
end

% average
if cfg.plotaverage
    av = cell2mat(data.trial);
    av = reshape(av,[length(data.trial{1}) length(data.trial)]);
    av = nanmean(av');
    av = rescale(av,cfg.yscale(1),cfg.yscale(2));
    plot(data.time{1},av,'b');
end

set(gca,'XLim',cfg.xlim,'YTick',(1:nTrials)*cfg.ystep,'YTickLabel',1:nTrials);
title(sprintf('%s',cfg.channel));