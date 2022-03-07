clear all;

% load preselected sessions
%cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
%load fd_isidro

% load all promoted sessions
cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-20');
load all_files
fd = fd'

%%
for iFD =1:1:length(fd)
% for iFD =1:2
  pushdir(fd{iFD});
  fd{iFD} = FindFile('*Events.nev','CheckSubdirs',1);
  %fd{iFD} = num2str(cell2mat(fd{iFD}));
end

fd=fd';

% fd = FindFiles('*Events.nev','CheckSubdirs',1);
fsz = 16;
writeOutput = 0;

%%
% for each candidate, get rat id, date, and type of session (value/risk)
clear ratID yyyy mm dd type
for iFD = length(fd):-1:1
   
    [fp,fn,fe] = fileparts(fd{iFD});
    
    fp_sep = regexp(fp,'\'); % note this is untested on non Windows machines
    session_id = fp(fp_sep(end)+1:end);
    
    ratID(iFD) = str2num(session_id(2:4));
    yyyy(iFD) = str2num(session_id(6:9));
    mm(iFD) = str2num(session_id(11:12));
    dd(iFD) = str2num(session_id(14:15));
    type{iFD} = session_id(17:end); % this is only for data with separate value and risk folders
    
end

% mapping from outcomes (stored in events file) to cue identity (for VALUE)
%what_cueV = {'c1','c3','c5'};
cue_to_toneV(7,:) = [2 5 1]; % means: 1 pellet cue is tone 2, 3 pellet cue is tone 5, 5 pellet cue is tone 1 (for rat 7)
cue_to_toneV(11,:) = [3 2 4];
cue_to_toneV(12,:) = [4 1 5];
cue_to_toneV(16,:) = [3 2 4];
cue_to_toneV(18,:) = [2 5 1];
cue_to_toneV(20,:) = [5 4 2];
cue_to_toneV(14,:) = [4 1 5];

% mapping from outcomes (stored in events file) to cue identity (for RISK)
%what_cueR = {'norisk','lorisk','hirisk'};
cue_to_toneR(7,:) = [5 4 3]; % means: no risk cue is tone 5, low risk is tone 4, high risk tone 3
cue_to_toneR(11,:) = [2 1 5];
cue_to_toneR(12,:) = [1 3 2];
cue_to_toneR(16,:) = [2 1 5];
cue_to_toneR(18,:) = [5 4 3];
cue_to_toneR(20,:) = [4 3 1];

%% find sessions to use -- rats
goodRats = [14 16 18 20];
keep = [];
for iRat = 1:length(goodRats)
   
    keep = cat(2,keep,find(ratID == goodRats(iRat)));
    
end

ratID = ratID(keep); mm = mm(keep); dd = dd(keep); type = type(keep); fd = fd(keep);

% days
%keep = find(dd == 8);

%ratID = ratID(keep); mm = mm(keep); dd = dd(keep); type = type(keep); fd = fd(keep);

%% go through sessions and get data

ALL_trial = []; ALL_rat = []; ALL_day = []; ALL_outcome = []; ALL_cue = []; ALL_type = [];index=1;fd_new = cell(length(fd),1);indx = 1;
for iFD = length(fd):-1:1
   
   %% io
   [fp,fn,fe] = fileparts(fd{iFD});
   
   session_id = fp(fp_sep(end)+1:end);
   
   cd(fp);
   session_id = fp(fp_sep(end)+1:end);
   fprintf('=== processing session %s ===\n',fp(end-14:end));
   
   % unzip vt file if not available
   if isempty(FindFiles('*.nvt'))
       zf = FindFile('*.zip');
       unzip(zf);
   end
    
   % load vt data
   [VTTimeStamps, x, y] = Nlx2MatVT('VT1.nvt',[1 1 1 0 0 0],0,1,[]); VTTimeStamps = VTTimeStamps * 10^-6;
   goodSamples = find(x ~= 0 & y ~= 0 & ~isnan(x) & ~isnan(y));
   t = VTTimeStamps(goodSamples);
   x = x(goodSamples);
   y = y(goodSamples);
   x = tsd(t',x'); y = tsd(t',y');
   
   s = GetLinSpd(x,y);
   
   %% get value events
  
   run(FindFile('*keys.m'));
   %trialfun_validate;
   
   cfg = [];
   cfg.trialfun = 'ft_trialfun_lineartracktone2';
   cfg.trialdef.pre = 2.5; cfg.trialdef.post = 5;
   
   cfg.trialdef.eventtype = 'cue'; % could be 'nosepoke', 'reward', 'cue'
   cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
   cfg.trialdef.block = 'value'; % could be 'value', 'risk', 'both'
   cfg.ExpKeys = ExpKeys;
   
   
   cfg.trialdef.cue = {'c1'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.c1 = evt.ts;
   
   cfg.trialdef.cue = {'c3'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.c3 = evt.ts;
   events.c3v = events.c3(events.c3 < ExpKeys.TimeOffTrack(1));
   events.c3r = events.c3(events.c3 > ExpKeys.TimeOffTrack(1));
   
   cfg.trialdef.cue = {'c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.c5 = evt.ts;
   
   cfg.trialdef.block = 'risk'; % could be 'value', 'risk', 'both'
   
   cfg.trialdef.cue = {'lo'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.lo = evt.ts;
   
   cfg.trialdef.cue = {'hi'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.hi = evt.ts;
   
   
   %%
   this_type = 0; % value is 0, risk 1
   
   switch this_type
       
       case 0
           
           what_cueV = {'c1','c3v','c5'};
           for iW = 1:length(what_cueV)
               
               [P,X0D,SD] = tsdPETH2(s, events.(what_cueV{iW}),'dt',0.1);
               ALL_tvec = Range(P);
               
               nTrials = size(X0D,1);
               
               ALL_trial = cat(1,ALL_trial,X0D);
               ALL_rat = cat(1,ALL_rat,repmat(ratID(iFD),[nTrials 1]));
               ALL_day = cat(1,ALL_day,repmat(dd(iFD),[nTrials 1]));
               ALL_type = cat(1,ALL_type,repmat(this_type,[nTrials 1]));
               ALL_outcome = cat(1,ALL_outcome,repmat(iW,[nTrials 1]));
               ALL_cue = cat(1,ALL_cue,repmat(cue_to_toneV(ratID(iFD),iW),[nTrials 1]));
               
           end
           
       case 1
           
           what_cueR = {'lo','c3r','hi'};
           for iW = 1:length(what_cueR)
               
               [P,X0D,SD] = tsdPETH2(s, events.(what_cueR{iW}),'dt',0.1);
               ALL_tvec = Range(P);
               
               nTrials = size(X0D,1);
               
               ALL_trial = cat(1,ALL_trial,X0D);
               ALL_rat = cat(1,ALL_rat,repmat(ratID(iFD),[nTrials 1]));
               ALL_day = cat(1,ALL_day,repmat(dd(iFD),[nTrials 1]));
               ALL_type = cat(1,ALL_type,repmat(this_type,[nTrials 1]));
               ALL_outcome = cat(1,ALL_outcome,repmat(iW,[nTrials 1]));
               ALL_cue = cat(1,ALL_cue,repmat(cue_to_toneR(ratID(iFD),iW),[nTrials 1]));
               
           end
           
   end
   
end

%% z-score speed data
ALL_trialZ = (ALL_trial-repmat(nanmean(ALL_trial,2),[1 size(ALL_trial,2)]))./repmat(nanstd(ALL_trial,[],2),[1 size(ALL_trial,2)]);

%% create peak speed output
out.peakSpeed = max(ALL_trial,[],2);
out.peakSpeedZ = max(ALL_trialZ,[],2);
t1 = 0; t2 = 1; t_idx = find(ALL_tvec > t1 & ALL_tvec <= t2);
out.avSpeed = nanmean(ALL_trial(:,t_idx),2);
out.avSpeedZ = nanmean(ALL_trialZ(:,t_idx),2);
out.cue = ALL_cue;
out.day = ALL_day;
out.rat = ALL_rat;
out.outcome = ALL_outcome;


%%
%%%%%%%%%%%%%%%%%%%
%%% VALUE BLOCK %%%
%%%%%%%%%%%%%%%%%%%

%value outcome plot (all days, all rats)

what = [1 2 3]; cols = 'rgb'; clear h1 h2;
nTrials = 0;
for iW = 1:length(what)
    
    keep = find(ALL_type == this_type & ALL_outcome == what(iW));
    
    nTrials = nTrials + length(keep);
    
    temp_data = ALL_trial(keep,:);
    temp_dataZ = ALL_trialZ(keep,:);
    
    figure(1); % raw
    subplot(2,4,[1 2 5 6]);
    h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
    errorbar(ALL_tvec,nanmean(temp_data),nanstd(temp_data)./sqrt(length(fd)),'Color',cols(iW)); hold on; % added dec 5
    figure(2); % z
    subplot(2,4,[1 2 5 6]);
    h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
    errorbar(ALL_tvec,nanmean(temp_dataZ),nanstd(temp_dataZ)./sqrt(length(fd)),'Color',cols(iW)); hold on;  %added dec 5
end
figure(1);
legend(h1,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',20,'XLim',[-2 5],'XTick',-2:5); xlabel('time (s)'); ylabel('speed (cm/s)');
title(sprintf('all rats, all days (%d trials)',nTrials));
vl = plot([0 0],ylim,'k:','LineWidth',2);

figure(2);
legend(h2,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',20,'XLim',[-2 5],'XTick',-2:5); xlabel('time (s)'); ylabel('speed (z-scored)');
title(sprintf('all rats, all days (%d trials)',nTrials));
vl = plot([0 0],ylim,'k:','LineWidth',2);


%% by rat
rhat = [14 16 18 20];
clear h1 h2;
subplots = [3 4 7 8];

for iR = 1:length(rhat)
    what = [1 2 3]; cols = 'rgb';
    nTrials = 0;
    for iW = 1:length(what)
        
        keep = find(ALL_type == 0 & ALL_outcome == what(iW) & ALL_rat == rhat(iR));
        
        nTrials = nTrials + length(keep);
        
        temp_data = ALL_trial(keep,:);
        temp_dataZ = ALL_trialZ(keep,:);
        
        figure(1); % raw
        subplot(2,4,subplots(iR));
        h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
        vl = plot([0 0],[-1 2.5],'k:','LineWidth',2);
        
        figure(2); % z
        subplot(2,4,subplots(iR));
        h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
        vl = plot([0 0],[-1 2.5],'k:','LineWidth',2);
        
    end
    figure(1);
    %legend(h1,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
    set(gca,'FontSize',20,'XLim',[-2 5],'YLim',[-1 2.5],'XTick',0); 
    xlabel('time (s)'); ylabel('speed (cm/s)');
    title(sprintf('R0%d (%d t)',rhat(iR),nTrials));
    
    figure(2);
    %legend(h2,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
    set(gca,'FontSize',20,'XLim',[-2 5],'YLim',[-1 2.5],'XTick',0); 
    xlabel('time (s)'); ylabel('speed (z-scored)');
    title(sprintf('R0%d (%d t)',rhat(iR),nTrials));
end

%% stats
idx1 = find(out.outcome == 1);
idx5 = find(out.outcome == 3);

fprintf('peakSpeedZ 1p (%.2f +/- %.2f) v 5p (%.2f +/- %.2f): p = %3.2e (Wilcoxon rank-sum)\n',nanmean(out.peakSpeedZ(idx1)),nanstd(out.peakSpeedZ(idx1))./2,nanmean(out.peakSpeedZ(idx5)),nanstd(out.peakSpeedZ(idx5))./2,ranksum(out.peakSpeedZ(idx1),out.peakSpeedZ(idx5)));
fprintf('AvSpeedZ 1p (%.2f +/- %.2f) v 5p (%.2f +/- %.2f): p = %3.2e (Wilcoxon rank-sum)\n',nanmean(out.avSpeedZ(idx1)),nanstd(out.avSpeedZ(idx1))./2,nanmean(out.avSpeedZ(idx5)),nanstd(out.avSpeedZ(idx5))./2,ranksum(out.avSpeedZ(idx1),out.avSpeedZ(idx5)));

%% cue plot (all days, all rats) for value trials
what = 1:5; cols = 'rgbcm';  clear h1 h2;
nTrials = 0;
for iW = 1:length(what)
    
    keep = find(ALL_type == 0 & ALL_cue == what(iW));
    
    nTrials = nTrials + length(keep);
    
    temp_data = ALL_trial(keep,:);
    temp_dataZ = ALL_trialZ(keep,:);
    
    figure(3); % raw
    subplot(2,4,[1 2 5 6]);
    h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
    
    figure(4); % z
    subplot(2,4,[1 2 5 6]);
    h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
    
end
figure(3);
legend(h1,{'s1','s2','s3','s4','s5'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',20,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (cm/s)');
title(sprintf('(V) all rats, all days (%d trials)',nTrials));

figure(4);
legend(h2,{'s1','s2','s3','s4','s5'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',20,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (z-scored)');
title(sprintf('(V) all rats, all days (%d trials)',nTrials));

%%

%%%%%%%%%%%%%%%%%%
%%% RISK BLOCK %%%
%%%%%%%%%%%%%%%%%%

% risk outcome plot (all days, all rats)
what = [1 2 3]; cols = 'rgb'; clear h1 h2;
nTrials = 0;
for iW = 1:length(what)
    
    keep = find(ALL_type == 1 & ALL_outcome == what(iW));
    
    nTrials = nTrials + length(keep);
    
    temp_data = ALL_trial(keep,:);
    temp_dataZ = ALL_trialZ(keep,:);
    
    figure(5); % raw
    subplot(2,4,[1 2 5 6]);
    h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
    
    errorbar(ALL_tvec,nanmean(temp_data),nanstd(temp_data)./sqrt(15),'Color',cols(iW)); hold on;
    figure(6); % z
    subplot(2,4,[1 2 5 6]);
    h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
    errorbar(ALL_tvec,nanmean(temp_dataZ),nanstd(temp_dataZ)./sqrt(15),'Color',cols(iW)); hold on;
end
figure(5);
legend(h1,{'lorisk','norisk','hirisk'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',10,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (cm/s)');
title(sprintf('all rats, all days (%d trials)',nTrials));

figure(6);
legend(h2,{'lorisk','norisk','hirisk'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',10,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (z-scored)');
title(sprintf('all rats, all days (%d trials)',nTrials));

%% by rat
rhat = [14 16 18 20];  clear h1 h2;
subplots = [3 4 7 8];
for iR = 1:length(rhat)
    what = [1 3]; cols = 'rgb';
    nTrials = 0;
    for iW = 1:length(what)
        
        keep = find(ALL_type == 1 & ALL_outcome == what(iW) & ALL_rat == rhat(iR));
        
        nTrials = nTrials + length(keep);
        
        temp_data = ALL_trial(keep,:);
        temp_dataZ = ALL_trialZ(keep,:);
        
        figure(5); % raw
        subplot(2,4,subplots(iR));
        h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
        
        figure(6); % z
        subplot(2,4,subplots(iR));
        h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
        
    end
    figure(5);
    legend(h1,{'lorisk','hirisk'},'Location','Northeast'); legend boxoff;
    set(gca,'FontSize',10,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (cm/s)');
    title(sprintf('rat %d, all days (%d trials)',rhat(iR),nTrials));
    
    figure(6);
    legend(h2,{'lorisk','hirisk'},'Location','Northeast'); legend boxoff;
    set(gca,'FontSize',10,'XLim',[-2 5]); xlabel('time (s)'); ylabel('speed (z-scored)');
    title(sprintf('rat %d, all days (%d trials)',rhat(iR),nTrials));
end



