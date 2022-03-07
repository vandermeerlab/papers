function out = cueSpeed_singleSession(varargin)

fsz = 16;
writeOutput = 0;

extract_varargin;

run(FindFile('*keys.m'));
fd{1} = FindFile('*.Nev');

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
what_cueV = {'c1','c3','c5'};
cue_to_toneV(7,:) = [2 5 1]; % means: 1 pellet cue is tone 2, 3 pellet cue is tone 5, 5 pellet cue is tone 1 (for rat 7)
cue_to_toneV(11,:) = [3 2 4];
cue_to_toneV(12,:) = [4 1 5];
cue_to_toneV(16,:) = [3 2 4];
cue_to_toneV(18,:) = [2 5 1];
cue_to_toneV(20,:) = [5 4 2];

% mapping from outcomes (stored in events file) to cue identity (for RISK)
what_cueR = {'norisk','lorisk','hirisk'};
cue_to_toneR(7,:) = [5 4 3]; % means: no risk cue is tone 5, low risk is tone 4, high risk tone 3
cue_to_toneR(11,:) = [2 1 5];
cue_to_toneR(12,:) = [1 3 2];
cue_to_toneR(16,:) = [2 1 5];
cue_to_toneR(18,:) = [5 4 3];
cue_to_toneR(20,:) = [4 3 1];

%% find sessions to use -- rats
goodRats = [14];
keep = [];
for iRat = 1:length(goodRats)
   
    keep = cat(2,keep,find(ratID == goodRats(iRat)));
    
end

ratID = ratID(keep); mm = mm(keep); dd = dd(keep); type = type(keep); fd = fd(keep);

% days
%keep = find(dd == 8);

%ratID = ratID(keep); mm = mm(keep); dd = dd(keep); type = type(keep); fd = fd(keep);

%% go through sessions and get data

ALL_trial = []; ALL_rat = []; ALL_day = []; ALL_outcome = []; ALL_cue = []; ALL_type = [];
for iFD = length(fd):-1:1
    
   %% io
   [fp,fn,fe] = fileparts(fd{iFD});
   cd(fp);
   
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
   
   cfg.trialdef.cue = {'c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
   [trl, evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
   events.c5 = evt.ts;
   
   %%
   this_type = 0; % value is 0, risk 1
   
   what_cueV = {'c1','c3','c5'};
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
   
end

%% z-score speed data
ALL_trialZ = (ALL_trial-repmat(nanmean(ALL_trial,2),[1 size(ALL_trial,2)]))./repmat(nanstd(ALL_trial,[],2),[1 size(ALL_trial,2)]);

%% create peak speed output
out.peakSpeed = max(ALL_trial,[],2);
out.peakSpeedZ = max(ALL_trialZ,[],2);
out.cue = ALL_cue;
out.day = ALL_day;
out.rat = ALL_rat;

%%
%%%%%%%%%%%%%%%%%%%
%%% VALUE BLOCK %%%
%%%%%%%%%%%%%%%%%%%

%value outcome plot (all days, all rats)

what = [1 2 3]; cols = 'rgb'; clear h1 h2;
nTrials = 0;
for iW = 1:length(what)
   
    keep = find(ALL_type == 0 & ALL_outcome == what(iW));
    
    nTrials = nTrials + length(keep);
    
    temp_data = ALL_trial(keep,:);
    temp_dataZ = ALL_trialZ(keep,:);
  
    figure(1); % raw
    subplot(2,4,[1 2 5 6]);
    h1(iW) = plot(ALL_tvec,nanmean(temp_data),'Color',cols(iW),'LineWidth',2); hold on;
    
    subplot(2,4,[3 4 7 8]);
    h2(iW) = plot(ALL_tvec,nanmean(temp_dataZ),'Color',cols(iW),'LineWidth',2); hold on;
    
end
subplot(2,4,[1 2 5 6]);
legend(h1,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',fsz,'XLim',[-2 5]); xlabel('cue time (s)'); ylabel('speed (pix/s??)');
th = title(sprintf('%s (%d trials)',session_id,nTrials)); set(th,'Interpreter','none');

subplot(2,4,[3 4 7 8]);
legend(h2,{'1p','3p','5p'},'Location','Northeast'); legend boxoff;
set(gca,'FontSize',fsz,'XLim',[-2 5]); xlabel('cue time (s)'); ylabel('z-scored speed');
th = title(sprintf('%s (%d trials)',session_id,nTrials)); set(th,'Interpreter','none');

if writeOutput
   
    fout = cat(2,pwd,'/images/',session_id,'_cueSpeed.png')
    print(gcf,'-dpng','-r300',fout);
    
end

pause;

