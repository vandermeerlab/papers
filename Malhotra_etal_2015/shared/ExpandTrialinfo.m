function ExpandTrialinfo(cfg_in)
% function ExpandTrialinfo(cfg_in)
%
% populate existing *trialinfo.mat skeleton with data fields

cfg = [];
cfg.writeOutput = 1;
cfg.doBaseline = 0; % is slow
ProcessConfig;

%% load
run(FindFile('*keys.m'));

fc = ExpKeys.goodGamma(1);
data = ft_read_neuralynx_interp(fc);
data.label{1} = 'vStr';

load(FindFile('*trialinfo.mat'));
load(FindFile('*vt.mat'))

%% get speed
spd = GetLinSpd(x,y);

%% z-score speed (all speeds are now zscored, 2015-08-05)
spdR = Range(spd); spdD = Data(spd);
spdZ = tsd(spdR,(spdD-nanmean(spdD)./nanstd(spdD)));

%% filter -- this is important to fill in NaNs in the data
data = ft_filterLFP(data,200,'fmode','lowpass','nan_mode','interp');

%% define time windows of interest (to loop over later)
DEF_cueWindowLabel = {'preCue','postCue'};
DEF_cueWindowTime = {[-0.5 0],[0 0.5]};

DEF_npWindowLabel = {'preNP','postNP'};
DEF_npWindowTime = {[-0.5 0],[0 1.25]};

DEF_freqBandLabel = {'g50','g80'};
DEF_freqBandFreq = {[50 65],[70 85]};

%% define subject ordering
DEF_ratID = {'R014','R016','R018','R020'};

%% define trialification cfg's (twin's should be bigger than windows of interest above)
cfg_cue = [];
cfg_cue.twin = [-1 1];
cfg_cue.mode = 'nlx';
cfg_cue.hdr = data.hdr;

cfg_np = [];
cfg_np.twin = [-1.75 1.75];
cfg_np.mode = 'nlx';
cfg_np.hdr = data.hdr;

%% define TFR cfg (other than cfg.toi)
cfg_TFR              = []; % start with empty cfg
cfg_TFR.output       = 'pow';
cfg_TFR.keeptrials   = 'yes'; % need this for stats later
cfg_TFR.channel      = {'vStr'};
cfg_TFR.method       = 'mtmconvol';
cfg_TFR.taper        = 'hanning';
cfg_TFR.foi          = 25:1:110; % frequencies of interest
cfg_TFR.t_ftimwin    = 19*(1./cfg_TFR.foi);

%% loop over trials
nTrials = length(trialinfo.cue_time);

for iT = 1:nTrials

    fprintf('** Processing trial %d/%d...\n',iT,nTrials);
    
    % first, make cue data
    cfg_cue.t = trialinfo.cue_time(iT) + ExpKeys.TimeOnTrack(1); % matching MakeTrialinfo where this is subtracted
    trl_cue = ft_maketrl(cfg_cue);
    
    cfg_temp = [];
    cfg_temp.trl = trl_cue;
    data_cue = ft_redefinetrial(cfg_temp,data);
    
    % TFR
    cfg_TFR.toi          = cfg_cue.twin(1):0.01:cfg_cue.twin(2);
    TFR = ft_freqanalysis(cfg_TFR, data_cue);
    
    % go through cue windows and frequencies to process
    nCues = length(DEF_cueWindowLabel);
    for iCueWin = 1:nCues
        
        cfg.temp = [];
        cfg_temp.t = DEF_cueWindowTime{iCueWin};

        nFreqs = length(DEF_freqBandLabel);
        for iFreq = 1:nFreqs
            
            cfg_temp.f = DEF_freqBandFreq{iFreq};
            p = ft_getAvgTFRPow(cfg_temp,TFR);
            
            if(isnan(p))
                error('NaN found.');
            end
            
            % write output
            fn = cat(2,DEF_cueWindowLabel{iCueWin},'_',DEF_freqBandLabel{iFreq});
            trialinfo.(fn)(iT) = p;

        end
        
        % speed
        t0 = cfg_cue.t + cfg_temp.t(1);
        t1 = cfg_cue.t + cfg_temp.t(2);
        fn = cat(2,DEF_cueWindowLabel{iCueWin},'_spd');
        trialinfo.(fn)(iT) = nanmean(Data(Restrict(spdZ,t0,t1)));
        
        fn = cat(2,DEF_cueWindowLabel{iCueWin},'_spdmax'); % added 2015-05-04
        trialinfo.(fn)(iT) = max(Data(Restrict(spdZ,t0,t1)));
    
    end
    
    % next, make nosepoke data
    if iT > length(trialinfo.pb_time) % no nosepoke recorded, enter NaN
        trialinfo.pb_time(iT) = NaN;
    end
    
    if isnan(trialinfo.pb_time(iT)) % alsp populate other fields with NaN
        nNP = length(DEF_npWindowLabel);
        for iNPWin = 1:nNP
            nFreqs = length(DEF_freqBandLabel);
            for iFreq = 1:nFreqs
                fn = cat(2,DEF_npWindowLabel{iNPWin},'_',DEF_freqBandLabel{iFreq});
                trialinfo.(fn)(iT) = NaN;
            end
            fn = cat(2,DEF_npWindowLabel{iNPWin},'_spd');
            trialinfo.(fn)(iT) = NaN;
        end
    else
        cfg_np.t = trialinfo.pb_time(iT) + ExpKeys.TimeOnTrack(1); % matching MakeTrialinfo where this is subtracted;
        trl_np = ft_maketrl(cfg_np);
        
        cfg_temp = [];
        cfg_temp.trl = trl_np;
        data_np = ft_redefinetrial(cfg_temp,data);
        
        % TFR
        cfg_TFR.toi          = cfg_np.twin(1):0.01:cfg_np.twin(2);
        TFR = ft_freqanalysis(cfg_TFR, data_np);
        
        % go through cue windows and frequencies to process
        nNP = length(DEF_npWindowLabel);
        for iNPWin = 1:nNP
            
            cfg.temp = [];
            cfg_temp.t = DEF_npWindowTime{iNPWin};
            
            nFreqs = length(DEF_freqBandLabel);
            for iFreq = 1:nFreqs
                
                cfg_temp.f = DEF_freqBandFreq{iFreq};
                p = ft_getAvgTFRPow(cfg_temp,TFR);
                
                if(isnan(p))
                    pause;
                end
                
                % write output
                fn = cat(2,DEF_npWindowLabel{iNPWin},'_',DEF_freqBandLabel{iFreq});
                trialinfo.(fn)(iT) = p;
                
            end
            
            % speed
            t0 = cfg_np.t + cfg_temp.t(1);
            t1 = cfg_np.t + cfg_temp.t(2);
            fn = cat(2,DEF_npWindowLabel{iNPWin},'_spd');
            trialinfo.(fn)(iT) = nanmean(Data(Restrict(spdZ,t0,t1)));
            
        end
    end

end

%% compute full-session baseline
if cfg.doBaseline
    
    fprintf('\n** Computing baseline...\n\n');
    
    cfg_baseline = [];
    cfg_baseline.t = mean([ExpKeys.TimeOnTrack(1) ExpKeys.TimeOffTrack(2)]);
    cfg_baseline.twin = [ExpKeys.TimeOnTrack(1) ExpKeys.TimeOffTrack(2)]-cfg_baseline.t;
    cfg_baseline.mode = 'nlx';
    cfg_baseline.hdr = data.hdr;
    
    trl_baseline = ft_maketrl(cfg_baseline);
    
    cfg_temp = [];
    cfg_temp.trl = trl_baseline;
    data_baseline = ft_redefinetrial(cfg_temp,data);
    
    cfg_TFR              = []; % start with empty cfg
    cfg_TFR.output       = 'pow';
    cfg_TFR.keeptrials   = 'yes'; % need this for stats later
    cfg_TFR.channel      = {'vStr'};
    cfg_TFR.method       = 'mtmconvol';
    cfg_TFR.taper        = 'hanning';
    cfg_TFR.foi          = 25:1:110; % frequencies of interest
    cfg_TFR.t_ftimwin    = 19*(1./cfg_TFR.foi);
    cfg_TFR.toi          = cfg_baseline.twin(1):0.5:cfg_baseline.twin(2);
    
    TFR = ft_freqanalysis(cfg_TFR, data_baseline);
    
    for iFreq = 1:nFreqs
        cfg_temp.f = DEF_freqBandFreq{iFreq};
        cfg_temp.t = [-Inf Inf]
        p = ft_getAvgTFRPow(cfg_temp,TFR);
        
        if(isnan(p))
            error('NaN found.');
        elseif p == 0
            error('zero power!');
        end
        
        % write output
        fn = cat(2,'baseline_',DEF_freqBandLabel{iFreq});
        trialinfo.(fn)(1:iT) = p;
    end
end
%% get subject ID
[fp fd fe] = fileparts(pwd);
rat_ID_str = fd(1:4);
trialinfo.subject(1:iT) = strmatch(rat_ID_str,DEF_ratID);

%% block
%trialinfo.block = double(trialinfo.cue_time < ExpKeys.TimeOffTrack(1));
trialinfo.block = double(trialinfo.cue_time < (ExpKeys.TimeOffTrack(1)-ExpKeys.TimeOnTrack(1))); % correct for time normalization to ExpKeys.TimeOnTrack(1) in MakeTrialinfo()

%% reward amount on previous trial & time since
r_string = {'1 pellet dispensed','2 pellet dispensed','3 pellet dispensed','4 pellet dispensed','5 pellet dispensed'};

fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
EVTimeStamps = EVTimeStamps * 10^-6;

for iR = 1:length(r_string)
    fn = cat(2,'p',num2str(iR));
    all_t.(fn) = EVTimeStamps(strmatch(r_string{iR},EventStrings));
end

cfg_prev_rew = []; cfg_prev_rew.mode = 'prev'; cfg_prev_rew.fields = {'p1','p2','p3','p4','p5'};

for iT = nTrials:-1:1
   
    if iT == 1
       trialinfo.prevRew(iT) = NaN;
       trialinfo.tSinceRew(iT) = NaN;
       continue;
    end
    
   curr_t = trialinfo.cue_time(iT) + ExpKeys.TimeOnTrack(1); % matching MakeTrialinfo where this is subtracted;
   
   [prev_r_t,prev_r_string] = FindFieldTime(cfg_prev_rew,all_t,curr_t);
   
   trialinfo.prevRew(iT) = strmatch(prev_r_string,cfg_prev_rew.fields);
   trialinfo.tSinceRew(iT) = curr_t - prev_r_t;
    
end

%% save
if cfg.writeOutput
   [~,fd,~] = fileparts(pwd);
   fn_out = cat(2,fd,'_trialinfo.mat');
   save(fn_out,'trialinfo');
    
end