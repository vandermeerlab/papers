%% ALL_BehavRegress.m

% fit various regression models relating sequence SWR content to behavioral choice

% requires ALL_Generate_DecSeq be run; use cfg.prefix setting below to
% specify which output files are to be loaded

%% list of potential regressors collected by this script


%%
clear all;
cfg = [];
cfg.prefix = 'S1_'; % prefix determining which decoding output files to load
cfg.whichSeq = 'all'; % {'all','fwd','bwd','either'}; % which sequences to process?
cfg.arms = {'left','right'};
cfg.sessions = {'food','water'};
cfg.writeDiary = 1; % keep a text file record of command window history
cfg.output_fd = 'C:\temp'; % where to place output files
cfg.saveData = 1;
cfg.cpbin = []; % if non-empty, restrict to sequences with start or end beyond this (relative to CP)
cfg.minActiveCells = 4;
cfg.minlen = 0.05; % otherwise, minimum length in s
cfg.SWRoverlap = 1; % if 1, only keep events detected as SWR candidates
cfg.SWRsuffix = ''; % filename suffix determining which SWR candidates to load
cfg.RemoveOverlappingEvents = 1; % remove events that are classified as both L and R -- default is 1
cfg.output_prefix = cat(2,cfg.prefix,cfg.SWRsuffix,'DecSeq_',cfg.whichSeq,'_BehavRegress');
cfg.processShuffles = 1;
cfg.shuffle_include_p = 0.05; % discard detected interval if more than this fraction of its bins is detected in shuffles

if ~cfg.RemoveOverlappingEvents
    cfg.output_prefix = cat(2,cfg.output_prefix,'_overlap');
end

%%
fd = getTmazeDataPath(cfg);
iWasHere = pwd;

if cfg.writeDiary % save command window text
    warning off
    cd(cfg.output_fd)
    diary([cfg.output_prefix,'stats','.txt'])
    cd(cfg.output_fd)
    disp(' ')
    disp(date)
    disp(' ')
end
fd([1 7 8 9 12]) = [];

%% init vars to place data into
ALL_sig_seq = [];
swapfun = @(x) -x+3; % utility function to access variable from other condition

%% collect
nSessions = 0;
ALL_t = table;
for iFD = 1:length(fd)
    nSessions = nSessions + 1;
    close all;
    
    %%% STEP 1: load data for this session %%%
    cd(fd{iFD});
    LoadExpKeys;
    cfg_cand = []; cfg_cand.suffix = cfg.SWRsuffix; LoadCandidates(cfg_cand);
    LoadMetadata;
    
    cd([fd{iFD},'\files'])
    [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
    this_file = FindFiles([cfg.prefix,sessionID,'-DecSeq_data.mat']);
    
    if isempty(this_file)
        fprintf('Session %s: no DecSeq file found, skipping...\n',fd{iFD});
        continue;
    end
    
    load(this_file{1}); cd .. % loads 'out' variable saved by ALL_Generate_DecSeq.m
    this_t = table;
    
    %%% STEP 3: obtain all eligible sequences for analysis %%%
    cfg_seq = cfg;
    cfg_seq.evt = evt; % pass in candidate events
    this_seq = getEligibleSequences(cfg_seq,out); % KEY STEP: obtain eligible events from Generate_DecSeq output
    
    %%% STEP 4: populate table of regressors and dependent variables %%%
    clear trial_nseq_pre trial_nseq_post trial_pleft_pre trial_pleft_post
    clear trial_modal_pre trial_modal_post trial_last_pre trial_first_post
    clear trial_no trial_side trial_forced trial_type
    clear trial_Spleft trial_Sleft_taken trial_Sleft_chosen trial_Sped_side trial_SratID trial_Sdep

    nTrials = length(metadata.taskvars.sequence);
    for iTrial = 1:nTrials
       
        % session-level sequence proportions
        trial_Spleft(iTrial) = sum(this_seq.usr.side == 1)./length(this_seq.tend);
        
        % get preceding and following sequences for this trial
        if iTrial == 1 % preceding segment is prerecord
            
            pre_seq = restrict(this_seq,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
            
            cfg_ivf = []; cfg_ivf.mode = 1;
            post_idx = FindNearestIV(cfg_ivf,metadata.taskvars.rest_iv,metadata.taskvars.trial_iv.tstart(iTrial));
            if ~isempty(post_idx)
                post_seq = restrict(this_seq,metadata.taskvars.rest_iv.tstart(post_idx),metadata.taskvars.rest_iv.tend(post_idx));
            else
                post_seq = iv;
            end
            
            trial_type(iTrial) = 1;
            
        elseif iTrial == nTrials % following segment includes postrecord
                
            cfg_ivf = []; cfg_ivf.mode = -1;
            pre_idx = FindNearestIV(cfg_ivf,metadata.taskvars.rest_iv,metadata.taskvars.trial_iv.tstart(iTrial));
            if ~isempty(pre_idx)
                pre_seq = restrict(this_seq,metadata.taskvars.rest_iv.tstart(iTrial-1),metadata.taskvars.rest_iv.tend(iTrial-1));
            else
                pre_seq = iv;
            end
                        
            cfg_ivf = []; cfg_ivf.mode = 1;
            post_idx = FindNearestIV(cfg_ivf,metadata.taskvars.rest_iv,metadata.taskvars.trial_iv.tstart(iTrial));
            if ~isempty(post_idx)
                post_seq = restrict(this_seq,metadata.taskvars.rest_iv.tstart(post_idx),metadata.taskvars.rest_iv.tend(post_idx));
            else
                post_seq = iv;
            end
            
            post_seq = UnionIV([],post_seq,restrict(this_seq,ExpKeys.postrecord(1),ExpKeys.postrecord(2)));
            
            trial_type(iTrial) = 2;
            
        else % can use rest_iv
                
            cfg_ivf = []; cfg_ivf.mode = -1;
            pre_idx = FindNearestIV(cfg_ivf,metadata.taskvars.rest_iv,metadata.taskvars.trial_iv.tstart(iTrial));
            if ~isempty(pre_idx)
                pre_seq = restrict(this_seq,metadata.taskvars.rest_iv.tstart(iTrial-1),metadata.taskvars.rest_iv.tend(iTrial-1));
            else
                pre_seq = iv;
            end
                        
            cfg_ivf = []; cfg_ivf.mode = 1;
            post_idx = FindNearestIV(cfg_ivf,metadata.taskvars.rest_iv,metadata.taskvars.trial_iv.tstart(iTrial));
            if ~isempty(post_idx)
                post_seq = restrict(this_seq,metadata.taskvars.rest_iv.tstart(post_idx),metadata.taskvars.rest_iv.tend(post_idx));
            else
                post_seq = iv;
            end
            trial_type(iTrial) = 0;
            
        end
        
        % add variables to table
        % number of sequences, proportion left, modal direction left/right,
        % last sequence, proportion fwd/bwd
        trial_nseq_pre(iTrial) = length(pre_seq.tstart);
        trial_nseq_post(iTrial) = length(post_seq.tstart);
        trial_pleft_pre(iTrial) = sum(pre_seq.usr.side == 1)./length(pre_seq.tend);
        trial_pleft_post(iTrial) = sum(post_seq.usr.side == 1)./length(post_seq.tend);
        trial_modal_pre(iTrial) = mode(pre_seq.usr.side);
        trial_modal_post(iTrial) = mode(post_seq.usr.side);
      
        if length(pre_seq.tstart) > 0
            trial_last_pre(iTrial) = pre_seq.usr.side(end);
        else
            trial_last_pre(iTrial) = NaN;
        end
        
        if length(post_seq.tstart) > 0
            trial_first_post(iTrial) = post_seq.usr.side(1);
        else
            trial_first_post(iTrial) = NaN;
        end
        
        trial_no(iTrial) = iTrial; 
        
        switch metadata.taskvars.sequence{iTrial}
            case 'L'
                trial_side(iTrial) = 1;
            case 'R'
                trial_side(iTrial) = 2;
        end
        
        trial_forced(iTrial) = ismember(iTrial,ExpKeys.forcedTrials);
        
        % add the session level stuff, like left-chosen and left-taken,
        % pedestal side, session number, rat ID, food/water dep, etc.
        
        trl_idx = 1:nTrials;
        trl_seq = metadata.taskvars.sequence(trl_idx);
        trial_Sleft_taken(iTrial) = length(strmatch('L',trl_seq))./length(trl_seq); % proportion of lefts chosen for session (non forced)
        
        trl_idx = setxor(trl_idx,ExpKeys.forcedTrials); % exclude forced trials
        trl_seq = metadata.taskvars.sequence(trl_idx);
        
        trial_Sleft_chosen(iTrial) = length(strmatch('L',trl_seq))./length(trl_seq); % proportion of lefts taken for session
        
        switch ExpKeys.Pedestal
            case 'L'
                trial_Sped_side(iTrial) = 1;
            case 'R'
               trial_Sped_side(iTrial) = 2; % 1 = left, 2 = right
        end
        
        switch ExpKeys.goodSWR{1}(1:4)
            case 'R042'
                trial_SratID(iTrial) = 1;
            case 'R044'
                trial_SratID(iTrial) = 2;
            case 'R050'
                trial_SratID(iTrial) = 3;
            case 'R064'
                trial_SratID(iTrial) = 4;
        end
        
        switch ExpKeys.RestrictionType
            case 'food'
                trial_Sdep(iTrial) =  1;% 1 = food, 2 = water
            case 'water'
                trial_Sdep(iTrial) =  2;
        end
        
    end % of trial loop
    
    % add all trial data for this session to master table
    vn = whos('trial_*');
    vn = {vn.name};
    
    for ivn = 1:length(vn)
        eval(sprintf('this_t.%s = %s'';',vn{ivn},vn{ivn}));
    end
        
    ALL_t = [ALL_t; this_t];
    
end % of sessions loop

%% ANALYSIS 1: predict behavior on free choice trials
ALL_t.trial_SratID = categorical(ALL_t.trial_SratID); % rat ID (1-4)
ALL_t.trial_side = categorical(ALL_t.trial_side); % behavioral choice (left/right, 1/2)
ALL_t.trial_Sdep = categorical(ALL_t.trial_Sdep); % motivational state (food/water restricted, 1/2)
ALL_t.trial_modal_pre = categorical(ALL_t.trial_modal_pre); % modal sequence content prior to trial (left/right arm, 1/2)
ALL_t.trial_last_pre = categorical(ALL_t.trial_last_pre); % content of last sequence prior to trial (left/right arm, 1/2)

% trial selection
%ALL_t.trial_side = double(ALL_t.trial_side)-1; % check if categorical vs. numeric response variable makes a difference (doesn't)
rows = ALL_t.trial_type == 0 & ALL_t.trial_forced == 0; % exclude first and last trials for now
ALL_t2 = ALL_t(rows,:);

% fit baseline model
modelspec = 'trial_side ~ 1 + trial_Sdep'; % motivational state only
%modelspec = 'trial_side ~ 1'; % sanity check, intercept only
glm = fitglm(ALL_t2,modelspec,'Link','logit','Distribution','binomial')
hist(glm.Residuals.Raw,100) % plot residuals -- note that model predictions aren't 0/1 for some reason

out = glm.Fitted;
out = categorical(out.Response > 0.5); data = categorical(double(ALL_t2.trial_side)-1 > 0.5);
sum(out == data) % count number of correct predictions

% test if replay content helps
modelspec = 'trial_side ~ 1 + trial_Sdep + trial_Spleft'; % add session-wide replay content
glm2 = fitglm(ALL_t2,modelspec,'Link','logit','Distribution','binomial')
hist(glm2.Residuals.Raw,100)

out = glm2.Fitted;
out = categorical(out.Response > 0.5); data = categorical(double(ALL_t2.trial_side)-1 > 0.5);
sum(out == data) % note, same number of correct predictions

% report some model comparison info
glm.ModelCriterion
glm2.ModelCriterion
p = 1-chi2cdf(glm.Deviance-glm2.Deviance,1) % note model 2 fits better

%% improve with trial by trial variables? no.
modelspec = 'trial_side ~ 1 + trial_Sdep + trial_Spleft + trial_modal_pre'; % add session-wide replay content
glm3 = fitglm(ALL_t2,modelspec,'Link','logit','Distribution','binomial')