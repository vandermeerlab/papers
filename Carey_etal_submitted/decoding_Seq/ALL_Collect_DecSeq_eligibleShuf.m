%%
clear all;
todo = {'all','prerecord','taskrest','postrecord'}; % which task epochs to analyze
%todo = {'track-preCP'}; % on-track events

for iDo = 1:length(todo) % this is pretty inefficient -- should add inner within-session loop for task epochs
    
    cfg = [];
    cfg.prefix = 'S1_'; % prefix determining which decoding output files to load
    cfg.whichEvents = todo{iDo}; %{'prerecord','taskrest','taskrun','postrecord'}; % which events to process? can only select one
    cfg.whichSeq = 'all'; % {'all','fwd','bwd','either'}; % which sequences to process?
    cfg.sessions = {'food','water'};
    cfg.arms = {'left','right'};
    cfg.writeDiary = 1; % keep a text file record of command window history
    %cfg.output_fd = 'D:\projects\AlyssaTmaze\resultsFiles'; % where to place output files
    cfg.output_fd = 'C:\temp'; % where to place output files
    cfg.saveData = 1;
    cfg.rats = {'R042','R044','R050','R064'};
    cfg.cpbin = []; % if non-empty, restrict to sequences with start or end beyond this (relative to CP)
    cfg.minActiveCells = 4;
    cfg.minlen = 0.05; % otherwise, minimum length in s
    cfg.SWRoverlap = 1; % if 1, only keep events detected as SWR candidates
    cfg.SWRsuffix = ''; % filename suffix determining which SWR candidates to load
    cfg.RemoveOverlappingEvents = 1; % remove events that are classified as both L and R -- default is 1
    cfg.output_prefix = cat(2,cfg.prefix,cfg.SWRsuffix,'DecSeq_',cfg.whichEvents,'_',cfg.whichSeq,'_eligible_');
    cfg.processShuffles = 1;
    cfg.shuffle_include_p = 0.05; % discard detected interval if more than this fraction of its bins is detected in shuffles
    cfg.paperSessions = 1;
    
    if ~cfg.RemoveOverlappingEvents
        cfg.output_prefix = cat(2,cfg.output_prefix,'overlap_');
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
    
    disp(['Script name: ',mfilename])
    disp(' ')
    disp('You have selected: ')
    disp(cfg.whichEvents)
    disp(cfg.whichSeq)
    disp(cfg.rats)
    disp(' ')
    
    %% GET COMBINED DATA
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%                                                               %%%')
    disp('%%%               Combined sequences data (all rats):             %%%')
    disp('%%%                                                               %%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    %% init vars to place data into
    ALL_sig_seq = [];
    
    ALL_sig_seq.count = []; % number of sig events
    ALL_sig_seq.countN = []; % proportions of sig events
    ALL_sig_seq.arm = []; % left (1), right (2)
    ALL_sig_seq.type = []; % restriction type (food/water)
    ALL_sig_seq.sess = []; % session ID in case we want to restrict later
    ALL_sig_seq.firstChoice = []; % left (1), right (2)
    ALL_sig_seq.choice = []; % number of trials free choice
    ALL_sig_seq.choiceN = []; % proportions of trials
    ALL_sig_seq.allTrials = []; % number of trials (all, so including blocked)
    ALL_sig_seq.allTrialsN = []; % proportion of trials (all, so including blocked)
    
    ALL_sig_seqI.rsq = []; % track event-by-event properties
    ALL_sig_seqI.pval = [];
    ALL_sig_seqI.beta = [];
    ALL_sig_seqI.len = [];
    
    ALL_sig_seq.rat = {}; % keep track of which rat session data comes from, {'R042','R044','R050','R064'}
    ALL_sig_seqI.rat = {}; % as above, but for single event data
    
    swapfun = @(x) -x+3; % utility function to access variable from other condition
    
    
    %% collect
    nSessions = 0;
    for iFD = 1:length(fd)
        nSessions = nSessions +1;
        close all;
        
        %%% load data for this session %%%
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
        
        % track some useful things about this session
        this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
        switch metadata.taskvars.sequence{1}
            case 'L'
                firstChoice = 1;
            case 'R'
                firstChoice = 2;
        end
        seq_tosave = []; % for exporting as "candidates" for later use
        
        % obtain eligible sequences
        cfg.evt = evt; % include SWR candidates into selection
        this_seqR = getEligibleSequences(cfg,out); % KEY STEP: obtain eligible sequences from decoded output
        
        % restrict events to epoch of interest
        switch cfg.whichEvents
            case 'prerecord'
                this_seqR = restrict(this_seqR,0,ExpKeys.prerecord(2));
            case 'postrecord'
                this_seqR = restrict(this_seqR,ExpKeys.postrecord(1),Inf);
            case 'taskrest'
                this_seqR = restrict(this_seqR,metadata.taskvars.rest_iv);
            case 'taskrun'
                this_seqR = restrict(this_seqR,metadata.taskvars.trial_iv);
            case 'all'
                % do nothing
            case 'track-preCP' % on track, but before choice point
                pos = LoadPos([]);
                cfg_track = []; cfg_track.verbose = 0; cfg_track.startbuffer = -10;
                PreCPTrackIV = FindPreCPTrackIV(cfg_track,pos,metadata,ExpKeys);
                this_seqR = restrict(this_seqR,PreCPTrackIV);
            otherwise
                error('Unknown events %s!',cfg.whichEvents);
        end
        fprintf('%d sequences remaining after epoch selection.\n',length(this_seqR.tstart));
        
        % restrict to having at least some piece past a certain point on the track, if requested
        % should use idx output of SelectIV() to merge idxs of events that have start or end past requested point
        % could be the job of getEligibleSequences()
        %if ~isempty(cfg.cpbin)
        %    cfg_s = []; cfg_s.operation = '>'; cfg_s.threshold = 0.5;
        %    this_seqR = SelectIV(cfg_s,this_seqR,'loc');
        %end
        
        % select fwd/bwd if requested; could be the job of getEligbleSequences()
        if ~isempty(this_seqR.tstart)
            switch cfg.whichSeq
                
                case 'all'
                    keep_idx = 1:length(this_seqR.usr.pval);
                case 'fwd'
                    keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta > 0);
                case 'bwd'
                    keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta < 0);
                case 'either'
                    keep_idx = (this_seqR.usr.pval < 0.05);
                    
            end
            this_seqR = SelectIV([],this_seqR,keep_idx);
        end
        
        % count events for each arm
        nEventsTotal = length(this_seqR.tstart);
        for iS = 1:2 % loop over arms
            
            cfg_arm = []; cfg_arm.operation = '='; cfg_arm.threshold = iS;
            [this_arm_seqR,this_arm_idx] = SelectIV(cfg_arm,this_seqR,'side');
            
            nEvents = length(this_arm_idx);
            
            ALL_sig_seq.count = cat(1,ALL_sig_seq.count,nEvents);
            ALL_sig_seq.countN = cat(1,ALL_sig_seq.countN,nEvents./nEventsTotal);
            ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iS);
            ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
            ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
            ALL_sig_seq.firstChoice = cat(1,ALL_sig_seq.firstChoice,firstChoice);
            
            this_trials = length(strmatch(upper(cfg.arms{iS}(1)),metadata.taskvars.sequence));
            ALL_sig_seq.allTrials = cat(1,ALL_sig_seq.allTrials,this_trials);
            ALL_sig_seq.allTrialsN = cat(1,ALL_sig_seq.allTrialsN,this_trials./length(metadata.taskvars.sequence));
            
            choice_trial_idx = setdiff(1:length(metadata.taskvars.sequence),ExpKeys.forcedTrials);
            choice_trials = metadata.taskvars.sequence(choice_trial_idx);
            this_choice = length(strmatch(upper(cfg.arms{iS}(1)),choice_trials));
            ALL_sig_seq.choice = cat(1,ALL_sig_seq.choice,this_choice);
            ALL_sig_seq.choiceN = cat(1,ALL_sig_seq.choiceN,this_choice./length(choice_trials));
            
            ALL_sig_seq.rat = cat(1,ALL_sig_seq.rat,{sessionID(1:4)});
            
        end % of loop over arms
        
        % process candidates & save
        %evt = MergeIV([],seq_tosave);
        %S = LoadSpikes([]);
        %out_fn = cat(2,S.cfg.SessionID,'-DecSeqCand.mat');
        %save(out_fn,'evt');
        
    end
    
    %% now break out the data by rat
    data.all.ALL_sig_seq = ALL_sig_seq;
    %data.all.ALL_sig_seqI = ALL_sig_seqI;
    
    for iRat = 1:length(cfg.rats)
        
        this_rat = cfg.rats{iRat};
        
        % extract counts for each rat
        this_idx = strmatch(this_rat,data.all.ALL_sig_seq.rat);
        sfun = @(x) x(this_idx);
        
        data.(this_rat).ALL_sig_seq = structfun(sfun,data.all.ALL_sig_seq,'UniformOutput',false);
        
        % extract event stats for each rat
        %this_idx = strmatch(this_rat,data.all.ALL_sig_seqI.rat);
        %sfun = @(x) x(this_idx);
        
        %data.(this_rat).ALL_sig_seqI = structfun(sfun,data.all.ALL_sig_seqI,'UniformOutput',false);
        
    end
    
    %% compute some summary statistics
    toDo2 = cat(2,{'all'},cfg.rats);
    
    for iDo2 = 1:length(toDo2)
        
        this_do = toDo2{iDo2};
        
        food_left = data.(this_do).ALL_sig_seq.count(data.(this_do).ALL_sig_seq.arm == 1 & data.(this_do).ALL_sig_seq.type == 1);
        food_right = data.(this_do).ALL_sig_seq.count(data.(this_do).ALL_sig_seq.arm == 2 & data.(this_do).ALL_sig_seq.type == 1);
        food_leftN = data.(this_do).ALL_sig_seq.countN(data.(this_do).ALL_sig_seq.arm == 1 & data.(this_do).ALL_sig_seq.type == 1);
        food_rightN = data.(this_do).ALL_sig_seq.countN(data.(this_do).ALL_sig_seq.arm == 2 & data.(this_do).ALL_sig_seq.type == 1);
        
        water_left = data.(this_do).ALL_sig_seq.count(data.(this_do).ALL_sig_seq.arm == 1 & data.(this_do).ALL_sig_seq.type == 2);
        water_right = data.(this_do).ALL_sig_seq.count(data.(this_do).ALL_sig_seq.arm == 2 & data.(this_do).ALL_sig_seq.type == 2);
        water_leftN = data.(this_do).ALL_sig_seq.countN(data.(this_do).ALL_sig_seq.arm == 1 & data.(this_do).ALL_sig_seq.type == 2);
        water_rightN = data.(this_do).ALL_sig_seq.countN(data.(this_do).ALL_sig_seq.arm == 2 & data.(this_do).ALL_sig_seq.type == 2);
        
        d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];
        
        data.(this_do).food_left = nansum(food_left); data.(this_do).food_leftN = nanmean(food_leftN);
        data.(this_do).food_right = nansum(food_right); data.(this_do).food_rightN = nanmean(food_rightN);
        data.(this_do).water_left = nansum(water_left); data.(this_do).water_leftN = nanmean(water_leftN);
        data.(this_do).water_right = nansum(water_right); data.(this_do).water_rightN = nanmean(water_rightN);
        
        disp(' ')
        disp([this_do, ', number of significant sequences: '])
        disp(['Left, food restr, ',num2str(data.(this_do).food_left)])
        disp(['Right, food restr ',num2str(data.(this_do).food_right)])
        disp(['Left, water restr ',num2str(data.(this_do).water_left)])
        disp(['Right, water restr ',num2str(data.(this_do).water_right)])
        disp(' ')
        
        % p_foodseq = tmaze_coin_pval(max([d(1) d(2)]),d(1)+d(2));
        % p_waterseq = tmaze_coin_pval(max([d(3) d(4)]),d(3)+d(4));
        % disp(' ');
        % disp('***************************')
        % disp('** binomial test results **')
        % disp('***************************')
        % disp(['Food day pval: ',num2str(p_foodseq)])
        % disp(['Water day pval: ',num2str(p_waterseq)])
        
        temp = data.(this_do); % because I don't want to type as much down there
        
        % [nleft on food day, nright on food day, nleft on water day, nright on water day]
        obs = [temp.food_left temp.food_right temp.water_left temp.water_right]; % the observed, real data
        
        nFood = obs(1) + obs(2); % total number of food day trials. could also do temp.nFood
        nWater = obs(3) + obs(4); % total number of water day trials
        ratioL = (obs(1)+obs(3)) / sum(obs); % ratio of L trials compared to total trials
        ratioR = (obs(2)+obs(4)) / sum(obs); % ratio of R trials
        
        exp = [nFood*ratioL nFood*ratioR nWater*ratioL nWater*ratioR]; % the "expected" values
        
        bins = 0:3;
        edges = -0.5:3.5;
        disp(' ');
        disp('***************************************')
        disp('**                                   **')
        disp('**      chi square test results      **')
        disp('**                                   **')
        disp('***************************************')
        [h,p,stats] = chi2gof(bins,'Alpha',0.01,'Edges',edges,'freq',obs,'expected',exp,'Emin',1) % no semicolon b/c want to see output
        
    end % of local toDo2 loop over rats
    
    %%
    data.date = datestr(now);
    data.script = mfilename;
    data.cfg = cfg;
    
    if cfg.writeDiary; diary off; end
    
    if cfg.saveData; cd(cfg.output_fd); save([cfg.output_prefix,'out.mat'],'data'); end
    
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~                      End of script run                          ~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
end % of mega toDo loop over epochs
cd(iWasHere)
warning on