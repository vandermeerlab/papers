%% MASTER_ClassifyGammaEvents.m
% classifies previously gamma events based on spiking data
%
% Julien Catanese & Matthijs van der Meer

%% set paths
restoredefaultpath;
cd('D:\My_Documents\GitHub\fieldtrip');
ft_defaults;

rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal\');

addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\tasks\Julien_linear_track')); % Detect events, CountCycles live here

%% load gamma events (use MASTER_CollectGammaEvents_adrlab.m to obtain) -- puts ALL_evt variable in workspace
cd('D:\My_Documents\Dropbox\projects\Julien_multiLFP\2016-01-10');
%load(FindFile('gamma*.mat'));

%% define what to run
rats = {'R117','R119','R131','R132'};

PARAM_twin = 0.2; % half-width of time window
PARAM_control_dt = 5; % offset for control classifier
PARAM_nIter = 1000; % how many iterations to run for each cell count (random subsets)
PARAM_opts = statset('UseParallel','always','UseSubstreams','never'); % classifier options'
PARAM_minEventsWithSpikes = 10; % only include neurons with this 
debug = 0;

%% set up classifiers
linearDiscriminant = @(xtr,ytr,xte,yte) classify(xte,xtr,ytr,'linear');
knn9 = @(xtr,ytr,xte,yte) knnclassify(xte,xtr,ytr,9);
naiveBayes = @(xtr,ytr,xte,yte) classify(xte,xtr,ytr,'diagLinear');

classifier_list =  {'linearDiscriminant','naiveBayes','knn9'}; 
%classifier_list =  {'linearDiscriminant'}; 
%%
available_rats = fieldnames(ALL_evt);
nRats = 0; nSessions = 0;

for iRat = 1:length(rats)
    
    this_rat = rats{iRat};
    
    if ~strmatch(this_rat,available_rats)
       warning('Rat %s not available -- skipping...',rats{iRat});
       continue;
    end
    
    available_sessions = fieldnames(ALL_evt.(this_rat));
    
    for iSession = 1:length(available_sessions)
    
        this_session = available_sessions{iSession};
        this_session_data = ALL_evt.(this_rat).(this_session);
        
        fprintf('\n\nProcessing session %s...\n',this_session);
        
        this_fd = this_session_data.fd;
        cd(this_fd);
        
        % load the data
        LoadExpKeys;
        
        please = []; please.fc = ExpKeys.goodGamma_vStr(1);
        csc = LoadCSC(please);
        
        S = LoadSpikes([]);
        
        % remove crossed events
        this_lg = DifferenceIV([],this_session_data.lg,this_session_data.hg);
        this_hg = DifferenceIV([],this_session_data.hg,this_session_data.lg);
        
        % get times at feeder
        at_feeder = Julien_atFeederTimes();
        
        % only keep events where rat is on pedestal or at feeder
        keep_iv = UnionIV([],at_feeder,iv(-Inf,ExpKeys.TimeOnTrack-5));
        keep_iv = UnionIV([],keep_iv,iv(ExpKeys.TimeOffTrack+5,Inf));
        
        orig_nlg = length(this_lg.tstart); this_lg = restrict(this_lg,keep_iv);
        fprintf('\n MASTER_ClassifyGammaEvents.m: %d/%d lg events kept.\n',length(this_lg.tstart),orig_nlg);
        orig_hlg = length(this_hg.tstart);  this_hg = restrict(this_hg,keep_iv);
        fprintf('\n MASTER_ClassifyGammaEvents.m: %d/%d hg events kept.\n',length(this_hg.tstart),orig_hlg);
               
        % create fixed-size events
        lg_t = IVcenters(this_lg); hg_t = IVcenters(this_hg);
        this_lg = iv(lg_t-PARAM_twin,lg_t+PARAM_twin); %this_lg_control = iv(lg_t-PARAM_twin+PARAM_control_dt,lg_t+PARAM_twin+PARAM_control_dt);
        this_hg = iv(hg_t-PARAM_twin,hg_t+PARAM_twin); %this_hg_control = iv(hg_t-PARAM_twin+PARAM_control_dt,hg_t+PARAM_twin+PARAM_control_dt);
        
        % find matched events, as a control
        cfg_match = []; cfg_match.twin = [-PARAM_control_dt PARAM_control_dt]; cfg_match.evt_twin = [-PARAM_twin PARAM_twin];
        ctrl_evt = MatchGammaEvents(cfg_match,csc,lg_t,hg_t);
        this_lg_control = iv(ctrl_evt.lg-PARAM_twin,ctrl_evt.lg+PARAM_twin);
        this_hg_control = iv(ctrl_evt.hg-PARAM_twin,ctrl_evt.hg+PARAM_twin);
        
        if debug
            cfg_plot = []; cfg_plot.display = 'tsd';
            PlotTSDfromIV(cfg_plot,this_lg,csc);
            cfg_plot.iv_only = 1; cfg_plot.fgcol = 'g';
            PlotTSDfromIV(cfg_plot,this_lg_control,csc);
        end
        
        % create spike count vector for events
        cfg_sc = [];     
        cfg_sc.iv = this_lg; sc_lg = getSpikeCount(cfg_sc,S);
        cfg_sc.iv = this_hg; sc_hg = getSpikeCount(cfg_sc,S);
        sc_merged = cat(2,sc_lg,sc_hg);
        
        cfg_sc.iv = this_lg_control; sc_lg_control = getSpikeCount(cfg_sc,S);
        cfg_sc.iv = this_hg_control; sc_hg_control = getSpikeCount(cfg_sc,S);
        sc_merged_control = cat(2,sc_lg_control,sc_hg_control);
        
        % remove cells with same number of spikes in all events (trips up some classifiers)
        nCells = size(sc_merged,1);
        keep_idx = ones(nCells,1);
        for iC = nCells:-1:1
            this_sc = sc_merged(iC,:); this_scC = sc_merged_control(iC,:);
            if length(unique(this_sc)) == 1 | length(unique(this_scC)) == 1
                keep_idx(iC) = 0;
            end
        end
        sc_merged = sc_merged(logical(keep_idx),:);
        sc_merged_control = sc_merged_control(logical(keep_idx),:);
        nCells = size(sc_merged,1);
        
        % remove cells with fewer than 10 spikes across all events (breaks
        % classifiers, esp when crossvalidating)
        keep_idx = sum(sc_merged > 0,2) > PARAM_minEventsWithSpikes & sum(sc_merged_control > 0,2) > PARAM_minEventsWithSpikes;
        sc_merged = sc_merged(keep_idx,:);
        sc_merged_control = sc_merged_control(keep_idx,:);
        nCells = size(sc_merged,1);
        
        % create labels (0: lg, 1: hg)
        nlg = length(this_lg.tstart); nhg = length(this_hg.tstart);
        label = cat(1,zeros(nlg,1),ones(nhg,1));
        
        % classify
        nCount = 1:5:nCells;
        sc_merged = sc_merged'; sc_merged_control = sc_merged_control';
 
        vals = nan(length(classifier_list),length(nCount),PARAM_nIter);
        valsC = nan(length(classifier_list),length(nCount),PARAM_nIter);
        valsS = nan(length(classifier_list),length(nCount),PARAM_nIter);
        for iClassifier = 1:length(classifier_list)
            for iCount = 1:length(nCount)
                fprintf('\nClassifier %d/%d, count %d/%d...\n',iClassifier,length(classifier_list),iCount,length(nCount));
                for iI = PARAM_nIter:-1:1 % iterate over number of neurons included
                    
                    % select subset of neurons
                    neuron_idx = randperm(nCells); neuron_idx = neuron_idx(1:nCount(iCount));
                    
                    % select subset of data points from class with most
                    % data
                    if nlg > nhg
                        keep_idx = randperm(nlg);
                        keep_idx = keep_idx(1:nhg);
                        keep_idx = [keep_idx nlg+1:nlg+nhg];
                        this_sc_merged = sc_merged(keep_idx,:);
                        this_sc_merged_control = sc_merged_control(keep_idx,:);
                        this_label = label(keep_idx);
                        this_nlg = nhg; this_nhg = nhg;
                    elseif nhg > nlg
                        keep_idx = randperm(nhg);
                        keep_idx = keep_idx(1:nlg);
                        keep_idx = [1:nlg nlg+keep_idx];
                        this_sc_merged = sc_merged(keep_idx,:);
                        this_sc_merged_control = sc_merged_control(keep_idx,:);
                        this_label = label(keep_idx);
                        this_nhg = nlg; this_nlg = nlg;
                    else
                        this_sc_merged = sc_merged;
                        this_sc_merged_control = sc_merged_control;
                        this_label = label;
                        this_nhg = nhg; this_nlg = nlg;
                    end
                    
                    try 
                        vals(iClassifier,iCount,iI) = crossval('mcr',this_sc_merged(:,neuron_idx),this_label,'Predfun',eval(classifier_list{iClassifier}),'options',PARAM_opts); 
                    catch
                        disp(sprintf('WARNING: classifier %s failed for neuron count %d',classifier_list{iClassifier},nCount(iCount)));
                        %pause;
                    end
                    
                    try
                        valsC(iClassifier,iCount,iI) = crossval('mcr',this_sc_merged_control(:,neuron_idx),this_label,'Predfun',eval(classifier_list{iClassifier}),'options',PARAM_opts);
                    catch
                        disp(sprintf('WARNING: classifier %s failed for neuron count %d (control)',classifier_list{iClassifier},nCount(iCount)));
                        %pause;
                    end
                    
                    % for within-label shuffle, reassign spike counts to
                    % random events for each neuron independently (breaks
                    % between neuron correlations)
                    this_sc_shuf = this_sc_merged;
                    for iC = 1:nCells
                       this_sc_shuf(1:this_nlg,iC) = this_sc_shuf(randperm(this_nlg),iC); 
                       this_sc_shuf(this_nlg+1:end,iC) = this_sc_shuf(this_nlg+randperm(this_nhg),iC); 
                    end
                    
                    try
                        valsS(iClassifier,iCount,iI) = crossval('mcr',this_sc_shuf(:,neuron_idx),this_label,'Predfun',eval(classifier_list{iClassifier}),'options',PARAM_opts);
                    catch
                        disp(sprintf('WARNING: classifier %s failed for neuron count %d (shuf)',classifier_list{iClassifier},nCount(iCount)));
                    end
                    
                end
            end
        end % of loop over classifier
        
        % average data and add to aggregator
        ALL_evt.(this_rat).(this_session).classify.nCells = nCount;
        ALL_evt.(this_rat).(this_session).classify.clist = classifier_list;
        ALL_evt.(this_rat).(this_session).classify.cdata = sq(nanmean(vals,3));
        ALL_evt.(this_rat).(this_session).classify.cdataC = sq(nanmean(valsC,3));
        ALL_evt.(this_rat).(this_session).classify.cdataS = sq(nanmean(valsS,3));
        ALL_evt.(this_rat).(this_session).classify.cdata_min = sq(min(vals,[],3));
        ALL_evt.(this_rat).(this_session).classify.cdataC_min = sq(min(valsC,[],3));
        ALL_evt.(this_rat).(this_session).classify.cdataS_min = sq(min(valsS,[],3));
        ALL_evt.(this_rat).(this_session).classify.bestSingleNeuron = min(vals(:,1,:),[],3); % easier to handle later, should be redundant
        ALL_evt.(this_rat).(this_session).classify.classifier_list = classifier_list;
        ALL_evt.(this_rat).(this_session).classify.ctrl_evt = ctrl_evt; 
        
    end % sessions
    
end % rats

