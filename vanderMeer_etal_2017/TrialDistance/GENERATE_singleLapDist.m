%% preamble
clear all; pack

%MASTER_root = 'D:\My_Documents\Dropbox\projects\Alyssa'; % replace this with home folder of your project
MASTER_root = 'C:\Users\mvdm\Dropbox\projects\Alyssa'; % replace this with home folder of your project
cd(MASTER_root);
MASTER_path; % reset and then set up path for this project

%% master config
cfg_master = [];
cfg_master.dt = 0.025;
cfg_master.TCsmooth = 0.5; % bins SD, choose one
cfg_master.QsmoothSD = 0.005; % time SD, choose one
cfg_master.minSpikes = 25; % remove cells with less than this number of spikes during run
cfg_master.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
cfg_master.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
cfg_master.tc_baseline = 0.1; % baseline firing rate, replaces zeros in TC when unsmoothed with smoothed Q
cfg_master.load_questionable_cells = 1;
cfg_master.includeAllCells = 1; % otherwise, only place cells
cfg_master.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
cfg_master.trackExcludeEnd = 15;
cfg_master.nSpikesHist = -0.5:105;
cfg_master.matchTrials = 0;

cfg_master.nEncLaps = -10:10; % lap distances to find

%% find data folders
fd_cfg = []; fd_cfg.requireCandidates = 0; % don't need previously saved sharp wave-ripple (SWR) candidates here
fd_cfg.requireEvents = 1;
fd = getTmazeDataPath(fd_cfg);

%% init some vars
clear expCond;

expCond(1).label = 'left';
expCond(2).label = 'right';
    
for iCond = 1:length(expCond)
    ALL_decErr.(expCond(iCond).label).meanErr = nan(length(fd),length(cfg_master.nEncLaps));
    ALL_decErr.(expCond(iCond).label).spaceErr = nan(length(fd),length(cfg_master.nEncLaps),cfg_master.nBins-1);
end

%%
for iFD = 1:length(fd)
    
    cd(fd{iFD});
    fprintf('Entering session %d/%d...\n',iFD,length(fd));
    
    %% load data
    please = []; please.load_questionable_cells = cfg_master.load_questionable_cells;
    S = LoadSpikes(please);
   
    LoadExpKeys;
    LoadMetadata;
    
    cfg_pos = []; cfg_pos.convFact = ExpKeys.convFact;
    pos = LoadPos(cfg_pos); % pos is now in cm
        
    %% set up data structs for L, R
    clear expCond;
    expCond(1).label = 'left';
    expCond(2).label = 'right';
    
    % match trials
    if cfg_master.matchTrials
        [left,right,left_indices,right_indices] = GetMatchedTrials([],metadata,ExpKeys);
        expCond(1).t = left; 
        expCond(1).t_idx = left_indices; % indices into combined lap sequence (both L, R)
        expCond(2).t = right; 
        expCond(2).t_idx = right_indices; % indices into combined lap sequence (both L, R)
    else
        metadata = RemoveBadTrials([],metadata,ExpKeys);
        
        left.tstart = metadata.taskvars.trial_iv.tstart(strcmp('L',metadata.taskvars.sequence));
        left.tend = metadata.taskvars.trial_iv.tend(strcmp('L',metadata.taskvars.sequence));
        left = iv(left.tstart,left.tend);
        
        right.tstart = metadata.taskvars.trial_iv.tstart(strcmp('R',metadata.taskvars.sequence));
        right.tend = metadata.taskvars.trial_iv.tend(strcmp('R',metadata.taskvars.sequence));
        right = iv(right.tstart,right.tend);
        
        expCond(1).t = left; % indices into combined lap sequence (both L, R)
        expCond(1).t_idx = find(strcmp('L',metadata.taskvars.sequence));
        expCond(2).t = right; % indices into combined lap sequence (both L, R)
        expCond(2).t_idx = find(strcmp('R',metadata.taskvars.sequence));
    end
    
    nLapsMax = max(length(expCond(1).t.tstart),length(expCond(2).t.tstart));
    
    % tighter run boundaries: start run based on center photobeam break
    evt = getEvents_Tmaze();
    for iCond = 1:length(expCond)
       
        this_runStart = expCond(iCond).t.tstart;
        for iT = 1:length(this_runStart)
            next_centerpb_break_idx = nearest_idx3(this_runStart,evt.center_pb,1); 
            expCond(iCond).t.tstart = evt.center_pb(next_centerpb_break_idx)'+cfg_master.trialStartOffset;
        end
    end
    
    expCond(1).coord = metadata.coord.coordL_cm;
    expCond(2).coord = metadata.coord.coordR_cm;
    
    expCond(1).S = S;
    expCond(2).S = S;
    
    %% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
    fprintf('Linearizing...');
    
    nCond = length(expCond);
    for iCond = 1:nCond
        
        cfg_linpos = []; cfg_linpos.Coord = expCond(iCond).coord;
        expCond(iCond).linpos = LinearizePos(cfg_linpos,pos);
        
        % ensure that linpos is now in cm
        expCond(iCond).linpos.data = (expCond(iCond).linpos.data ./ length(cfg_linpos.Coord)).*ExpKeys.pathlength;
        
    end
    
    %% find intervals where rat is running
    spd = getLinSpd([],pos); % get speed (in "camera pixels per second")
    
    cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 5;
    run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 5 pix/s
    
    %% exclude positions at beginning and end of track
    cfg_track1 = []; cfg_track1.method = 'raw'; cfg_track1.threshold = cfg_master.trackExcludeStart;
    cfg_track2 = []; cfg_track2.method = 'raw'; cfg_track2.dcn = '<'; cfg_track2.threshold = ExpKeys.pathlength - cfg_master.trackExcludeEnd;
    
    for iCond = 1:nCond
        track_iv1 = TSDtoIV(cfg_track1,expCond(iCond).linpos);
        
        fh = @(x) restrict(x,track_iv1);
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        
        track_iv2 = TSDtoIV(cfg_track2,expCond(iCond).linpos);
        
        fh = @(x) restrict(x,track_iv2);
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        
    end
    
    %% deal with some rat specific oddities
    [fp fn fe] = fileparts(fd{iFD});
    switch fn(1:4)
        case 'R042'
            % deal with non-spike sorted chewing intervals
            tf = FindFile('*times.mat');
            if isempty(tf)
               error('No *times.mat file found for R042'); 
            else
                load(tf);
                chew_iv = iv(t_start*10^-6,t_end*10^-6); % times are in neuralynx timestamps, so convert to s
            end
            
        case 'R044'
            % deal with HS detachments -- remove extended times with no spikes
            cfg_Q = []; cfg_Q.dt = 1;
            Q = MakeQfromS(cfg_Q,S);
            spk_count = tsd(Q.tvec,sum(Q.data));
            
            cfg_det = []; cfg_det.threshold = 0.5; cfg_det.dcn = '>'; cfg_det.method = 'raw';
            chew_iv = TSDtoIV(cfg_det,spk_count);
            
        otherwise
            chew_iv = [];
            
    end
    
    %% restrict (linearized) position data and spike data to desired intervals
    fprintf('Restricting data...');
    for iCond = 1:nCond
        
        fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        
        fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        
        if ~isempty(chew_iv)
            fh = @(x) restrict(x,chew_iv); % restrict S and linpos to non-detached times
            expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        end
        
        % also remove cells with insufficient spikes
        expCond(iCond).S = removeEmptyCells(expCond(iCond).S);
        
        spk_count = getSpikeCount([],expCond(iCond).S);
        cell_keep_idx = spk_count >= cfg_master.minSpikes;
        
        expCond(iCond).S = SelectTS([],expCond(iCond).S,cell_keep_idx);
    end
    
    %% main loop
    for iCond = 1:nCond
        
        nLaps = length(expCond(iCond).t.tstart);
        
        % first, initialize P and trueZ
        for iDist = 1:length(cfg_master.nEncLaps)
            
            this_Psplice{iDist} = tsd;
            this_Psplice{iDist}.usr.nActiveNeurons = []; % need to initialize this for UnionTSD() to work
            this_Psplice{iDist}.usr.nActiveNeuronsPassed = [];
            
            this_trueZsplice{iDist} = tsd;
            
        end
        
        %%
        for iLap = 1:nLaps % lap to decode
            
            fprintf('.lap %d/%d...\n',iLap,nLaps);
                        
            %%
            for iDist = 1:length(cfg_master.nEncLaps) % laps to base encoding model on
                
                this_lap = expCond(iCond).t_idx(iLap);
                this_dist = cfg_master.nEncLaps(iDist);
                fprintf('..Dist %d/%d...\n',iDist,length(cfg_master.nEncLaps));
                
                %% construct tuning curves
                cfg_tc = [];
                cfg_tc.binEdges{1} = linspace(0,ExpKeys.pathlength,cfg_master.nBins); % ~3cm bins for full maze (smaller for R042...)
                
                if cfg_master.TCsmooth ~= 0
                    cfg_tc.smoothingKernel = gausskernel(51,cfg_master.TCsmooth);
                end
                
                % restrict to the encoding model for this lap
                %nAvailableLaps = length(expCond(iCond).t.tstart);
                
                % goal: find set of laps that doesn't include current lap
                %cfg_encLaps = []; cfg_encLaps.mode = 'dist';
                %enc_laps = get_LOOCV_laps(cfg_encLaps,expCond(iCond).t_idx,iLap,this_dist);
                enc_laps = find(expCond(iCond).t_idx == (this_lap + this_dist));
                
                if isempty(enc_laps)
                    fprintf('Dist unavailable, skipping...\n');
                    continue;
                end
                
                cfg_select = []; cfg_select.verbose = 0;
                lap_iv = SelectIV(cfg_select,expCond(iCond).t,enc_laps);
                
                enc_S = restrict(expCond(iCond).S,lap_iv);
                enc_linpos = restrict(expCond(iCond).linpos,lap_iv);
                
                % could be empty -- no cells for this lap's encoding model
                if all(cellfun(@isempty,enc_S.t))
                   continue; 
                end
                
                expCond(iCond).tc = TuningCurves(cfg_tc,enc_S,enc_linpos);
                           
                %% Q-mat
                this_Qsd = cfg_master.QsmoothSD;
                
                cfg_Q = [];
                cfg_Q.dt = cfg_master.dt;
                
                if this_Qsd == 0
                    cfg_Q.smooth = [];
                else
                    cfg_Q.smooth = 'gauss';
                    cfg_Q.gausswin_sd = this_Qsd;
                end
                    
                expCond(iCond).Q = MakeQfromS(cfg_Q,expCond(iCond).S);
                
                % restrict Q-mat with same intervals (run, position, etc)
                % as everything else... track_iv, run_iv, trial times
                
                expCond(iCond).Q = restrict(expCond(iCond).Q,track_iv1);
                expCond(iCond).Q = restrict(expCond(iCond).Q,track_iv2);
                expCond(iCond).Q = restrict(expCond(iCond).Q,run_iv);
                expCond(iCond).Q = restrict(expCond(iCond).Q,expCond(iCond).t);
                
                if ~isempty(chew_iv)
                    expCond(iCond).Q = restrict(expCond(iCond).Q,chew_iv);
                end
                
                %% decode
                cfg_temp.verbose = 0;
                dec_lap_iv = SelectIV(cfg_temp,expCond(iCond).t,iLap);
                dec_Q = restrict(expCond(iCond).Q,dec_lap_iv);
                
                dec_S = restrict(expCond(iCond).S,dec_lap_iv);
                dec_linpos = restrict(expCond(iCond).linpos,dec_lap_iv);
                    
                % could be that at this point there is no data left!
                if isempty(dec_Q.data) | isempty(dec_linpos.data)
                    fprintf('--> WARNING: skipped lap, no decoding spikes left!\n');
                    continue;
                end
                
                dec_tc = TuningCurves(cfg_tc,enc_S,dec_linpos); % to get true Z pos later, not for TCs
                
                cfg_decode = [];
                cfg_decode.nMinSpikes = cfg_master.dt;
                this_lap_P = DecodeZ(cfg_decode,dec_Q,expCond(iCond).tc.tc);
                
                % true position
                this_trueZ = tsd(dec_linpos.tvec,dec_tc.pos_idx);
                
                % append to tsd
                this_Psplice{iDist} = UnionTSD([],this_Psplice{iDist},this_lap_P);
                this_trueZsplice{iDist} = UnionTSD([],this_trueZsplice{iDist},this_trueZ); % could move outside Q-loop
                                    
            end % of encoding lap loop
            
        end % of decoding lap lap loop
        
        % now need to compute error over all laps (for each iEncLap)
        cfg_err = [];
        cfg_err.mode = 'max';
        
        for iDist = 1:length(cfg_master.nEncLaps)
            
            if isempty(this_Psplice{iDist}.tvec)
                % no data available, skip
                continue;
            end
            
            [expCond(iCond).Perr,expCond(iCond).confMat] = DecodeErrorZ(cfg_err,this_Psplice{iDist},this_trueZsplice{iDist});
            
            % decoding error is in bins, so convert to cm
            this_binSize = ExpKeys.pathlength./cfg_master.nBins;
            expCond(iCond).Perr.data = expCond(iCond).Perr.data*this_binSize + this_binSize/2;
            
            % find out if any time bins that should have been decoded
            % failed
            goodPbins = logical(this_Psplice{iDist}.usr.nActiveNeuronsPassed);
            expCond(iCond).Pnan = sum(isnan(nansum(this_Psplice{iDist}.data(:,goodPbins))));
            expCond(iCond).Pnan = expCond(iCond).Pnan./sum(goodPbins);
            expCond(iCond).Ppassed = sum(goodPbins)./length(goodPbins);
            
            
            % keep track -- add to ALL_decErr
            ALL_decErr.(expCond(iCond).label).Pnan(iFD,iDist) = expCond(iCond).Pnan;
            ALL_decErr.(expCond(iCond).label).Ppassed(iFD,iDist) = expCond(iCond).Ppassed;
            
            ALL_decErr.(expCond(iCond).label).meanErr(iFD,iDist) = nanmean(expCond(iCond).Perr.data);
            
            ALL_decErr.(expCond(iCond).label).trueZ{iFD} = this_trueZsplice{iDist};
            
            cfg_tsdZ = []; cfg_tsdZ.edges = cfg_tc.binEdges{1};
            ALL_decErr.(expCond(iCond).label).spaceErr(iFD,iDist,:) = TSDbyZ(cfg_tsdZ,expCond(iCond).linpos,expCond(iCond).Perr);
            
            ALL_decErr.(expCond(iCond).label).nSpikesHist(iFD,iDist,:) = histc(this_Psplice{iDist}.usr.nActiveNeurons,cfg_master.nSpikesHist);
            
            ALL_decErr.(expCond(iCond).label).confMat{iFD,iDist} = expCond(iCond).confMat;
            
        end % of encoding lap loop
        
        % track S
        ALL_decErr.(expCond(iCond).label).S{iFD} = expCond(iCond).S;
        
    end % of iCond loop
        
end % of sessions

%% do some housekeeping    
% add configs etc
ALL_decErr.fd = fd;
ALL_decErr.cfg = cfg_master;