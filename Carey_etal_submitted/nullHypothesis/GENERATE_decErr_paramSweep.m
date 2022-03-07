%% preamble
clear all; pack

MASTER_root = 'C:\Users\mvdm\Dropbox\projects\Alyssa'; % replace this with home folder of your project
cd(MASTER_root);
MASTER_path; % reset and then set up path for this project

%% master config
cfg_master = [];
cfg_master.dt = 0.025; % time bin (tau) in s
cfg_master.TCsmooth = 1; % SD of Gaussian smoothing kernel for tuning curves
cfg_master.QsmoothSD = 0.002; % SD of Gaussian kernel for SDFs
cfg_master.minSpikes = 25; % remove cells with less than this number of spikes during run
cfg_master.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
cfg_master.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
cfg_master.load_questionable_cells = 1;
cfg_master.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
cfg_master.trackExcludeEnd = 15;
cfg_master.nSpikesHist = -0.5:105; % used to keep histogram of number of spikes per time bin

nMaxLaps = 20;
cfg_master.encdecmat = 1-eye(nMaxLaps); % leave-one-out

%% find data folders eligible for analysis
fd_cfg = []; fd_cfg.requireCandidates = 0; % don't need previously saved sharp wave-ripple (SWR) candidates here
fd_cfg.requireEvents = 1;
fd = getTmazeDataPath(fd_cfg);

%% init some vars collecting output across sessions
clear expCond;

expCond(1).label = 'left';
expCond(2).label = 'right';
    
for iCond = 1:length(expCond)
    ALL_decErr.(expCond(iCond).label).meanErr = nan(length(fd),length(cfg_master.TCsmooth),length(cfg_master.QsmoothSD));
    ALL_decErr.(expCond(iCond).label).spaceErr = nan(length(fd),length(cfg_master.TCsmooth),length(cfg_master.QsmoothSD),cfg_master.nBins-1);
end

%% main loop across sessions
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
        
    %% set up data structs for L, R data    
    [left,right] = GetMatchedTrials([],metadata,ExpKeys);
    expCond(1).t = left;
    expCond(2).t = right;
    
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
    cfg_track2 = []; cfg_track2.method = 'raw'; cfg_track2.operation = '<'; cfg_track2.threshold = ExpKeys.pathlength - cfg_master.trackExcludeEnd;
    
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
            
            cfg_det = []; cfg_det.threshold = 0.5; cfg_det.operation = '>'; cfg_det.method = 'raw';
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
    
    %% run through conditions (L,R) and parameters requested
    for iCond = 1:nCond
        
        nLaps = length(expCond(iCond).t.tstart);
        
        % first, initialize P (decoding posterior) and trueZ (actual linearlized position)
        for iTCsmooth = 1:length(cfg_master.TCsmooth)
            for iQ = 1:length(cfg_master.QsmoothSD)
                this_Psplice{iTCsmooth,iQ} = tsd; 
                this_Psplice{iTCsmooth,iQ}.usr.nActiveNeurons = [];
                this_Psplice{iTCsmooth,iQ}.usr.nActiveNeuronsPassed = [];
                
                this_trueZsplice{iTCsmooth,iQ} = tsd;
            end
        end
        
        for iLap = 1:nLaps
            
            fprintf('.lap %d/%d...\n',iLap,nLaps);
                        
            %%
            for iTCsmooth = 1:length(cfg_master.TCsmooth)
                
                this_TCsmooth = cfg_master.TCsmooth(iTCsmooth);
                fprintf('..TC smoothing %d/%d...\n',iTCsmooth,length(cfg_master.TCsmooth));
                
                %% construct tuning curves
                cfg_tc = [];
                cfg_tc.binEdges{1} = linspace(0,ExpKeys.pathlength,cfg_master.nBins); % ~3cm bins for full maze (smaller for R042...)
                
                if this_TCsmooth ~= 0
                    cfg_tc.smoothingKernel = gausskernel(51,this_TCsmooth);
                end
                
                % restrict to the encoding model for this lap
                enc_laps = cfg_master.encdecmat(iLap,1:length(expCond(iCond).t.tstart));
                cfg_select = []; cfg_select.verbose = 0;
                lap_iv = SelectIV(cfg_select,expCond(iCond).t,find(enc_laps));
                
                enc_S = restrict(expCond(iCond).S,lap_iv);
                enc_linpos = restrict(expCond(iCond).linpos,lap_iv);
                
                % could be empty -- no cells for this lap's encoding model, then skip
                if all(cellfun(@isempty,enc_S.t))
                   continue; 
                end
                
                expCond(iCond).tc = TuningCurves(cfg_tc,enc_S,enc_linpos);
                
                
                %% Q smoothing loop
                for iQ = 1:length(cfg_master.QsmoothSD)
                    
                    this_Qsd = cfg_master.QsmoothSD(iQ);
                    fprintf('...Q smoothing %d/%d...\n',iQ,length(cfg_master.QsmoothSD));
                    
                    %% Q-mat
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
                    
                    % could be that at this point there is no data left,
                    % then skip
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
                    this_Psplice{iTCsmooth,iQ} = UnionTSD([],this_Psplice{iTCsmooth,iQ},this_lap_P);
                    this_trueZsplice{iTCsmooth,iQ} = UnionTSD([],this_trueZsplice{iTCsmooth,iQ},this_trueZ); % could move outside Q-loop
                            
                    
                end % of Q-mat smoothing
                
            end % of TC smoothing
            
        end % of lap loop
        
        % now need to compute error over all laps (for each iTCsmooth, iQ pair)
        cfg_err = [];
        cfg_err.mode = 'max';
            
        for iTCsmooth = 1:length(cfg_master.TCsmooth)
            for iQ = 1:length(cfg_master.QsmoothSD)
                
                [expCond(iCond).Perr,expCond(iCond).confMat] = DecodeErrorZ(cfg_err,this_Psplice{iTCsmooth,iQ},this_trueZsplice{iTCsmooth,iQ});
                
                % decoding error is in bins, so convert to cm
                this_binSize = ExpKeys.pathlength./cfg_master.nBins;
                expCond(iCond).Perr.data = expCond(iCond).Perr.data*this_binSize + this_binSize/2;
                
                % find out if any time bins that should have been decoded
                % failed
                goodPbins = logical(this_Psplice{iTCsmooth,iQ}.usr.nActiveNeuronsPassed);
                expCond(iCond).Pnan = sum(isnan(nansum(this_Psplice{iTCsmooth,iQ}.data(:,goodPbins))));
                expCond(iCond).Pnan = expCond(iCond).Pnan./sum(goodPbins);
                expCond(iCond).Ppassed = sum(goodPbins)./length(goodPbins);
                
                % keep track -- add to ALL_decErr
                ALL_decErr.(expCond(iCond).label).Pnan(iFD,iTCsmooth,iQ) = expCond(iCond).Pnan;
                ALL_decErr.(expCond(iCond).label).Ppassed(iFD,iTCsmooth,iQ) = expCond(iCond).Ppassed;
                ALL_decErr.(expCond(iCond).label).meanErr(iFD,iTCsmooth,iQ) = nanmean(expCond(iCond).Perr.data);
                ALL_decErr.(expCond(iCond).label).trueZ{iFD} = this_trueZsplice{iTCsmooth,iQ};
                
                cfg_tsdZ = []; cfg_tsdZ.edges = cfg_tc.binEdges{1};
                ALL_decErr.(expCond(iCond).label).spaceErr(iFD,iTCsmooth,iQ,:) = TSDbyZ(cfg_tsdZ,expCond(iCond).linpos,expCond(iCond).Perr);
                ALL_decErr.(expCond(iCond).label).nSpikesHist(iFD,iQ,:) = histc(this_Psplice{iTCsmooth,iQ}.usr.nActiveNeurons,cfg_master.nSpikesHist);
                ALL_decErr.(expCond(iCond).label).confMat{iFD,iTCsmooth,iQ} = expCond(iCond).confMat;
                
            end
        end
        
        % track S
        ALL_decErr.(expCond(iCond).label).S{iFD} = expCond(iCond).S;
        
    end % of iCond loop
        
end % of sessions

%% do some housekeeping    
% add configs etc
ALL_decErr.fd = fd;
ALL_decErr.cfg = cfg_master;