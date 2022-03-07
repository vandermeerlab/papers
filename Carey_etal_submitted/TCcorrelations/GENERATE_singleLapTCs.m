% GENERATE_singleLapTCs
%
% build single lap tuning curves, place into ALL_TC variable
%
% analyze output with COLLECT_singleLapTCs.m

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
cfg_master.fd_skip = [1 7 8 9 12]; % sessions to skip because not enough cells; cf. van der Meer et al. (2017) Hippocampus

%% find data folders eligible for analysis
fd_cfg = []; fd_cfg.requireCandidates = 0; % don't need previously saved sharp wave-ripple (SWR) candidates here
fd_cfg.requireEvents = 1;
fd = getTmazeDataPath(fd_cfg);

%% main loop across sessions
for iFD = 1:length(fd)
    
    clear expCond;
    
    expCond(1).label = 'left';
    expCond(2).label = 'right';
    
    % skip sessions that aren't included in final analysis anyway
    if any(cfg_master.fd_skip == iFD)
        fprintf('Session skipped as per cfg.\n');
        continue;
    end
    
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
    %[left,right] = GetMatchedTrials([],metadata,ExpKeys);
    %expCond(1).t = left;
    %expCond(2).t = right;
    
    metadata = RemoveBadTrials([],metadata,ExpKeys);
    expCond(1).t = metadata.taskvars.trial_iv_L;
    expCond(2).t = metadata.taskvars.trial_iv_R;
    
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
        
        %if ~isempty(chew_iv)
        %    fh = @(x) restrict(x,chew_iv); % restrict S and linpos to non-detached times
        %    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
        %end
        
        % also remove cells with insufficient spikes
        expCond(iCond).S = removeEmptyCells(expCond(iCond).S);
        
        spk_count = getSpikeCount([],expCond(iCond).S);
        cell_keep_idx = spk_count >= cfg_master.minSpikes;
        
        expCond(iCond).S = SelectTS([],expCond(iCond).S,cell_keep_idx);
    end
    
    %% run through conditions (L,R) and parameters requested
    for iCond = 1:nCond
        
        nLaps = length(expCond(iCond).t.tstart);
        
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
                cfg_select = []; cfg_select.verbose = 0;
                lap_iv = SelectIV(cfg_select,expCond(iCond).t,iLap);
                
                enc_S = restrict(expCond(iCond).S,lap_iv);
                enc_linpos = restrict(expCond(iCond).linpos,lap_iv);

                expCond(iCond).tc{iLap} = TuningCurves(cfg_tc,enc_S,enc_linpos);
                expCond(iCond).tc{iLap}.label = enc_S.label;
                
            end % of TC smoothing
            
        end % of lap loop
        
        % track TCs across sessions
        ALL_TC{iFD}.(expCond(iCond).label) = expCond(iCond).tc;
        
    end % of iCond loop
        
end % of sessions