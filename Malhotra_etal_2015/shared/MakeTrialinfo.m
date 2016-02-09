function trialinfo = MakeTrialinfo(cfg_in)
% function trialinfo = MakeTrialinfo(cfg)

cfg.writeOutput = 1;
ProcessConfig;

%% parameters
x_center = 320;
x_left = 150;
x_right = 490;

%% load some data

% keys
run(FindFile('*keys.m'));

% evt

fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
EVTimeStamps = EVTimeStamps * 10^-6;
EVTimeStamps = EVTimeStamps - ExpKeys.TimeOnTrack(1); % get rid of absolute event times

% vt

load(FindFile('*vt.mat'));
xd = Data(x);
xr = Range(x)-ExpKeys.TimeOnTrack(1);

%% determine event strings to detect

% first check if an unused photobeam is on (from Rob's task, which ran in
% parallel with Sushant's LinearTrackTone)

idx = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004)',EventStrings);

if ~isempty(idx)
    fprintf('NOTICE: photobeam 3 detected. Verifying data consistency...\n');
    
    idx2 = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0001)',EventStrings);
    idx3 = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0002)',EventStrings);
    
    if ~isempty(idx2) | ~isempty(idx3)
        error('0x0001 or 0x0002 events found; should not occur if photobeam 3 is on!');
    else
        fprintf('NOTICE: events OK, no 0x0001 or 0x0002 events found.\n');
    end
    
    pboff_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004)';
    pb0_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0005)';
    pb1_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0006)';
else
    pboff_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000)';
    pb0_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0001)';
    pb1_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0002)';
end

foff_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0000)';
f0_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0001)';
f1_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0002)';

c1_string = '1 pellet cue';
c3_string = '3 pellet cue';
c5_string = '5 pellet cue';
lo_string = '2 or 4 pellet cue';
hi_string = '1 or 5 pellet cue';

%% step 1: assemble raw event times
all_t.pb0 = EVTimeStamps(strmatch(pb0_string,EventStrings));
all_t.pb1 = EVTimeStamps(strmatch(pb1_string,EventStrings));

all_t.f0 = EVTimeStamps(strmatch(f0_string,EventStrings));
all_t.f1 = EVTimeStamps(strmatch(f1_string,EventStrings));

all_t.c1 = EVTimeStamps(strmatch(c1_string,EventStrings));
all_t.c3 = EVTimeStamps(strmatch(c3_string,EventStrings));
all_t.c5 = EVTimeStamps(strmatch(c5_string,EventStrings));
all_t.lo = EVTimeStamps(strmatch(lo_string,EventStrings));
all_t.hi = EVTimeStamps(strmatch(hi_string,EventStrings));

%% step 2: detect x-coordinate crossings (virtual photobeams)
in_left_zone = xd < x_left;
all_t.l_enter = xr(find(diff(in_left_zone) == 1));
all_t.l_exit = xr(find(diff(in_left_zone) == -1));

in_right_zone = xd > x_right;
all_t.r_enter = xr(find(diff(in_right_zone) == 1));
all_t.r_exit = xr(find(diff(in_right_zone) == -1));

%% some useful definitions

all_cues = {'c1','c3','c5','lo','hi'}; 
cfg_next_cue = []; cfg_next_cue.mode = 'next'; cfg_next_cue.fields = all_cues;

all_feeders = {'f0','f1'}; 
cfg_next_reward = []; cfg_next_reward.mode = 'next'; cfg_next_reward.fields = all_feeders;
cfg_prev_reward = []; cfg_prev_reward.mode = 'prev'; cfg_prev_reward.fields = all_feeders;

all_pb = {'pb0','pb1'}; 
cfg_next_pb = []; cfg_next_pb.mode = 'next'; cfg_next_pb.fields = all_pb;

all_enters = {'l_enter','r_enter'}; 
cfg_next_enter = []; cfg_next_enter.mode = 'next'; cfg_next_enter.fields = all_enters;

all_exits = {'l_exit','r_exit'}; 
cfg_prev_exit = []; cfg_prev_exit.mode = 'prev'; cfg_prev_exit.fields = all_exits;

%%%%%%%%%%%%%%%%%
%%% MAIN LOOP %%%
%%%%%%%%%%%%%%%%%
%%

trial_count = 0;
warning_count = 0;
prev_lr = '';

[next_cue_t,next_cue_ID] = FindFieldTime(cfg_next_cue,all_t,0);

clear trialinfo;
while ~isempty(next_cue_t)
    
    trial_count = trial_count + 1;
    
    % pick up current cue
    curr_t = next_cue_t; curr_ID = next_cue_ID{1};
    fprintf('** Trial %d (%s, t = %.2f)\n',trial_count,curr_ID,curr_t);
    
    trialinfo.cue_time(trial_count) = curr_t;
    trialinfo.cue_ID{trial_count} = curr_ID;
   
    % find x-coordinate
    curr_x = interp1(xr,xd,curr_t,'linear');
    
    if curr_x < x_center
        curr_lr = 'R';
    else
        curr_lr = 'L';
    end
    
    [next_enter_t,next_enter_ID] = FindFieldTime(cfg_next_enter,all_t,curr_t);
    if ~isempty(next_enter_ID)
        fprintf('  xcoord is %.1f, so trial side is %s (next enter: %s)\n',curr_x,curr_lr,next_enter_ID{1});
    else
        fprintf('  xcoord is %.1f, so trial side is %s (next enter: n/a)\n',curr_x,curr_lr);
    end
    
    if strcmp(curr_lr,prev_lr)
        y = input(sprintf('!! Previous trial is same (%s); override current trial? (Y/N) ',prev_lr),'s');
        if strcmp(y,'Y')
            switch curr_lr
                case 'R'
                    curr_lr = 'L';
                case 'L'
                    curr_lr = 'R';
            end
            fprintf('Overridden.\n');
        end
    end
    
    % set x-dist
    switch curr_lr
        case 'L'
            trialinfo.xcoord(trial_count) = curr_x - x_left;
        case 'R'
            trialinfo.xcoord(trial_count) = x_right - curr_x;
    end
    
    % find next pb
    [next_pb_t,next_pb_ID] = FindFieldTime(cfg_next_pb,all_t,curr_t);
    if isempty(next_pb_t)
        fprintf('   ...no next pb found.\n');
        trialinfo.pb_time(trial_count) = NaN;
    else
        switch curr_lr
            case 'R'
                if strcmp(next_enter_ID{1},'r_enter')
                    fprintf('  expected zone entered (%s triggered next)\n',next_pb_ID{1});
                    trialinfo.pb_time(trial_count) = next_pb_t;
                else
                    warning('!! unexpected zone entered (%s triggered next)\n',next_pb_ID{1});
                    warning_count = warning_count + 1;
                    trialinfo.pb_time(trial_count) = NaN;
                end
            case 'L'
                if strcmp(next_enter_ID{1},'l_enter')
                    fprintf('  expected zone entered (%s triggered next)\n',next_pb_ID{1});
                    trialinfo.pb_time(trial_count) = next_pb_t;
                else
                    warning('!! unexpected zone entered (%s triggered next)\n',next_pb_ID{1});
                    warning_count = warning_count + 1;
                    trialinfo.pb_time(trial_count) = NaN;
                end
        end
    end
    
    % get previous exit and do some checks
    [prev_exit_t,prev_exit_ID] = FindFieldTime(cfg_prev_exit,all_t,curr_t);
    if isempty(prev_exit_t)
       fprintf('   ...no previous exit found.\n'); 
       trialinfo.prev_exit_t(trial_count) = NaN;
    else
        switch curr_lr
            case 'R'
                if strcmp(prev_exit_ID{1},'l_exit')
                    fprintf('  expected zone exited.\n');
                    trialinfo.prev_exit_t(trial_count) = prev_exit_t;
                else
                    warning('!! unexpected zone exited');
                    warning_count = warning_count + 1;
                    trialinfo.prev_exit_t(trial_count) = NaN;
                end
            case 'L'
                if strcmp(prev_exit_ID{1},'r_exit')
                    fprintf('  expected zone exited.\n');
                    trialinfo.prev_exit_t(trial_count) = prev_exit_t;
                else
                    warning('!! unexpected zone exited');
                    warning_count = warning_count + 1;
                    trialinfo.prev_exit_t(trial_count) = NaN;
                end
        end
    end
    
    % check if reward happens before the next cue (correct trial if so)
    [next_cue_t,next_cue_ID] = FindFieldTime(cfg_next_cue,all_t,curr_t);
    [next_reward_t,next_reward_ID] = FindFieldTime(cfg_next_reward,all_t,curr_t);
    
    if isempty(next_cue_t) % last cue
        
        if isempty(next_reward_t)
            correct_trial = 0; next_reward_t = NaN;
        else
            correct_trial = 1;
        end
        
    else % not the last cue
        
        if next_reward_t < next_cue_t
            correct_trial = 1;
        else
            correct_trial = 0; next_reward_t = NaN;
        end
        
    end
    
    if correct_trial
        trialinfo.reward_time(trial_count) = next_reward_t;
    else
        trialinfo.reward_time(trial_count) = NaN;
    end
        
    trialinfo.correct(trial_count) = correct_trial;
    
    % get previous reward
    prev_reward_t = FindFieldTime(cfg_prev_reward,all_t,curr_t);
    if isempty(prev_reward_t) % first reward
        trialinfo.prev_reward_t(trial_count) = NaN;
    else
        trialinfo.prev_reward_t(trial_count) = prev_reward_t;
    end
    
    prev_lr = curr_lr;
end

%% save
if cfg.writeOutput
   [~,fd,~] = fileparts(pwd);
   fn_out = cat(2,fd,'_trialinfo.mat');
   save(fn_out,'trialinfo');
    
end