%%
rats = {'R026','R032','R033','R039'};
task_time = 0; rest_time = 0;

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
        
        cd(this_session_data.fd);
        LoadExpKeys;
        
        t_track = ExpKeys.TimeOffTrack-ExpKeys.TimeOnTrack;
        t_rest = (ExpKeys.TimeOffPre-ExpKeys.TimeOnPre) + (ExpKeys.TimeOffPost-ExpKeys.TimeOnPost);
        
        task_time = task_time + t_track;
        rest_time = rest_time + t_rest;
        
    end
    
end