%% for Julien's data set
rats = {'R026','R032','R033','R039'};
available_rats = fieldnames(ALL_evt);

%% for ADRLAB
rats = {'R117','R119','R131','R132'};
available_rats = fieldnames(ALL_evt);

%%
total_lg = 0;
total_hg = 0;
nSession = 0;
for iRat = 1:length(rats)
    
    this_rat = rats{iRat};
    
    if ~strmatch(this_rat,available_rats)
       warning('Rat %s not available -- skipping...',rats{iRat});
       continue;
    end
    
    available_sessions = fieldnames(ALL_evt.(this_rat));
    
    for iSession = 1:length(available_sessions)
        nSession = nSession + 1;
        
        this_session = available_sessions{iSession};
        this_session_data = ALL_evt.(this_rat).(this_session);
        
        fprintf('Session %s: %d lg, %d hg events.\n',this_session,length(this_session_data.lg.tstart),length(this_session_data.hg.tstart));
        total_lg = total_lg + length(this_session_data.lg.tstart);
        total_hg = total_hg + length(this_session_data.hg.tstart);
    
    end % over sessions
    
end

fprintf('TOTAL %d sessions, %d lg, %d hg events.\n',nSession,total_lg,total_hg);


%%