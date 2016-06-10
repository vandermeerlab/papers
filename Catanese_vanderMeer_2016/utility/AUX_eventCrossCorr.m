%% load ALL_evt variable first

%%
rats = {'R026','R032','R033','R039'};
cfg_ccf = [];
cfg_ccf.binsize = 0.025;

ALL_ccf = [];
%%
available_rats = fieldnames(ALL_evt);
nSession = 1;

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

        % remove overlap
        this_lg = DifferenceIV([],this_session_data.lg,this_session_data.hg);
        this_hg = DifferenceIV([],this_session_data.hg,this_session_data.lg);
        
        lg_t = IVcenters(this_lg);
        hg_t = IVcenters(this_hg);

        [this_ccf,xbins] = ccf(cfg_ccf,lg_t,hg_t);
        
        if nSession == 1
            ALL_ccf = this_ccf;
        else
            ALL_ccf = ALL_ccf + this_ccf;
        end
        
        subplot(3,3,nSession);
        bar(xbins,this_ccf);
        
        nSession = nSession + 1;
    end
    
end

subplot(3,3,9);
bar(xbins,ALL_ccf/(nSession-1));