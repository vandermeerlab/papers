% batch script to generate decoding output for each session

%%
cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'prerecord','postrecord','goodTheta','pathlength'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'coord','taskvars'};
cfg.requireCandidates = 1;
cfg.requireVT = 1;
cfg.requireTimes = 1; %R042 only
cfg.requireFiles = 0;
cfg.rats = {'R042','R044','R050','R064'};
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));
    
    cfg_decSeq = [];
    cfg_decSeq.removeInterneurons = 1;
    
    for iFD = length(fd):-1:1
        
        cd(fd{iFD});
        
        out{iFD} = Generate_DecSeqCombined(cfg_decSeq);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')