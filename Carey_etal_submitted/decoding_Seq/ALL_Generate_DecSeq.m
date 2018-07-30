% basic batch script to generate sequence detection based on decoded data

%%
cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'prerecord','postrecord','goodTheta','pathlength'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'coord','taskvars'};
cfg.requireCandidates = 1;
cfg.requireVT = 1;
cfg.requireTimes = 1; %R042 only
cfg.requireFiles = 1;
cfg.rats = {'R042','R044','R050','R064'};
cfg.excludeSessions = [1 7 8 9 12]; % sessions with insufficient cells
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));
    fd(cfg.excludeSessions) = [];
    
    cfg_decSeq = [];
    cfg_decSeq.minSeqLength = 10;
    cfg_decSeq.removeInterneurons = 0;
    cfg_decSeq.nMinNeurons = 4;
    cfg_decSeq.nMaxNanSkipSequential = 0;
    cfg_decSeq.output_file_prefix = 'S0_';
    cfg_decSeq.nShuffles = 1000;
    
    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        
        Generate_DecSeqShuf(cfg_decSeq);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')