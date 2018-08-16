% basic batch script to generate decoding for sequence detection

%%
cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'prerecord','postrecord','goodTheta','pathlength'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'coord','taskvars'};
cfg.requireCandidates = 1;
cfg.requireVT = 1;
cfg.requireTimes = 1;
cfg.requireFiles = 1;
cfg.rats = {'R042','R044','R050','R064'};
cfg.excludeSessions = [1 7 8 9 12]; % sessions with insufficient cells
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));
    fd(cfg.excludeSessions) = [];
    
    cfg_decSeq = [];
    cfg_decSeq.minSeqLength = 10; % minimum sequence length in time steps (cfg_def.Qdt in Generate_DecSeqShuf)
    cfg_decSeq.removeInterneurons = 0;
    cfg_decSeq.nMinNeurons = 4; % minimum number of neurons that must be active for posterior to be retained (set to NaN otherwise)
    cfg_decSeq.nMaxNanSkipSequential = 0; % number of consecutive bins that can be NaN
    cfg_decSeq.output_file_prefix = 'S1_';
    cfg_decSeq.nShuffles = 500;
    
    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        
        Generate_DecSeqShuf(cfg_decSeq);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')