% batch script to generate non-sequence decoding output for each session
%
% the only configuration setting you may want to change here is:
%
% cfg_decSeq.removeInterneurons: set to 0 (default) to include all cells for
% decoding, set to 1 to remove putative interneurons (> 5 Hz firing rate)

%%
cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'prerecord','postrecord','goodTheta','pathlength'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'coord','taskvars'};
cfg.requireCandidates = 1;
cfg.requireVT = 1;
cfg.requireTimes = 1;
cfg.requireFiles = 0;
cfg.rats = {'R042','R044','R050','R064'};
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));
    
    cfg_decSeq = [];
    cfg_decSeq.removeInterneurons = 0;
    
    for iFD = length(fd):-1:1
        
        cd(fd{iFD});
        
        out{iFD} = Generate_DecSeqCombined(cfg_decSeq);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')