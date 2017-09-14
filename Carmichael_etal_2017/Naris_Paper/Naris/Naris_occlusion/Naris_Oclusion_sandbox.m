%% Naris Oclusion sandbox
% Naris.R053.sessions = {'R053-2014-12-27', 'R053-2014-12-28', 'R053-2014-12-30', 'R053-2014-12-31'};
% Naris.R054.sessions = {'R054-2014-12-06', 'R054-2014-12-07', 'R054-2014-12-09', 'R054-2014-12-10'};
% Naris.R060.session_list = {'R053-2014-12-27', 'R053-2014-12-28', 'R053-2014-12-30', 'R053-2014-12-31'};
load('G:\Naris\Naris_psd_white.mat')

ids = fieldnames(Naris)

%%
for id = 1:length(ids)
    sessions = fieldnames(Naris.(ids{id}).data);
    for iSess = 1:length(sessions)
        if ispc
            cd(['G:\Naris\' (ids{id}) '\' Naris.(ids{id}).sessions{iSess}])
        elseif isunix
            cd(['/Users/jericcarmichael/Documents/Nairs_data/' Naris.(ids{id}).sessions{iSess}])
        end
        cfg = [];
        
        cfg.tetrodes = Naris.(ids{id}).data.(sessions{iSess}).pre.cfg.tetrodes;
        Naris_gamma_count(cfg)
    end
end