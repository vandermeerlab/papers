function [ All_naris ] = Naris_comp_plot(cfg)
%NARIS_COMP_PLOT(): this will take the "naris" data files generated in the
%initial processing of the naris data and then normalize the data to the
%pre record.  A figure is then made with the four subplots that
%corresponding to the four recording days.

%INPUTS:
%     naris [struct]: output from Naris_analysis_fast or
%     Naris_analysis_sandbox
%
%OUTPUTS:
%     all_naris_data [struct]: This is a structre with all the naris data
%     compiled together.

%% Initialize which sessions to use
og_cd = cd;

Naris.R054.ses_dir = 'G:\Naris\R054';
Naris.R054.sessions = {'R054-2014-12-06', 'R054-2014-12-07','R054-2014-12-09','R054-2014-12-10'};
Naris.R054.phases = {'pre', 'ipsi', 'contra', 'post'};
Naris.R054.chan_to_plot = 34;

Naris.R053.ses_dir = 'G:\Naris\R053';
Naris.R053.sessions = {'R053-2014-12-27', 'R053-2014-12-28','R053-2014-12-30','R053-2014-12-31'};
Naris.R053.phases = {'pre', 'right', 'left', 'post'};

ids = fieldnames(Naris);
%% load the data into the Naris structure.
if exist('G:\Naris\all_Naris.mat') ==0
    for id = 1%:length(fieldnames(Naris)) % loop over the rat numbers
        Naris.(ids{id}).data = [];
        for nSess = 1:length(Naris.(ids{id}).sessions) % loop over the session files
            Naris.(ids{id}).data.(strrep(Naris.(ids{id}).sessions{nSess}, '-', '_')) = [];
            for iphase = 1:length(Naris.(ids{id}).sessions) % loop over the
                cd([Naris.(ids{id}).ses_dir '\' Naris.(ids{id}).sessions{nSess}])
                if strcmp(Naris.(ids{id}).ses_dir(end-3:end), 'R054') ==0
                    temp_data = load('Naris_data.mat');
                    Naris.(ids{id}).data.(strrep(Naris.(ids{id}).sessions{nSess}, '-', '_')).(Naris.(ids{id}).phases{iphase}) = temp_data.naris.(Naris.(ids{id}).phases{iphase});
                else
                    temp_data = load(['Naris_' (Naris.(ids{id}).phases{iphase}) '_data_preprocess.mat']);
                    Naris.(ids{id}).data.(strrep(Naris.(ids{id}).sessions{nSess}, '-', '_')).(Naris.(ids{id}).phases{iphase}) =  temp_data.data;
                end
                disp([Naris.(ids{id}).sessions{nSess} ' ' (Naris.(ids{id}).phases{iphase})])
            end
        end
    end
    
    disp('Finished')
end
%     [cfg] = Naris_cfgs(cfg.file_name(1:4));
%     if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R060')
%     cfg.Naris_exp = {'pre', 'right', 'left', 'post'};
%     tetrodes = tts(cfg);
%     cfg.tetrodes = tetrodes.channels_array;

%% get the mean and std for each phase per rat. 


%% Extract and excompile all the data
% for id = 1:length(fieldnames(Naris))
%     for nSess = 1:length(Naris.(ids{id}).sessions)
%         cd([Naris.(ids{id}).ses_dir '\' Naris.(ids{id}).sessions{nSess}])
%         Naris.(ids{id}).sessions{nSess}
%         all_naris_data.(strrep(Naris.(ids{id}).sessions{nSess}, '-', '_')) = load(['Naris_' Naris.(ids{id}).sessions '_preprocessed.mat']);
%         disp([sessions{nSess} ' loaded'])
%     end
%
%     %% normalize the data to the pre session
%     for nSess = 1:length(sessions)
%         temp_data = all_naris_data.(strrep(sessions{nSess}, '-', '_'));
%         %     norm_data.pre = (temp_data.naris.pre.data{1}.Data)./temp_data.naris.pre.data{1}.Data;
%         for nPhase = 2:4
%             norm_data.(phases{nPhase}) = (temp_data.naris.(phases{nPhase}).data{1}.Data)./temp_data.naris.pre.data{1}.Data;
%         end
%
%     end
%
%     %% plot everything
%     figure(100)
%     hold on
%     maximize
%     freq= temp_data.naris.(phases{nPhase}).data{1}.Frequencies;
%     plot(freq, 10*log10(norm_data.right), freq, 10*log10(norm_data.left), freq, 10*log10(norm_data.post))
%     legend ('right', 'left', 'post')
%     xlim([0 100])
%     vline([45 60], {'g','g'}, {'Gamma50', ' '})
%     ylim([-6 6])
%     %%
%
%
%


end

