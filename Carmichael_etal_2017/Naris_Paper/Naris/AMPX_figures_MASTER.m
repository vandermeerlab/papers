%% AMPX_Naris_generate_figures: Generates all the figures used in Carmichael et al. 
%   
%     - figure 1
%         - C: raw traces of a random low gamma event
%         - D: heat plot of the gamma power during the same random low gamma event
%         - E: 4 x 2 of the average high and low gamma power across
%              sessions per rat (1-4)
%         - S1: same as E but the color axis scale is the min/max across
%              all sessions/rats
%
%     - figure 2 [AMPX_Naris_fig_2_task.m]
%         - A: generates figure 2 which contains the heatmaps
%              for the average rewarded/approach gamma 50 and gamma 80 
%              events for the two rats that reached task criterion
%         - B: plane fitting histograms and statistics (R^2 and
%              angles/directions)
%
%     - figure 3 [AMPX_Naris_fig_3_phase.m]
%         - A: creates 4 x 2 figues of the phase differences
%              for both across entire gamma events 
%         - B: 4 x 2 for phase differences across the central three cycles.
%     
%     - figure 4 [AMPX_Naris_fig_4_CSD.m]
%         - A: kernal current source density for the example event used in
%              fig 1E
%         - B: 4 x 2 of the average kCSD for the central three cycles in
%              each event per session per rat
%
%     - figure 5 [AMPX_Naris_fig_5_gamma_count.m]
%         - bar plot of the normalized count of gamma events across all
%           Naris occlusion recordings
%         - figure
%
%
% EC - 2016-06-14

global PARAMS

%% load the pre and post data:
load([PARAMS.intermediate_dir '\Naris_all_data_pre.mat']);
load([PARAMS.intermediate_dir '\Naris_all_data_post.mat']);

%% Figure 1: individual gamma example (raw and heat map), averages across all rats
cfg_fig1 = [];
cfg_fig1.example = 63;
cfg_fig1.example2 = 42; 
cfg_fig1.example3 = 109;
cfg_fig1.session_name = 'R061_2014_09_26';
cfg_fig1.session_name2 = 'R054_2014_10_10'; %'R061_2014_09_26';
cfg_fig1.version = 1; % this means it will exclude the corrupted sessions for R054 and just use sessions from 2016_10_10 and 2016_10_13
AMPX_Naris_fig_1_example(cfg_fig1, all_data, all_data_post);

%% Figure 2: 2x2x2 task
load([PARAMS.intermediate_dir '\Naris_all_data_task.mat']);
AMPX_Naris_fig_2_task(all_data, all_data_task, 'save_fig', 'yes');

%% Figure 3: creates 4 x 2 figues of the phase differences for both across entire gamma events and within the triplets
AMPX_Naris_fig_3_phase([], all_data, all_data_post, 'save_fig', 'yes');

%% Figure 4: CSD for same event used in figure 1 C/D and average for each rat (50/80)
cfg_csd = [];
cfg_csd.example = cfg_fig1.example;
cfg_csd.session_name = cfg_fig1.session_name;
AMPX_Naris_fig_4_CSD(cfg_csd, all_data)

%% figure 6 power by distance
AMPX_Naris_pow_dist_corr(all_data)

%% get the spindle plot
cfg = [];
AMPX_Naris_PCA(cfg, all_data)

%% get all event stats
load([PARAMS.intermediate_dir '\Naris_all_data_task.mat']);
stats = Naris_get_all_event_stats(all_data,all_data_post,all_data_task);