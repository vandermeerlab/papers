function [csd_struct, csd_triplet_struct, cfg_fig4a, cfg_fig4_c_f] = AMPX_Naris_fig_4_CSD(cfg_in, all_data)
%% AMPX_Naris_fig_4_kCSD: uses the csd output from AMPX_CSD_example and
% AMPX_cycle_CSD to produce figure 4 for Carmichael, Gmaz, and van der
% Meer.  
%   
% 
% 
%          Inputs: 
%           - 
%           - 
%           - 
%          Outputs: 
%           - 
%           - 
%           -  
% 
% EC - 2016-07-09

%% Collect inputs/defaults
cfg_def.session_name = 'R061_2014_09_26';
cfg_def.example = 29;
cfg_def.chan_to_plot =  [1 10 19 28 37 46 55 64]; 
cfg = ProcessConfig(cfg_def, cfg_in);


%% generate a CSD for the same event used in figure 1c
% load('C:\temp\Naris_all_data_pre.mat');
[csd_struct, cfg_fig4a] = AMPX_CSD_example(cfg, all_data);

close all
%% generate a CSD for each triplet event and then average across all events to produce fig4_c-f

[csd_triplet_struct] = AMPX_cycle_CSD(cfg, all_data);