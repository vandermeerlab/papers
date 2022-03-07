%% list of sessions to analyze
% Session_list = {'R054-2014-10-11', 'R054-2014-10-12', 'R054-2014-10-13', 'R054-2014-10-14'};
Session_list = {'R054-2014-10-10', 'R054-2014-10-13', ...
    'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',... % 'R049-2014-02-09',...
    'R061-2014-09-26', 'R061-2014-09-27', 'R061-2014-09-28',...
    'R045-2014-04-16', 'R045-2014-04-15', 'R045-2014-04-17'};

type = 'pre';

% load('C:\temp\Naris_all_data_pre_Paper_spin_nov.mat')
%% loop over each session to get: events, power, event_phase, middle cycles, cycle_phase

for iSess =1:length(Session_list)
    disp(['Running Session: ' (strrep(Session_list{iSess},'-', '_'))])
    if strcmp(type, 'pre') == 1
        all_data.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}],all_data.(strrep(Session_list{iSess},'-', '_')),'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'post') == 1
        all_data_post.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}],all_data_post.(strrep(Session_list{iSess},'-', '_')), 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'task') ==1
        all_data_task.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}],all_data_task.(strrep(Session_list{iSess},'-', '_')), 'session_type', type, 'plane_plot', 'yes');
    end
end

%% save the all_data struct
if strcmp(type, 'pre') == 1
    save('C:\temp\Naris_all_data_pre_Paper_spin.mat', 'all_data', '-v7.3')
elseif strcmp(type, 'task') ==1
    save('C:\temp\Naris_all_data_task.mat', 'all_data_task', '-v7.3')
elseif strcmp(type, 'post') ==1
    save('C:\temp\Naris_all_data_post.mat', 'all_data_post', '-v7.3')
end

clear all; close all
%% list of sessions to analyze
% Session_list = {'R054-2014-10-11', 'R054-2014-10-12', 'R054-2014-10-13', 'R054-2014-10-14'};
Session_list = {'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',... % 'R049-2014-02-09',...
    'R045-2014-04-16', 'R045-2014-04-17', 'R045-2014-04-15'};

type = 'post';
% load('C:\temp\Naris_all_data_post_nov.mat')
%% loop over each session to get: events, power, event_phase, middle cycles, cycle_phase

for iSess =1:length(Session_list)
    disp(['Running Session: ' (strrep(Session_list{iSess},'-', '_'))])
    if strcmp(type, 'pre') == 1
        all_data.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'post') == 1
        all_data_post.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'task') ==1
        all_data_task.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    end
end

%% save the all_data struct
if strcmp(type, 'pre') == 1
    save('C:\temp\Naris_all_data_pre_Paper_spin.mat', 'all_data', '-v7.3')
elseif strcmp(type, 'task') ==1
    save('C:\temp\Naris_all_data_task.mat', 'all_data_task', '-v7.3')
elseif strcmp(type, 'post') ==1
    save('C:\temp\Naris_all_data_post.mat', 'all_data_post', '-v7.3')
end

clear all; close all
%% get the task sessions
Session_list = {'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',...
    'R045-2014-04-16', 'R045-2014-04-17', 'R045-2014-04-15'};
type = 'task';

%% loop over each session to get: events, power, event_phase, middle cycles, cycle_phase

for iSess =1:length(Session_list)
    disp(['Running Session: ' (strrep(Session_list{iSess},'-', '_'))])
    if strcmp(type, 'pre') == 1
        all_data.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'post') == 1
        all_data_post.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'task') ==1
        all_data_task.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    end
end

%% save the all_data struct
if strcmp(type, 'pre') == 1
    save('C:\temp\Naris_all_data_pre_Paper_spin.mat', 'all_data', '-v7.3')
elseif strcmp(type, 'task') ==1
    save('C:\temp\Naris_all_data_task.mat', 'all_data_task', '-v7.3')
elseif strcmp(type, 'post') ==1
    save('C:\temp\Naris_all_data_post.mat', 'all_data_post', '-v7.3')
end

clear all; close all

addpath(genpath('C:\Users\mvdmlab\Dropbox\Matlab\BuzCSD'))
AMPX_Naris_generate_figures()