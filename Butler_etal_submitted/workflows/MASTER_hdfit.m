%%
clear all; clear global; pack
%addpath('D:\My_Documents\GitHub\hdmodelfit\shared');
addpath('C:\Users\mvdm\Documents\GitHub\hdmodelfit\shared'); % isidro

%datapath = 'D:\My_Documents\Dropbox\projects\HDfit\data';
datapath = 'C:\Users\mvdm\Dropbox\projects\HDfit\data';
cd(datapath);

%% test
global param; % to access internals of crossval();
global param_count;
global param_hist;
param_count = 0;

%out = hdfit_crossval_func([]);

%% multirun

clear ALL_out;

sessions_to_run = 3;%1:12;
targets_to_run = {'laser'};
%targets_to_run = {'std'};
%filters_to_run = {'smooth','kalmanwrapped'};
filters_to_run = {'kalmanwrapped'}%,'kalman','kalmanwrapped'};

for iS = length(sessions_to_run):-1:1
    
    fprintf('\n\n*** SESSION %d/%d ***\n\n',iS,length(sessions_to_run));
    
    for iT = 1:length(targets_to_run)
        
        for iF = 1:length(filters_to_run)
            
            fprintf('\n\n--> %s\n\n',filters_to_run{iF});

            this_cfg = [];
            this_cfg.session = sessions_to_run(iS);
            this_cfg.target_session = targets_to_run{iT};
            this_cfg.mode = filters_to_run{iF};
            this_cfg.datapath = datapath;
            this_cfg.gainbin_centers = 0.9:0.025:1.1;
            this_cfg.driftbin_centers = -1.5:0.25:1.5;
            this_cfg.simparams = [1 1 0.5]; % leave empty to use real data, spec [g_l g_r d] for simulated data
            this_cfg.models = 4;
            
            this_out = hdfit_crossval_func(this_cfg);
            
            ALL_out(iS).(targets_to_run{iT}).(filters_to_run{iF}) = this_out;
        end
    end
    
    save ALL_out;
end

   