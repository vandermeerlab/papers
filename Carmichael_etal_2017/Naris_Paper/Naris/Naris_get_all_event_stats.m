function all_stats = Naris_get_all_event_stats(all_data, all_data_post, all_data_task)

global PARAMS;

bands = {'lg', 'hg', 'lg_ran', 'hg_ran', 'spindles'};
all_stats = [];
for iband = 1:length(bands)
    all_stats.(bands{iband}).length = [];
    all_stats.(bands{iband}).nCycles = [];
    all_stats.(bands{iband}).rate = [];
    all_stats.(bands{iband}).total = [];
end

% across all 'pre' sessions
sess_list = fieldnames(all_data);
type = 'pre';
for iSess = 1:length(sess_list)
    cd(cat(2,PARAMS.data_dir,'\',sess_list{iSess}(1:4),'\',strrep(sess_list{iSess},'_', '-')));
    if strcmp(type, 'pre')
        fname = [strrep(sess_list{iSess},'_', '-') '-pre.dat'];
    elseif strcmp(type, 'post')
        fname = [strrep(sess_list{iSess},'_', '-') '-post.dat'];
    else
        fname = [strrep(sess_list{iSess},'_', '-') '.dat'];
    end
    data_in = AMPX_loadData(fname, 1, 10); % only used to get the length of the recording to determine the rate of gamma events
    all_stats.(type).(sess_list{iSess}) = AMPX_Naris_get_event_stats(all_data.(sess_list{iSess}), data_in);
    
    for iband = 1:length(bands)
        all_stats.(bands{iband}).length = [all_stats.(bands{iband}).length, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.length];
        all_stats.(bands{iband}).rate = [all_stats.(bands{iband}).rate, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.rate];
        all_stats.(bands{iband}).total = [all_stats.(bands{iband}).total, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.total];
        if isfield(all_stats.(type).(sess_list{iSess}).(bands{iband}).all, 'nCycles')
            all_stats.(bands{iband}).nCycles = [all_stats.(bands{iband}).nCycles, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.nCycles];
        end
    end
end

%% across all 'post' sessions
sess_list = fieldnames(all_data_post);
type = 'post';
for iSess = 1:length(sess_list)
    cd(cat(2,PARAMS.data_dir,'\',sess_list{iSess}(1:4),'\',strrep(sess_list{iSess},'_', '-')));
    
    fname = [strrep(sess_list{iSess},'_', '-') '-post.dat'];
    data_in = AMPX_loadData(fname, 1, 10); % only used to get the length of the recording to determine the rate of gamma events
    all_stats.(type).(sess_list{iSess}) = AMPX_Naris_get_event_stats(all_data_post.(sess_list{iSess}), data_in);
    
    for iband = 1:length(bands)
        all_stats.(bands{iband}).length = [all_stats.(bands{iband}).length, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.length];
        all_stats.(bands{iband}).rate = [all_stats.(bands{iband}).rate, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.rate];
        all_stats.(bands{iband}).total = [all_stats.(bands{iband}).total, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.total];
        if isfield(all_stats.(type).(sess_list{iSess}).(bands{iband}).all, 'nCycles')
            all_stats.(bands{iband}).nCycles = [all_stats.(bands{iband}).nCycles, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.nCycles];
        end
    end
end
%% task
%% across all 'task' sessions
sess_list = fieldnames(all_data_task);
type = 'task';
for iSess = 1:length(sess_list)
    cd(cat(2,PARAMS.data_dir,'\',sess_list{iSess}(1:4),'\',strrep(sess_list{iSess},'_', '-')));
    
    fname = [strrep(sess_list{iSess},'_', '-') '.dat'];
    data_in = AMPX_loadData(fname, 1, 10); % only used to get the length of the recording to determine the rate of gamma events
    all_stats.(type).(sess_list{iSess}) = AMPX_Naris_get_event_stats(all_data_post.(sess_list{iSess}), data_in);
    
    for iband = 1:length(bands)
        all_stats.(bands{iband}).length = [all_stats.(bands{iband}).length, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.length];
        all_stats.(bands{iband}).rate = [all_stats.(bands{iband}).rate, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.rate];
        all_stats.(bands{iband}).total = [all_stats.(bands{iband}).total, all_stats.(type).(sess_list{iSess}).(bands{iband}).event.total];
        if isfield(all_stats.(type).(sess_list{iSess}).(bands{iband}).all, 'nCycles')
            all_stats.(bands{iband}).nCycles = [all_stats.(bands{iband}).nCycles, all_stats.(type).(sess_list{iSess}).(bands{iband}).all.nCycles];
        end
    end
end
%% print all the stats
fprintf('\n\n\nGamma Event Stats\n')
fileID = fopen(cat(2,PARAMS.stats_dir,'\Naris_stats_events.txt'),'w');
fprintf(fileID,'Gamma Event Stats\n')
fprintf(fileID, ['_________________________________________\n'])
fprintf(fileID, [date '\n'])

for iband = 1:length(bands)
    all_stats.(bands{iband}).avg_len = nanmean(all_stats.(bands{iband}).length);
    all_stats.(bands{iband}).std_len = nanstd(all_stats.(bands{iband}).length);
    all_stats.(bands{iband}).total = length(all_stats.(bands{iband}).length);
    all_stats.(bands{iband}).avg_rate = nanmean(all_stats.(bands{iband}).rate); % divide by the number of minutes in the session.
    all_stats.(bands{iband}).avg_nCycle = nanmean(all_stats.(bands{iband}).nCycles);
    fprintf(fileID, [bands{iband}  ':       mean length: ' num2str(all_stats.(bands{iband}).avg_len*1000)...
        ' std: ' num2str(all_stats.(bands{iband}).std_len*1000) ...
        ' nCycle: ' num2str(all_stats.(bands{iband}).avg_nCycle)...
        ' rate: ' num2str(all_stats.(bands{iband}).avg_rate)...
        ' total: ' num2str(all_stats.(bands{iband}).total) '\n']);
end

fclose(fileID);

