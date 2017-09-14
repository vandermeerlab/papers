function stats = AMPX_Naris_get_event_stats(data_out, data_remap)
%  AMPX_Naris_get_event_stats: gets the statstics on all the events in a
%
%          Inputs:
%           - data_out: [struct] main data struct for either 'pre',
%           'task',or 'post'
%
%          Outputs:
%           - stats: [struct] contains the basic stats for events within
%           and across rats/sessions
%
% EC - 2016-12-03
%% setup some events structs
bands = {'lg', 'hg', 'lg_ran', 'hg_ran', 'spindles'};

%% collect all the events across sessions
for iband = 1:length(bands)
    events_sess.(bands{iband}).evt_length = [];
    events_sess.(bands{iband}).evt_cycles = [];
end

for iband = 1:length(bands)
    if isfield(data_out, bands{iband})
        for iEvt = 1:length(data_out.(bands{iband}).evts.tstart)
            events_sess.(bands{iband}).evt_length = [events_sess.(bands{iband}).evt_length, data_out.(bands{iband}).evts.tend(iEvt) - data_out.(bands{iband}).evts.tstart(iEvt)];
            if ~isempty(data_out.(bands{iband}).evts.usr)
                events_sess.(bands{iband}).evt_cycles = [ events_sess.(bands{iband}).evt_cycles, data_out.(bands{iband}).evts.usr.nCycles(iEvt)];
            end
        end
    else
        events_sess.(bands{iband}).evt_cycles =[];
        events_sess.(bands{iband}).evt_length = [];
    end
end
%% get the average length of the events

for iband = 1:length(bands)
    stats.(bands{iband}).event.avg_len = nanmean(events_sess.(bands{iband}).evt_length);
    stats.(bands{iband}).event.std_len = nanstd(events_sess.(bands{iband}).evt_length);
    stats.(bands{iband}).event.total = length(events_sess.(bands{iband}).evt_length);
    stats.(bands{iband}).event.rate = length(events_sess.(bands{iband}).evt_length)/(data_remap.tvec(end)/60); % divide by the number of minutes in the session.
    stats.(bands{iband}).event.avg_nCycle = nanmean(events_sess.(bands{iband}).evt_cycles);
    fprintf(['\n' bands{iband} 'mean length: ' num2str(stats.(bands{iband}).event.avg_len)...
        ' std: ' num2str(stats.(bands{iband}).event.std_len) ...
        ' nCycle: ' num2str(stats.(bands{iband}).event.avg_nCycle) '\n']);
    stats.(bands{iband}).all.length = events_sess.(bands{iband}).evt_length;
    if isfield(data_out, bands{iband}) ~=1 || ~isempty(data_out.(bands{iband}).evts.usr) 
        stats.(bands{iband}).all.nCycles = events_sess.(bands{iband}).evt_cycles;
    end
end


