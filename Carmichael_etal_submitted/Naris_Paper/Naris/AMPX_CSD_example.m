function [csd_struct, cfg_out] = AMPX_CSD_example(cfg_in, all_data)
%% AMPX_CSD_example: takes the data for the specified example (should be 
% the same one used in figure 1 (AMPX_Naris_fig1_example, as cfg.example 
% for the session given in cfg._session_name). 
%

% Computes the 5ptCSD for a single event using the methods from
global PARAMS
%% Preamble
cfg_def.buffer = 0.5;
cfg_def.f = [40 55]; % default filter range in the low gamma band.
cfg_def.c_axis_for_scale = 1; % sets the colorbar axis to reflect the values relative to a full 180 degree phase inversion between two pairs. 
cfg_def.save = 1; % save the figures. 
cfg_def.save_dir = PARAMS.figure_dir; % where they are saved
cfg = ProcessConfig2(cfg_in, cfg_def);

ExpKeys = all_data.(cfg.session_name).lg.cycles.ExpKeys;
evts = all_data.(cfg.session_name).lg.evts;

addpath(genpath(PARAMS.CSD_dir))

%% extract the data for the example gamma event
fname = strrep(cfg.session_name, '_', '-');
cd([PARAMS.data_dir '\' fname(1:4) '\' fname(1:15) ])
fname = strrep(cfg.session_name, '-', '_');
[data, ~] = AMPX_Naris_preprocess([],fname,'pre');
data_remap_AMPX = AMPX_remap_sites(data, ExpKeys);
ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys, data_remap_AMPX);
cfg.chan_to_plot = diag(ExpKeys_remap.Probe_layout);

data_tsd = AMPX2tsd(data_remap_AMPX);
clear data
clear data_remap_AMPX

Remove_ft()
cfg_filter = [];
cfg_filter.f = cfg.f;
cfg_filter.type = 'cheby1';
cfg_filter.order = 5;
event_data = FilterLFP(cfg_filter, data_tsd);

event_data_filt = restrict(event_data, evts.tstart(cfg.example)-cfg.buffer/2, evts.tend(cfg.example)+cfg.buffer/2);
event_data_not_filt = restrict(data_tsd, evts.tstart(cfg.example)-cfg.buffer/2, evts.tend(cfg.example)+cfg.buffer/2);


rm_chan = ExpKeys.BadChannels;
for iChan = rm_chan
    event_data_filt.data(iChan,:) = NaN;
end
% convert to simpler variables
filt_data = event_data_filt.data;
filt_tvec = event_data_filt.tvec;
%% interpolate over missing channels using griddata.
[X, Y] = meshgrid(1:8);
temp_data = filt_data;
for ii = length(temp_data):-1:1
    temp = reshape(temp_data(:,ii),8,8);
    x_known = X(~isnan(temp));  % find NaN values and remove them from the X / Y inputs
    y_known = Y(~isnan(temp));
    z_grid = griddata(x_known,y_known,temp(~isnan(temp)),X(isnan(temp)),Y(isnan(temp))); % query for the NaN values in the data to be interpolated.
    temp(isnan(temp)) = z_grid; % put the interpolated points in for the corresponding NaN values.
    %     temp = inpaintn(temp); % this is another method but it was more difficult to explain than grid data.
    temp = reshape(temp,1,64);
    temp_data(:,ii) = temp;
end


%% Compute the CSD using the 5pt method

%create sites by time matrix (input to csd5pt.m)
xcpots = temp_data(cfg.chan_to_plot,:); %extract the channels on the dorsomedial to ventrolateral diagonal of the probe.  
elec_sep_mm = sqrt(0.2^2+0.2^2); % probe spacing for a NeuroNexus A8x8 probe wii .2mm regular spacing between sites

% CSD
[csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);
label = num2cell(1:8);
% legend('show')


%% plot the filtered trace of the example event along with the output from the CSD.
% If you want to set the scale to be relative to a 180 degree phase inversion use cfg.c_axis_for_scale = 1;
exfig = figure(100);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [36 70 1824 936])
c_ord = linspecer(8); % nice colours from linspecer.m (https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap/content/linspecer.m)
subplot(2,1,1)
hold on
for iChan = 1:8
    plot(filt_tvec',temp_data(cfg.chan_to_plot(iChan),:), 'color', c_ord(iChan, :)); % plot the filtered trace
end
xlim([filt_tvec(1) filt_tvec(end)]);
% legend(num2str(cfg.chan_to_plot))
ylim([-250 250])
y_lim = ylim; % needs to be done or the next line doesn't work 
rectangle('position',[evts.tstart(cfg.example), y_lim(1)+5, evts.tend(cfg.example) - evts.tstart(cfg.example), 10], 'facecolor', [.6 .6 .6], 'edgecolor', [.6 .6 .6])
set(gca, 'ytick', [-250 0 250])
cfg_set_fig.resize = 0;
ylabel('amplitude')
SetFigure(cfg_set_fig, gcf) % sets the figure properties to mvdmlab default. 

subplot(212) % plot the CSD
imagesc(filt_tvec,CSDelecinds,csd);
pos = get(gca, 'position');
colorbar('location', 'eastoutside'); % note units now in uV/mm^2... but numbers here seem very small?
set(gca,'position',pos) % sets figure size
set(gca, 'ytick',[2.5 3 4 5 6 6.5], 'yticklabel', {'DM','3', '4', '5', '6', 'VL'})
ylabel('diagonal electrode')
xlabel('time (ms)')

if cfg.c_axis_for_scale
    % scale the CSD axis to be relative to a 180 phase inversion of the highest
    % amplitude channel.
    xcpots_for_scale = xcpots;
    xcpots_for_scale(7,:) = -xcpots_for_scale(8,:);
        xcpots_for_scale(2,:) = -xcpots_for_scale(8,:);
    [csd_for_scale, CSDelecinds_for_scale] = csd5pt(xcpots_for_scale, elec_sep_mm);
    min_max_inverse = [min(min(csd_for_scale)) max(max(csd_for_scale))];
    caxis( min_max_inverse)  % set the colorbar to be relative to a full 180 degree phase inversion.
    cfg.min_max_inverse = min_max_inverse;
end
% format the figure
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
set(axesHandles, 'xtick', [filt_tvec(1) filt_tvec(end)], 'xticklabel', [0 500]);

cfg_set_fig.resize = 0;
SetFigure(cfg_set_fig, gcf) % sets the figure properties to mvdmlab default. 


%% save the figure
if cfg.save == 1
    saveas(gcf, [PARAMS.figure_dir '\Fig4_A'], 'epsc')

    saveas(gcf, [PARAMS.figure_dir '\Fig4_A'], 'fig')
end
% output the CSD and config structs
csd_struct.CSDelecinds = CSDelecinds; 
csd_struct.csd = csd; 
csd_struct.xcpot = xcpots;
csd_struct.tvec = filt_tvec;

cfg_out = cfg;

% %% triplet example
% % low
% figure(11)
% hold on
% subplot(211)
%  rm_chan = ExpKeys.BadChannels;
%  temp_data = all_data.(cfg.session_name).lg.cycles.data{cfg.example};
%  for iChan = rm_chan
%      temp_data(iChan,:) = NaN; % ensure all the bad channels are NaN'd out before using grid data to interp over them correctly.
%  end
%  for ii = length(temp_data):-1:1
%      temp = reshape(temp_data(:,ii),8,8);
%      x_known = X(~isnan(temp));  % find NaN values and remove them from the X / Y inputs
%      y_known = Y(~isnan(temp));
%      z_grid = griddata(x_known,y_known,temp(~isnan(temp)),X(isnan(temp)),Y(isnan(temp))); % query for the NaN values in the data to be interpolated.
%      temp(isnan(temp)) = z_grid; % put the interpolated points in for the corresponding NaN values.
%      %     temp = inpaintn(temp); % this is another method but it was more difficult to explain than grid data.
%      temp = reshape(temp,1,64);
%      temp_data(:,ii) = temp;
%  end
% for iChan = 1:8
%     hold on
%     plot(temp_data(cfg.chan_to_plot(iChan),:), 'color', c_ord(iChan, :), 'linewidth', 2); % plot the filtered trace
% end
% xlim([1 length(all_data.(cfg.session_name).lg.cycles.data{cfg.example})])
% set(gca, 'ytick', [], 'xtick', [1:(round(length(all_data.(cfg.session_name).lg.cycles.data{cfg.example})/6)):length(all_data.(cfg.session_name).lg.cycles.data{cfg.example}),length(all_data.(cfg.session_name).lg.cycles.data{cfg.example}) ])
% SetFigure([], gcf)
% set(gca, 'xticklabels', {'-3p' '-2p', '-p', '0', 'p', '2p' '3p'}, 'fontname', 'symbol')
% subplot(212)
% xcpots = temp_data(cfg.chan_to_plot,:); %extract the channels on the dorsomedial to ventrolateral diagonal of the probe.  
% elec_sep_mm = sqrt(0.2^2+0.2^2); % probe spacing for a NeuroNexus A8x8 probe wii .2mm regular spacing between sites
% % CSD
% [csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);
% imagesc(1:length(csd), CSDelecinds, csd)
% if cfg.save == 1
% %     saveas(gcf, [cfg.save_dir 'Fig4_low_trip'], 'epsc')
% 
%     saveas(gcf, [cfg.save_dir 'Fig4_A'], 'fig')
% end
%%
close()
%     cfg.example = randi([1 length(all_data.(cfg.session_name).hg.cycles.data)],1,1);
% figure(22)
% hold on
% for iChan = 1:8
%     plot(all_data.(cfg.session_name).hg.cycles.data{cfg.example}(cfg.chan_to_plot(iChan),:)', 'color', c_ord(iChan, :), 'linewidth',2); % plot the filtered trace
% end
% xlim([1 length(all_data.(cfg.session_name).hg.cycles.data{cfg.example})])
% set(gca, 'ytick', [], 'xtick', [1:(round(length(all_data.(cfg.session_name).hg.cycles.data{cfg.example})/6)):length(all_data.(cfg.session_name).hg.cycles.data{cfg.example}),length(all_data.(cfg.session_name).lg.cycles.data{cfg.example}) ])
% SetFigure([], gcf)
% set(gca, 'xticklabels', {'-3p' '-2p', '-p', '0', 'p', '2p' '3p'}, 'fontname', 'symbol')
% 
% if cfg.save == 1
%     saveas(gcf, [cfg.save_dir 'Fig4_A_high_trip'], 'epsc')
% 
%     saveas(gcf, [cfg.save_dir 'Fig4_A_high_trip'], 'fig')
% end
