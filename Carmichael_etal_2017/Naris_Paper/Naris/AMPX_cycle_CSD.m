function [CSD] = AMPX_cycle_CSD(cfg_in, all_data)
%% AMPX_cycle_CSD: computes the current source density using the 5pt method by E.W.Schomburg, (March 2014) for
%  each three cycle event in all_cycles.
%
%
%          Inputs:
%           - cfg_in: [struct] contains paramters
%           - all_cycles: [struct] output from AMPX_extract_max_cycle
%          Outputs:
%           - cycle_csd: [struct] contains the
%
% EC - 2016-05-19

global PARAMS;

%% parameters
cfg_def.interp_factor = 0.1; % new spacing for interpolation
cfg_def.elec_sep_mm = sqrt(0.2^2+0.2^2); % probe spacing for a NeuroNexus A8x8 probe wii .2mm regular spacing between sites
cfg_def.chan_to_plot = diag(all_data.R054_2014_10_10.lg.cycles.ExpKeys.Probe_layout);
cfg_def.c_axis_for_scale = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);

%% cycle all the events for each rat
session_list = fieldnames(all_data);
bands = {'lg', 'hg'};
peaks = {'prev', 'center', 'next'};

for iSess =1:length(session_list)
    for iband = 1:2
        all_cycles = all_data.(session_list{iSess}).(bands{iband}).cycles;
        for ievt =length(all_cycles.peaks.center.data):-1:1;
            %             for iPeak = 1:length(peaks)
            if isempty(all_cycles.data{ievt})
                continue
            end
            % interp missing channels
            [X, Y] = meshgrid(1:8);
            temp_data = all_cycles.data{ievt};
            rm_chan = all_cycles.ExpKeys.BadChannels;
            for iChan = rm_chan
                temp_data(iChan,:) = NaN; % ensure all the bad channels are NaN'd out before using grid data to interp over them correctly.
            end
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
            
            % instead of trying to interpolate or extrapolate this value, we will instead only use the 7 remaining diagonal channels.
            if isnan(temp(1)) || isnan(temp(64))
                xcpots = temp_data(cfg.chan_to_plot(~isnan(diag(reshape(temp,8,8)))),:); %extract the channels on the dorsomedial to ventrolateral diagonal of the probe.
                
                %  z_grid = interp1(diag(X), diag(temp),1:8,'spline'); % DON"T DO THIS!!!! extrapolate using spline method for NaN on corners of probe.
            else
                % just to make sure nothing gets wierd.
                xcpots = temp_data(cfg.chan_to_plot,:); %extract the channels on the dorsomedial to ventrolateral diagonal of the probe.
                
            end
            
            if cfg.c_axis_for_scale ==1
                % scale the CSD axis to be relative to a 180 phase inversion of the highest
                % amplitude channel.
                xcpots_for_scale = xcpots;
                xcpots_for_scale(end-2,:) = -xcpots_for_scale(end-1,:);
                [csd_for_scale, ~] = csd5pt(xcpots_for_scale, cfg.elec_sep_mm);
                min_max_inverse = [min(min(csd_for_scale)) max(max(csd_for_scale))];
                cfg.min_max_inverse = min_max_inverse;
            end
            % CSD
            [csd, CSDelecinds] = csd5pt(xcpots, cfg.elec_sep_mm);
            %                 if iSess >=6 && iSess < 9
            figure(iband)
            %                 subplot(15,15,ievt)
            %                 plot(all_cycles.tvec{ievt}, xcpots)
            %                 close(gcf)
            imagesc(all_cycles.tvec{ievt},CSDelecinds,csd);
            caxis( min_max_inverse)  % set the colorbar to be relative to a full 180 degree phase inversion.
            %                 end
            CSD.(session_list{iSess}).(bands{iband}).all_csd_2d(:,:,ievt) = csd;
            CSD.(session_list{iSess}).(bands{iband}).min_max(:,ievt) = min_max_inverse;
            CSD.(session_list{iSess}).(bands{iband}).CSDelecinds = CSDelecinds;
        end
        %         end
    end
end
%     if CSD.(cfg.band).failures(1, ievt) ==0
%         mkdir([datestr(date, 'yyyy-mm-dd') '/CSD/'])
%         saveas(h_cycle, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_' num2str(ievt) '_CSD_cycle.fig'])
%         saveas(h_cycle, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_' num2str(ievt) '_CSD_cycle.png'])
%     end
close all

%% remove any event matrices built from a missing cycle.  This can happen if the cycle was excluded during the cycle extraction, likely if it overlapped with another or was at the start or end of a session
for iSess =1:length(session_list)
    for iband = 1:2
        for ievt =length(all_data.(session_list{iSess}).(bands{iband}).cycles.data):-1:1;
            if isempty(all_data.(session_list{iSess}).(bands{iband}).cycles.data{ievt})
                CSD.(session_list{iSess}).(bands{iband}).all_csd_2d(:,:,ievt) = [];
            end
        end
    end
end

%% collect all the cycles per subject/band/cycle for averaging
for iband = 1:2
    for iPeak = 1:3
        CSD.R061.(bands{iband}).csd = nanmean(cat(3,CSD.R061_2014_09_26.(bands{iband}).all_csd_2d,CSD.R061_2014_09_27.(bands{iband}).all_csd_2d,CSD.R061_2014_09_28.(bands{iband}).all_csd_2d),3) ;
        CSD.R054.(bands{iband}).csd = nanmean(cat(3,CSD.R054_2014_10_10.(bands{iband}).all_csd_2d,CSD.R054_2014_10_13.(bands{iband}).all_csd_2d),3) ;
        CSD.R049.(bands{iband}).csd = nanmean(cat(3,CSD.R049_2014_02_07.(bands{iband}).all_csd_2d,CSD.R049_2014_02_08.(bands{iband}).all_csd_2d,CSD.R049_2014_02_10.(bands{iband}).all_csd_2d),3) ;
        CSD.R045.(bands{iband}).csd = nanmean(cat(3,CSD.R045_2014_04_16.(bands{iband}).all_csd_2d,CSD.R045_2014_04_17.(bands{iband}).all_csd_2d,CSD.R045_2014_04_15.(bands{iband}).all_csd_2d),3) ;
    end
end

%% collect all the CSD min_max_inverse values
for iband = 1:2
    for iPeak = 1:3
        CSD.R061.(bands{iband}).min_max = [min([CSD.R061_2014_09_26.(bands{iband}).min_max(1,:),CSD.R061_2014_09_27.(bands{iband}).min_max(1,:),CSD.R061_2014_09_28.(bands{iband}).min_max(1,:)]), max([CSD.R061_2014_09_26.(bands{iband}).min_max(2,:),CSD.R061_2014_09_27.(bands{iband}).min_max(2,:),CSD.R061_2014_09_28.(bands{iband}).min_max(2,:)])] ;
        CSD.R054.(bands{iband}).min_max = [min([CSD.R054_2014_10_10.(bands{iband}).min_max(1,:),CSD.R054_2014_10_13.(bands{iband}).min_max(1,:)]), max([CSD.R054_2014_10_10.(bands{iband}).min_max(2,:),CSD.R054_2014_10_13.(bands{iband}).min_max(2,:)])] ;
        CSD.R049.(bands{iband}).min_max = [min([CSD.R049_2014_02_07.(bands{iband}).min_max(1,:),CSD.R049_2014_02_08.(bands{iband}).min_max(1,:),CSD.R049_2014_02_10.(bands{iband}).min_max(1,:)]), max([CSD.R049_2014_02_07.(bands{iband}).min_max(2,:),CSD.R049_2014_02_08.(bands{iband}).min_max(2,:),CSD.R049_2014_02_10.(bands{iband}).min_max(2,:)])] ;
        CSD.R045.(bands{iband}).min_max = [min([CSD.R045_2014_04_16.(bands{iband}).min_max(1,:),CSD.R045_2014_04_17.(bands{iband}).min_max(1,:),CSD.R045_2014_04_15.(bands{iband}).min_max(1,:)]), max([CSD.R045_2014_04_16.(bands{iband}).min_max(2,:),CSD.R045_2014_04_17.(bands{iband}).min_max(2,:),CSD.R045_2014_04_15.(bands{iband}).min_max(2,:)])] ;
    end
end
%% plot the average on a new figure and save the data.
% tvec = 1:length(all_data.(session_list{iSess}).(bands{iband}).cycles.tvec{10}); % since all of the cycles were of the same length the tvec length will be the same number of samples.
for iband = 1:2
    % h_avg = figure(10);
    %     set(gcf, 'position', [677 526 736 215])
    %     hold on
    h1 = figure(1);
    imagesc(1:length(csd),CSD.R061_2014_09_26.(bands{iband}).CSDelecinds, CSD.R061.(bands{iband}).csd)
    p = get(h1, 'pos');
    p(3) = p(3) + 0.05;
    set(h1, 'pos', p);
    caxis(CSD.R061.(bands{iband}).min_max)
    xlim([1 length(csd)])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [CSD.R061_2014_09_26.(bands{iband}).CSDelecinds(1), CSD.R061_2014_09_26.(bands{iband}).CSDelecinds(end)], 'yticklabel', {}, 'fontname', 'helvetica')
%     xlabel('triplet phase'); ylabel('electrode number');
    SetFigure([], gcf)
    %         set(gca, 'xticklabel', {'-3pi', '-2pi', '-pi', '0', 'p', '2p', '3p'})
    
    % text(30, 5, 'prev')
    %     axis off
    print(cat(2,PARAMS.figure_dir,'\CSD_R061_',bands{iband},'_cycles'), '-depsc');
    saveas(gcf,cat(2,PARAMS.figure_dir,'\CSD_R061_',datestr(date, 'yyyy-mm-dd'),bands{iband},'_avg_CSD_cycle.fig'));
    
    h2 = figure(2);
    imagesc(1:length(csd),CSD.R054_2014_10_10.(bands{iband}).CSDelecinds,CSD.R054.(bands{iband}).csd)
    % text(30, 5, 'center')
    caxis(CSD.R054.(bands{iband}).min_max)
    %     axis off
    p = get(h2, 'pos');
    p(3) = p(3) + 0.05;
    set(h2, 'pos', p);
    xlim([1 length(csd)])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [CSD.R054_2014_10_10.(bands{iband}).CSDelecinds(1), CSD.R054_2014_10_10.(bands{iband}).CSDelecinds(end)], 'yticklabel', {}, 'fontname', 'helvetica')
%     xlabel('triplet phase'); ylabel('electrode number');
    SetFigure([], gcf)
    print(cat(2,PARAMS.figure_dir,'\CSD_R054_',bands{iband},'_cycles'), '-depsc');
    saveas(gcf,cat(2,PARAMS.figure_dir,'\CSD_R054_',datestr(date, 'yyyy-mm-dd'),bands{iband},'_avg_CSD_cycle.fig'));
    
    
    
    h3 = figure(3);
    imagesc(1:length(csd),CSD.R049_2014_02_07.(bands{iband}).CSDelecinds,CSD.R049.(bands{iband}).csd)
    caxis(CSD.R049.(bands{iband}).min_max)
    % text(30, 5, 'next')
    %     axis off
    p = get(h3, 'pos');
    p(3) = p(3) + 0.05;
    set(h3, 'pos', p);
    xlim([1 length(csd)])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [CSD.R049_2014_02_07.(bands{iband}).CSDelecinds(1), CSD.R049_2014_02_07.(bands{iband}).CSDelecinds(end)], 'yticklabel', {}, 'fontname', 'helvetica')
%     xlabel('triplet phase'); ylabel('electrode number');
    SetFigure([], gcf)
    SetFigure([], gcf)
    set(gca, 'xtick', [])
    print(cat(2,PARAMS.figure_dir,'\CSD_R049_',bands{iband},'_cycles'), '-depsc');
    saveas(gcf,cat(2,PARAMS.figure_dir,'\CSD_R049_',datestr(date, 'yyyy-mm-dd'),bands{iband},'_avg_CSD_cycle.fig'));
    
    h4 = figure(4);
    imagesc(1:length(csd),CSD.R045_2014_04_16.(bands{iband}).CSDelecinds,CSD.R045.(bands{iband}).csd)
    % text(30, 5, 'next')
    %     axis off
    p = get(h4, 'pos');
    p(4) = p(4) + 0.05;
    set(h4, 'pos', p);
    caxis(CSD.R045.(bands{iband}).min_max)
    xlim([1 length(csd)])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [CSD.R045_2014_04_16.(bands{iband}).CSDelecinds(1), CSD.R045_2014_04_16.(bands{iband}).CSDelecinds(end)], 'yticklabel', {}, 'fontname', 'helvetica')
%     xlabel('triplet phase'); ylabel('electrode number');
    SetFigure([], gcf)
    SetFigure([], gcf)
    axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
    set(axesHandles, 'fontsize', 18, 'fontname', 'helvetica')
    print(cat(2,PARAMS.figure_dir,'\CSD_R045_',bands{iband},'_cycles'), '-depsc');
    saveas(gcf,cat(2,PARAMS.figure_dir,'\CSD_R045_',datestr(date, 'yyyy-mm-dd'),bands{iband},'_avg_CSD_cycle.fig'));
    
    close all
end

%% have a look at all the CSDs

% for iband = 1:2
%     figure(iband)
%     if length(CSD.R0
%     evts = length(
%     for ievt =
%     subplot(8,8,)