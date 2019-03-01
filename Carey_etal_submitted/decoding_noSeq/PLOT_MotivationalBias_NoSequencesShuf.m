function fh = PLOT_MotivationalBias_NoSequencesShuf(cfg_in,data)
% function fh = PLOT_Motivational_NoSequencesShuf(cfg_in,data)
%
% plots SWR content averaged over sessions as thresholded z-scores
% (MotivationalT data set)
%
% INPUTS:
%
% cfg_in: options that control display properties
%
% data: struct with .all, .pre, .task, .post fields containing data output
% from PLOT_DecSeqCombinedShuf.m
%
% CONFIGS:
%
% the `what` variable below can be set to {'all'} to plot & analyze all 
% data together, or to {'pre','task','post'} to plot epochs separately as 
% done in the paper figures

cfg_def = [];
cfg_def.showAllRatsText = 0;
cfg_def.colormode = 'inventory3';
cfg_def.colors = TmazeColors(cfg_def.colormode);
cfg_def.fs = 18; % font size
cfg_def.ylim = [1 2];
cfg_def.ylimtick = [0.5 0.5];
cfg_def.writeOutput = 0;
cfg_def.output_fn = 'temp';
cfg_def.what = {'all'}; % {'all'}, or {'pre', 'task', 'post'} to analyze pre-task, task and post-task epochs separately

cfg = ProcessConfig(cfg_def, cfg_in);

biasfun = @(d) (d(1)-d(2)); % computes bias measure as (food-water)

%%
what = cfg.what;
what_idx = {[1 2 7 8], [3 4 9 10], [5 6 11 12]};

% first plot data for all rats
rats = {'all'}; iRat = 1;
col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w [0 0 0]}; % {foodColor waterColor diffColor}

location = [3 2 1]; % where on y-axis to place the data (food restriction, water restriction, difference respectively)
recColor = [0.8 0.8 0.8]; %[0.5 0.5 0.5];
ylims = [0.5 max(location)+0.5];
ylab = {'FR', 'WR', 'diff'};

% main plots
for iW = 1:length(what) % loop over epochs
    
    subplot(4, 6, what_idx{iW});
    
    % arrange the data
    this_data_all{1} = data.(what{iW}).(rats{iRat}).food_fracL_evt;
    this_data_all{2} = data.(what{iW}).(rats{iRat}).water_fracL_evt;
    this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
    this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];
    
    % arrange the shuffled data
    this_shuf_mean = [data.(what{iW}).(rats{iRat}).foodShuf_fracL_evt data.(what{iW}).(rats{iRat}).waterShuf_fracL_evt data.(what{iW}).(rats{iRat}).diffShuf_fracL_evt + 0.5];
    this_shuf_sd = [data.(what{iW}).(rats{iRat}).foodShufsd_fracL_evt data.(what{iW}).(rats{iRat}).waterShufsd_fracL_evt data.(what{iW}).(rats{iRat}).diffShufsd_fracL_evt];
    
    % some stats
    m1 = nanmean(this_data_all{1}); s1 = nanstd(this_data_all{1}) ./ sqrt(4); % SEM
    m2 = nanmean(this_data_all{2}); s2 = nanstd(this_data_all{2}) ./ sqrt(4); % SEM
    p = ranksum(this_data_all{1}, this_data_all{2});
    fprintf('%s: food %.2f +/ %.2f, water %.2f +/- %.2f, p = %.2e (ranksum)\n', what{iW}, m1, s1, m2, s2, p);
    fprintf('%s: shuf SDs, food %.2f, water %.2f, diff %.2f\n', what{iW}, this_shuf_sd(1), this_shuf_sd(2), this_shuf_sd(3));
    
    % plot the shuffled data
    for iRec = 1:length(this_shuf_mean) % [x y w h] 
        
       rec_h(iRec) = rectangle('Position', [this_shuf_mean(iRec) - 0.5*this_shuf_sd(iRec) location(iRec) - 0.25 this_shuf_sd(iRec) 0.5]);
       set(rec_h(iRec), 'FaceColor', recColor, 'EdgeColor', 'none'); hold on;
        
       %rec_h(iRec) = errorbar(this_shuf_mean(iRec), location(iRec), 0.5*this_shuf_sd(iRec), 'horizontal'); hold on;
       %set(rec_h(iRec), 'LineWidth', 1, 'Color', recColor);
       
    end
        
    % plot the data
    for iBar = 1:length(this_data)
        mean_h(iBar) = plot([this_data(iBar) this_data(iBar)], [location(iBar)-0.33 location(iBar)+0.33], 'Color', col{iBar}, 'LineWidth', 3);
        hold on
        if iBar ~= 3
            plot(this_data_all{iBar}, location(iBar), '.', 'Color', col{iBar}, 'MarkerSize', 15);
        end
        
        % stats
        this_z = (this_data(iBar) - this_shuf_mean(iBar)) ./ this_shuf_sd(iBar);
        this_p = 2*normcdf(abs(this_z), 0, 1, 'upper');
        fprintf('%s: z %.2f, p %.2e\n', ylab{iBar}, this_z, this_p);
        
        DrawStars(cfg, this_p, [this_data(iBar) location(iBar)+0.37]);
             
    end
    set(gca, 'XTick', 0:0.25:1, 'FontSize', cfg.fs, 'LineWidth', 1, 'XLim', [0 1], 'XTickLabel', {'0 (W)', '', '0.5', '', '1 (F)'}, 'XDir', 'reverse');
    set(gca, 'YLim' , ylims, 'YTick', 1:3, 'TickDir', 'out', 'YTickLabel', ylab(location));
    title(sprintf('%s', what{iW})); box off
    
    set(gca,'Layer','top')
    
    plot([0.5 0.5], ylim, '--', 'Color', [0.5 0.5 0.5]); % vertical line
    
end % of whats for main plots

% insets for individual rats
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 subplot start positions for PRE

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];

    col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w [0 0 0]}; % {foodColor waterColor diffColor}

    for iW = 1:length(what)
    
        subplot(4,6,position(iW))
        
        % arrange the data
        this_data_all{1} = data.(what{iW}).(rats{iRat}).food_fracL_evt;
        this_data_all{2} = data.(what{iW}).(rats{iRat}).water_fracL_evt;
        this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
        this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];
        
        % arrange the shuffled data
        this_shuf_mean = [data.(what{iW}).(rats{iRat}).foodShuf_fracL_evt data.(what{iW}).(rats{iRat}).waterShuf_fracL_evt data.(what{iW}).(rats{iRat}).diffShuf_fracL_evt + 0.5];
        this_shuf_sd = [data.(what{iW}).(rats{iRat}).foodShufsd_fracL_evt data.(what{iW}).(rats{iRat}).waterShufsd_fracL_evt data.(what{iW}).(rats{iRat}).diffShufsd_fracL_evt];
               
        % plot the shuffled data
        for iRec = 1:length(this_shuf_mean) % [x y w h]
            
            rec_h(iRec) = rectangle('Position', [this_shuf_mean(iRec) - 0.5*this_shuf_sd(iRec) location(iRec) - 0.25 this_shuf_sd(iRec) 0.5]);
            set(rec_h(iRec), 'FaceColor', recColor, 'EdgeColor', 'none'); hold on;
            
            %rec_h(iRec) = errorbar(this_shuf_mean(iRec), location(iRec), 0.5*this_shuf_sd(iRec), 'horizontal'); hold on; box off;
            %set(rec_h(iRec), 'LineWidth', 0.5, 'Color', recColor);
            
        end
        
        % plot the data
        for iBar = 1:length(this_data)
            mean_h(iBar) = plot([this_data(iBar) this_data(iBar)], [location(iBar)-0.33 location(iBar)+0.33], 'Color', col{iBar}, 'LineWidth', 2);
            hold on
            if iBar ~= 3
                plot(this_data_all{iBar}, location(iBar), '.', 'Color', col{iBar}, 'MarkerSize', 10);
            end
                       
        end
        set(gca, 'XTick', 0:0.25:1, 'FontSize', 8, 'LineWidth', 1, 'XLim', [0 1], 'XTickLabel', {'0 (W)', '', '0.5', '', '1 (F)'}, 'XDir', 'reverse');
        if iRat == 2
            set(gca, 'XLim', [-0.5 1.5]);
        end
        
        set(gca, 'YLim' , ylims, 'YTick', 1:3, 'TickDir', 'out', 'YTickLabel', ylab(location));
        %title(sprintf('%s', what{iW})); box off
        set(gca,'Layer','top')
        
        plot([0.5 0.5], ylim, '--', 'Color', [0.5 0.5 0.5]);
        
    end % of ratwhats
    
end % of rats


% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end




