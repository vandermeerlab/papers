%% Generate behavioral data figure (Figure 1B)

% first, run Behavior_GenData
% then invert water choices, because stored as correct/incorrect, not
% left/right
all.sessionChoice(all.sessionType == 2, :) = 1 - all.sessionChoice(all.sessionType == 2, :);

%%
rats = {'all','R042','R044','R050','R064'};
rat_idx = {1:24, 1:6, 7:12. 13:18, 19:24}; % session idxs to use for each rat

food_val = 1;
water_val = 2;

nShuf = 1000;

%% Use raw output to generate shuffled data

for iRat = 1:length(rats)
    
    %
    this_data = all.sessionChoice(rat_idx{iRat}, :);
    this_lat = all.sessionLatency(rat_idx{iRat}, :);
    this_sess = all.sessionType(rat_idx{iRat}, :);
    
    food_sess = this_sess == food_val; water_sess = this_sess == water_val;

    data.(rats{iRat}).food_fracL_evt = nanmean(this_data(food_sess, :), 2);
    data.(rats{iRat}).water_fracL_evt = nanmean(this_data(water_sess, :), 2);
    data.(rats{iRat}).food_lat = nanmean(this_lat(food_sess, :), 2);
    data.(rats{iRat}).water_lat = nanmean(this_lat(water_sess, :), 2);
        
    %
    for iShuf = nShuf:-1:1
        
        this_sess = all.sessionType(randperm(length(all.sessionType)));
        this_sess = this_sess(rat_idx{iRat}, :);
        
        food_sess = this_sess == food_val; water_sess = this_sess == water_val;
        
        this_food = nanmean(this_data(food_sess, :), 2); % shuffled food fracL for each food session
        this_water = nanmean(this_data(water_sess, :), 2); % shuffled water fracL for each water sessions
        
        foodShuf_fracL_evt(iShuf, :) = nanmean(this_food);
        waterShuf_fracL_evt(iShuf, :) = nanmean(this_water);
        diffShuf_fracL_evt(iShuf, :) = nanmean(this_food) - nanmean(this_water);
        
        foodShuf_lat(iShuf, :) = nanmean(nanmean(this_lat(food_sess, :), 2));
        waterShuf_lat(iShuf, :) = nanmean(nanmean(this_lat(water_sess, :), 2));
        diffShuf_lat(iShuf, :) = nanmean(nanmean(this_lat(food_sess, :), 2)) - nanmean(nanmean(this_lat(water_sess, :), 2));
        
    end
    
    data.(rats{iRat}).foodShuf_fracL_evt = nanmean(foodShuf_fracL_evt);
    data.(rats{iRat}).waterShuf_fracL_evt = nanmean(waterShuf_fracL_evt);
    data.(rats{iRat}).diffShuf_fracL_evt = nanmean(diffShuf_fracL_evt);
    
    data.(rats{iRat}).foodShufsd_fracL_evt = nanstd(foodShuf_fracL_evt);
    data.(rats{iRat}).waterShufsd_fracL_evt = nanstd(waterShuf_fracL_evt);
    data.(rats{iRat}).diffShufsd_fracL_evt = nanstd(diffShuf_fracL_evt);
    
    data.(rats{iRat}).foodShuf_lat = nanmean(foodShuf_lat);
    data.(rats{iRat}).waterShuf_lat = nanmean(waterShuf_lat);
    data.(rats{iRat}).diffShuf_lat = nanmean(diffShuf_lat);
    
    data.(rats{iRat}).foodShufsd_lat = nanstd(foodShuf_lat);
    data.(rats{iRat}).waterShufsd_lat = nanstd(waterShuf_lat);
    data.(rats{iRat}).diffShufsd_lat = nanstd(diffShuf_lat);

end

%% report some stats
this_z = (nanmean(data.all.food_lat) - data.all.foodShuf_lat) ./ data.all.foodShufsd_lat;
this_p = 2*normcdf(abs(this_z), 0, 1, 'upper');
fprintf('%s: %.2f +/- %.2f, z %.2f, p %.2e\n', 'food lat', nanmean(data.all.food_lat), nanstd(data.all.food_lat), this_z, this_p);

this_z = (nanmean(data.all.water_lat) - data.all.waterShuf_lat) ./ data.all.waterShufsd_lat;
this_p = 2*normcdf(abs(this_z), 0, 1, 'upper');
fprintf('%s: %.2f +/- %.2f, z %.2f, p %.2e\n', 'water lat', nanmean(data.all.water_lat), nanstd(data.all.water_lat), this_z, this_p);


%% plot

cfg_def = [];
cfg_def.showAllRatsText = 0;
cfg_def.colormode = 'inventory3';
cfg_def.colors = TmazeColors(cfg_def.colormode);
cfg_def.fs = 18; % font size
cfg_def.ylim = [1 2];
cfg_def.ylimtick = [0.5 0.5];
cfg_def.writeOutput = 0;
cfg_def.output_fn = 'temp';

cfg = ProcessConfig(cfg_def, []);

%%

what_idx = {[1 2 7 8]};

% first plot data for all rats
rats = {'all'}; iRat = 1;
col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w [0 0 0]}; % {foodColor waterColor diffColor}

location = [3 2 1]; % where on y-axis to place the data (food restriction, water restriction, difference respectively)
recColor = [0.8 0.8 0.8]; %[0.5 0.5 0.5];
ylims = [0.5 max(location)+0.5];
ylab = {'FR', 'WR', 'diff'};

% main plots

subplot(4, 6, what_idx{1});

% arrange the data
this_data_all{1} = data.(rats{iRat}).food_fracL_evt;
this_data_all{2} = data.(rats{iRat}).water_fracL_evt;
this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];

% arrange the shuffled data
this_shuf_mean = [data.(rats{iRat}).foodShuf_fracL_evt data.(rats{iRat}).waterShuf_fracL_evt data.(rats{iRat}).diffShuf_fracL_evt + 0.5];
this_shuf_sd = [data.(rats{iRat}).foodShufsd_fracL_evt data.(rats{iRat}).waterShufsd_fracL_evt data.(rats{iRat}).diffShufsd_fracL_evt];

% some stats
m1 = nanmean(this_data_all{1}); s1 = nanstd(this_data_all{1}) ./ sqrt(4); % SEM
m2 = nanmean(this_data_all{2}); s2 = nanstd(this_data_all{2}) ./ sqrt(4); % SEM
p = ranksum(this_data_all{1}, this_data_all{2});
fprintf('food %.2f +/ %.2f, water %.2f +/- %.2f, p = %.2e (ranksum)\n', m1, s1, m2, s2, p);

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
set(gca, 'XTick', 0:0.25:1, 'FontSize', cfg.fs, 'LineWidth', 1, 'XLim', [-0.5 1.5], 'XTickLabel', {'0 (W)', '', '0.5', '', '1 (F)'}, 'XDir', 'reverse');
set(gca, 'YLim' , ylims, 'YTick', 1:3, 'TickDir', 'out', 'YTickLabel', ylab(location));

set(gca,'Layer','top')

plot([0.5 0.5], ylim, '--', 'Color', [0.5 0.5 0.5]); % vertical line

% insets for individual rats
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 subplot start positions for PRE

for iRat = 1:length(rats)
    
    col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w [0 0 0]}; % {foodColor waterColor diffColor}

    subplot(4,6,start(iRat))
    
    % arrange the data
    this_data_all{1} = data.(rats{iRat}).food_fracL_evt;
    this_data_all{2} = data.(rats{iRat}).water_fracL_evt;
    this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
    this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];
    
    % arrange the shuffled data
    this_shuf_mean = [data.(rats{iRat}).foodShuf_fracL_evt data.(rats{iRat}).waterShuf_fracL_evt data.(rats{iRat}).diffShuf_fracL_evt + 0.5];
    this_shuf_sd = [data.(rats{iRat}).foodShufsd_fracL_evt data.(rats{iRat}).waterShufsd_fracL_evt data.(rats{iRat}).diffShufsd_fracL_evt];
    
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
    set(gca, 'XTick', 0:0.25:1, 'FontSize', 8, 'LineWidth', 1, 'XLim', [-0.5 1.5], 'XTickLabel', {'0 (W)', '', '0.5', '', '1 (F)'}, 'XDir', 'reverse');
    if iRat == 2
        %set(gca, 'XLim', [-0.25 1.25]);
    end
    
    set(gca, 'YLim' , ylims, 'YTick', 1:3, 'TickDir', 'out', 'YTickLabel', ylab(location));
    %title(sprintf('%s', what{iW})); box off
    set(gca,'Layer','top')
    
    plot([0.5 0.5], ylim, '--', 'Color', [0.5 0.5 0.5]);
    
end % of rats


%% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end