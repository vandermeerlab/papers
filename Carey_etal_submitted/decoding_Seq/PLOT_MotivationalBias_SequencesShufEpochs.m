function data = PLOT_MotivationalBias_SequencesShufEpochs(cfg_in)
% function data = PLOT_Motivational_SequencesShufEpochs(cfg_in)
%
% plots SWR sequence content averaged over sessions
% (MotivationalT data set)
%
% INPUTS:
%
% cfg_in: options that control display properties
%
% assumes that outputs of ALL_Collect_DecSeq_eligibleShuf.m have been saved in target folder
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
cfg_def.epochs = {'all'}; % {'all'} or {'pre', 'task', 'post'};

cfg = ProcessConfig(cfg_def, cfg_in);
rng('default');

%% load data somehow
temp = load('C:\temp\S1_DecSeq_prerecord_all_eligibleCP_out'); data.pre = temp.data;
temp = load('C:\temp\S1_DecSeq_taskrest_all_eligibleCP_out'); data.task = temp.data;
temp = load('C:\temp\S1_DecSeq_postrecord_all_eligibleCP_out'); data.post = temp.data;
temp = load('C:\temp\S1_DecSeq_all_all_eligibleCP_out'); data.all = temp.data;

%% first arrange data and do shuffles
cfg.rats = {'R042','R044','R050','R064'};
toDo = cat(2,{'all'},cfg.rats);

for iEpoch = 1:length(cfg.epochs)
    for iDo = 1:length(toDo)
        
        this_e = cfg.epochs{iEpoch};
        this_do = toDo{iDo};
        
        data.(this_e).(this_do).ALL_sig_seq = removeNaNs(data.(this_e).(this_do).ALL_sig_seq); % remove sessions with zero counts
        
        % actual data
        data.(this_e).(this_do).food_left = data.(this_e).(this_do).ALL_sig_seq.count(data.(this_e).(this_do).ALL_sig_seq.arm == 1 & data.(this_e).(this_do).ALL_sig_seq.type == 1);
        data.(this_e).(this_do).food_right = data.(this_e).(this_do).ALL_sig_seq.count(data.(this_e).(this_do).ALL_sig_seq.arm == 2 & data.(this_e).(this_do).ALL_sig_seq.type == 1);
        data.(this_e).(this_do).food_leftN = data.(this_e).(this_do).ALL_sig_seq.countN(data.(this_e).(this_do).ALL_sig_seq.arm == 1 & data.(this_e).(this_do).ALL_sig_seq.type == 1);
        data.(this_e).(this_do).food_rightN = data.(this_e).(this_do).ALL_sig_seq.countN(data.(this_e).(this_do).ALL_sig_seq.arm == 2 & data.(this_e).(this_do).ALL_sig_seq.type == 1);
        
        data.(this_e).(this_do).water_left = data.(this_e).(this_do).ALL_sig_seq.count(data.(this_e).(this_do).ALL_sig_seq.arm == 1 & data.(this_e).(this_do).ALL_sig_seq.type == 2);
        data.(this_e).(this_do).water_right = data.(this_e).(this_do).ALL_sig_seq.count(data.(this_e).(this_do).ALL_sig_seq.arm == 2 & data.(this_e).(this_do).ALL_sig_seq.type == 2);
        data.(this_e).(this_do).water_leftN = data.(this_e).(this_do).ALL_sig_seq.countN(data.(this_e).(this_do).ALL_sig_seq.arm == 1 & data.(this_e).(this_do).ALL_sig_seq.type == 2);
        data.(this_e).(this_do).water_rightN = data.(this_e).(this_do).ALL_sig_seq.countN(data.(this_e).(this_do).ALL_sig_seq.arm == 2 & data.(this_e).(this_do).ALL_sig_seq.type == 2);
        
        % shuffles
        nShuf = 1000;
        shuf_fl = []; shuf_fln = []; shuf_fr = []; shuf_frn = []; shuf_diffn = [];
        shuf_wl = []; shuf_wln = []; shuf_wr = []; shuf_wrn = [];
        for iShuf = nShuf:-1:1
            
            this_type = data.(this_e).(this_do).ALL_sig_seq.type; temp_idx = nan(size(this_type));
            
            nPairs = length(this_type) ./ 2; % need to shuffle sessions (2 arms per session)
            pair_idx = randperm(nPairs);
            
            temp_idx1 = 2 * pair_idx - 1; temp_idx2 = 2 * pair_idx;
            temp_idx(1:2:end) = temp_idx1; temp_idx(2:2:end) = temp_idx2;
            this_type = this_type(temp_idx);
            
            this_arm = data.(this_e).(this_do).ALL_sig_seq.arm;
            
            shuf_fl(iShuf) = nansum(data.(this_e).(this_do).ALL_sig_seq.count(this_arm == 1 & this_type == 1));
            shuf_fr(iShuf) = nansum(data.(this_e).(this_do).ALL_sig_seq.count(this_arm == 2 & this_type == 1));
            shuf_fln(iShuf) = nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 1 & this_type == 1));
            shuf_frn(iShuf) = nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 2 & this_type == 1));
            
            shuf_wl(iShuf) = nansum(data.(this_e).(this_do).ALL_sig_seq.count(this_arm == 1 & this_type == 2));
            shuf_wr(iShuf) = nansum(data.(this_e).(this_do).ALL_sig_seq.count(this_arm == 2 & this_type == 2));
            shuf_wln(iShuf) = nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 1 & this_type == 2));
            shuf_wrn(iShuf) = nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 2 & this_type == 2));
            
            shuf_diffn(iShuf) = nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 1 & this_type == 1)) - ...
                nanmean(data.(this_e).(this_do).ALL_sig_seq.countN(this_arm == 1 & this_type == 2));
            
        end
        
        data.(this_e).(this_do).food_leftS = nanmean(shuf_fl); data.(this_e).(this_do).food_leftSd = nanstd(shuf_fl);
        data.(this_e).(this_do).food_rightS = nanmean(shuf_fr); data.(this_e).(this_do).food_rightSd = nanstd(shuf_fr);
        data.(this_e).(this_do).water_leftS = nanmean(shuf_wl); data.(this_e).(this_do).water_leftSd = nanstd(shuf_wl);
        data.(this_e).(this_do).water_rightS = nanmean(shuf_wr); data.(this_e).(this_do).water_rightSd = nanstd(shuf_wr);
        
        data.(this_e).(this_do).food_leftNS = nanmean(shuf_fln); data.(this_e).(this_do).food_leftNSd = nanstd(shuf_fln);
        data.(this_e).(this_do).food_rightNS = nanmean(shuf_frn); data.(this_e).(this_do).food_rightNSd = nanstd(shuf_frn);
        data.(this_e).(this_do).water_leftNS = nanmean(shuf_wln); data.(this_e).(this_do).water_leftNSd = nanstd(shuf_wln);
        data.(this_e).(this_do).water_rightNS = nanmean(shuf_wrn); data.(this_e).(this_do).water_rightNSd = nanstd(shuf_wrn);
        
        data.(this_e).(this_do).diffNS = nanmean(shuf_diffn); data.(this_e).(this_do).diffNSd = nanstd(shuf_diffn);
        
    end
    
end
%%
what = cfg.epochs;
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
    this_e = what{iW};
    
    % arrange the data
    this_data_all{1} = data.(this_e).(rats{iRat}).food_leftN;
    this_data_all{2} = data.(this_e).(rats{iRat}).water_leftN;
    this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
    this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];
    
    % arrange the shuffled data
    this_shuf_mean = [data.(this_e).(rats{iRat}).food_leftNS data.(this_e).(rats{iRat}).water_leftNS data.(this_e).(rats{iRat}).diffNS + 0.5];
    this_shuf_sd = [data.(this_e).(rats{iRat}).food_leftNSd data.(this_e).(rats{iRat}).water_leftNSd data.(this_e).(rats{iRat}).diffNSd];
    
    % some stats
    m1 = nanmean(this_data_all{1}); s1 = nanstd(this_data_all{1}) ./ sqrt(4); % SEM
    m2 = nanmean(this_data_all{2}); s2 = nanstd(this_data_all{2}) ./ sqrt(4); % SEM
    p = ranksum(this_data_all{1}, this_data_all{2});
    fprintf('%s: food %.2f +/ %.2f, water %.2f +/- %.2f, p = %.2e (ranksum)\n', what{iW}, m1, s1, m2, s2, p);
    
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

%%

% insets for individual rats
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 subplot start positions for PRE

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];

    col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w [0 0 0]}; % {foodColor waterColor diffColor}

    for iW = 1:length(what)
    
        subplot(4,6,position(iW))
        this_e = what{iW};
        
        % arrange the data
        this_data_all{1} = data.(this_e).(rats{iRat}).food_leftN;
        this_data_all{2} = data.(this_e).(rats{iRat}).water_leftN;
        this_data_all{3} = []; % individual session "food minus water" differences are not (yet) defined
        this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2}) nanmean(this_data_all{1}) - nanmean(this_data_all{2}) + 0.5];
        
        % arrange the shuffled data
        this_shuf_mean = [data.(this_e).(rats{iRat}).food_leftNS data.(this_e).(rats{iRat}).water_leftNS data.(this_e).(rats{iRat}).diffNS + 0.5];
        this_shuf_sd = [data.(this_e).(rats{iRat}).food_leftNSd data.(this_e).(rats{iRat}).water_leftNSd data.(this_e).(rats{iRat}).diffNSd];
    
        % plot the shuffled data
        for iRec = 1:length(this_shuf_mean) % [x y w h]
            
            rec_h(iRec) = rectangle('Position', [this_shuf_mean(iRec) - 0.5*this_shuf_sd(iRec) location(iRec) - 0.25 this_shuf_sd(iRec) 0.5]);
            set(rec_h(iRec), 'FaceColor', recColor, 'EdgeColor', 'none'); hold on;
            
            %rec_h(iRec) = errorbar(this_shuf_mean(iRec), location(iRec), 0.5*this_shuf_sd(iRec), 'horizontal'); hold on;
            %set(rec_h(iRec), 'LineWidth', 0.5, 'Color', recColor); box off;
            
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
            set(gca, 'XLim', [-0.25 1.25]);
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

end % of main function

%%
function data = removeNaNs(data)

keep_idx = ~isnan(data.countN);

fn = fieldnames(data);
for iF = 1:length(fn)
    
    data.(fn{iF}) = data.(fn{iF})(keep_idx);
end

end