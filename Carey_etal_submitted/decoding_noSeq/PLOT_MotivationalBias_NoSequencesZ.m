function fh = PLOT_MotivationalBias_NoSequences(cfg_in,data)
% function fh = PLOT_Motivational_NoSequences(cfg_in,data)
%
% plots SWR content averaged over sessions as raw Z-scores
% (MotivationalT data set)
%
% INPUTS:
%
% cfg_in: options that control display properties
%
% data: struct with .all, .pre, .task, .post fields containing data output
% from PLOT_DecSeqCombined
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

cfg = ProcessConfig(cfg_def,cfg_in);

biasfun = @(d) (d(1)-d(2)); % computes bias measure as (food-water)

what = {'pre','task','post'}; % plot & analyze pre-task, task and post-task epochs separately
%what = {'all'}; % use this instead if you want to plot & analyze everything combined

what_idx = {[1 2 7 8],[3 4 9 10],[5 6 11 12]};

% first plot data for all rats
rats = {'all'}; iRat = 1;
col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w}; % {foodColor waterColor}

ylab = {'log odds p(left)/p(right) z-score'};
ylimsall = [-1 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [-2 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2]; % where to place the bar
xlims = [0 location(end)+1];

% main plots
for iW = 1:length(what)
    
    subplot(4,6,what_idx{iW});
    
    this_data_all{1} = data.(what{iW}).(rats{iRat}).food_median_z;
    this_data_all{2} = data.(what{iW}).(rats{iRat}).water_median_z;
    this_data = [nanmean(this_data_all{1}) nanmean(this_data_all{2})];
    
    % plot the data
    for iBar = 1:length(this_data)
        h(iBar) = bar(location(iBar),this_data(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
        plot(location(iBar),this_data_all{iBar},'.k');
    end
    set(gca,'XTick',location,'XTickLabel',{'F','W'},'FontSize',cfg.fs,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall,'TickDir','out')
    %xlabel('  food                  water','FontSize',cfg.fs)
    ylabel(ylab,'FontSize',cfg.fs)
    title(sprintf('%s %.2f',what{iW},biasfun(this_data)));
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',cfg.fs)
    end
    
    % some stats
    m1 = nanmean(this_data_all{1}); s1 = nanstd(this_data_all{1})./sqrt(4); % SEM
    m2 = nanmean(this_data_all{2}); s2 = nanstd(this_data_all{2})./sqrt(4); % SEM
    p = ranksum(this_data_all{1},this_data_all{2});
    fprintf('%s: food %.2f +/ %.2f, water %.2f +/- %.2f, p = %.2e\n',what{iW},m1,s1,m2,s2,p);
    
end % of whats for main plots

% insets for individual rats
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 subplot start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w cfg.colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    for iW = 1:length(what)
    
    subplot(4,6,position(iW))
    
    this_data = [nanmean(data.(what{iW}).(rats{iRat}).food_median_z) nanmean(data.(what{iW}).(rats{iRat}).water_median_z)];
    
    for iBar = 1:length(this_data)
        h(iBar) = bar(location(iBar),this_data(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        %ylabel(ylab,'FontSize',cfg.fs-2)
    else
        yticks = [];
    end
    
    %set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.fs-2,'LineWidth',1,'XLim',xlims);
    %set(gca,'YLim',ylimssing,'YTick',yticks)
    set(gca,'XTick',[],'XTickLabel',[],'FontSize',cfg.fs-2,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    %xlabel('  food   water','FontSize',cfg.fs-2)
    box off
    set(gca,'Layer','top')
    %txt = sprintf('%s %.2f',rats{iRat},biasfun(this_data));
    txt = sprintf('%.2f',biasfun(this_data));
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.fs-6)
    
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




