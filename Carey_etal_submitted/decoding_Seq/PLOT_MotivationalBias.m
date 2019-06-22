function fh = PLOT_MotivationalBias(cfg_in,data)
% function fh = PLOT_MotivationalBias(cfg_in,data)
%
% plots SWR content averaged over sessions
% (MotivationalT data set)
%
% INPUTS:
%
% cfg_in: options that control display properties
%
% data: struct with .all, .pre, .task, .post fields containing data output
% by ALL_CollectDecSeq.m, see ALL_Plot_DecSeq for an example of how to load
% these
%
% OUTPUTS:
%

cfg_def = [];
cfg_def.showAllRatsText = 0;
cfg_def.epochs = {'all'}; % {'all'} or {'pre','task','post'}

cfg = ProcessConfig(cfg_def,cfg_in);

biasfun = @(d) (d(1)-d(2))-(d(3)-d(4)); % computes bias measure as (food_left-food_right)-(water_left-water_right) sequence content proportions

what = cfg.epochs;
what_idx = {[1 2 7 8],[3 4 9 10],[5 6 11 12]};

% first plot data for all rats
rats = {'all'}; iRat = 1;
col = {cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).f cfg.colors.(rats{iRat}).w cfg.colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}

ylab = {'Number of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% main plots
for iW = 1:length(what)
    
    subplot(4,6,what_idx{iW});
    
    switch cfg.type
        case 'counts'
            this_data = [data.(what{iW}).data.(rats{iRat}).food_left data.(what{iW}).data.(rats{iRat}).food_right ...
                data.(what{iW}).data.(rats{iRat}).water_left data.(what{iW}).data.(rats{iRat}).water_right];
            session_data = IndivSessionData(data.(what{iW}).data.(rats{iRat}).ALL_sig_seq);
            
            fprintf('%s (all rats): %d sequences total\n', what{iW}, sum(this_data));
    
        case 'props'
            this_data = [data.(what{iW}).data.(rats{iRat}).food_leftN data.(what{iW}).data.(rats{iRat}).food_rightN ...
                data.(what{iW}).data.(rats{iRat}).water_leftN data.(what{iW}).data.(rats{iRat}).water_rightN];
            [~,session_data] = IndivSessionData(data.(what{iW}).data.(rats{iRat}).ALL_sig_seq);
        otherwise
            error('Unknown type %s',cfg.type);
    end
          

    % plot the data
    for iBar = 1:length(this_data)
        h(iBar) = bar(location(iBar),this_data(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
        plot(location(iBar),session_data(iBar,:),'.k');
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.fs,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall)
    xlabel('  food                  water','FontSize',cfg.fs)
    ylabel([ylab{1},' ',ylab{2}],'FontSize',cfg.fs)
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
    
    switch cfg.type
        case 'counts'
            this_data = [data.(what{iW}).data.(rats{iRat}).food_left data.(what{iW}).data.(rats{iRat}).food_right ...
                data.(what{iW}).data.(rats{iRat}).water_left data.(what{iW}).data.(rats{iRat}).water_right];
            
            fprintf('%s (%s): %d sequences\n', what{iW}, rats{iRat}, sum(this_data));
        case 'props'
            this_data = [data.(what{iW}).data.(rats{iRat}).food_leftN data.(what{iW}).data.(rats{iRat}).food_rightN ...
                data.(what{iW}).data.(rats{iRat}).water_leftN data.(what{iW}).data.(rats{iRat}).water_rightN];
        otherwise
            error('Unknown type %s',cfg.type);
    end
    
    for iBar = 1:length(this_data)
        h(iBar) = bar(location(iBar),this_data(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        ylabel(ylab,'FontSize',cfg.fs-2)
    else
        yticks = [];
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.fs-2,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',yticks)
    xlabel('  food   water','FontSize',cfg.fs-2)
    box off
    set(gca,'Layer','top')
    txt = sprintf('%s %.2f',rats{iRat},biasfun(this_data));
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.fs-2)
    
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


function [out,outN] = IndivSessionData(data)
% outputs 4 x nSessions individual session data for food-left, food-right, water-left, water-right
%
% input: ALL_sig_seq for data of interest
out = []; outN = [];
all_sessions = unique(data.sess);
for iS = 1:length(all_sessions)
    
    this_S = all_sessions(iS);
    
    this_out = nan(4,1); % food-left, food-right, water-left, water-right for this session
    this_outN = nan(4,1); % food-left, food-right, water-left, water-right for this session
    
    %for iArm = 1:2
        
        data_idx = find(data.sess == this_S);% & data.arm == iArm);
        this_data = data.count(data_idx); this_dataN = data.countN(data_idx);
        this_restrict = data.type(data_idx);
        
        if length(unique(this_restrict)) > 1, error('Session should not have more than one restriction type!'); end
        this_restrict = unique(this_restrict);
        
        switch this_restrict % check if food or water restriction
            case 1 % food
                this_out(1:2) = this_data; this_outN(1:2) = this_dataN;
            case 2 % water
                this_out(3:4) = this_data; this_outN(3:4) = this_dataN;
        end
        
    %end
    
    out = cat(2,out,this_out); outN = cat(2,outN,this_outN);
    
end

