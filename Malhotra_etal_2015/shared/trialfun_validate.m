function trialfun_validate(varargin)
% function trialfun_validate(varargin)
%
% write a figure with locations of cues, rewards etc. to test trialfun
%
% varargins:
%
% 'writeOutput': 0 or 1 to write .png image file in current folder
%
% MvdM 2013-07-03

writeOutput = 0;
extract_varargin;

evts = {'nosepoke','reward','cue'};
block = 'value';
locs = {'left','right','both'};
cues = {'c1','c3','c5'}; % {'c1','c3','c5','lo','hi'};

cl = 'rgbcm';
load(FindFile('*vt.mat'))
run(FindFile('*keys.m'));

x = Restrict(x,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2));
y = Restrict(y,ExpKeys.TimeOnTrack(1),ExpKeys.TimeOffTrack(2));

figure; iPlot = 1;

for iE = 1:numel(evts)
    for iL = 1:numel(locs)
        
        subaxis(numel(evts),numel(locs),iPlot, 'Spacing', 0.03, 'Padding', 0.03, 'Margin', 0.03);
        plot(Data(x),Data(y),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]); hold on; axis tight;
        axis off
        
        nTrials = 0;
        
        clear h;
        for iC = 1:numel(cues)
            
            cfg.trialdef.eventtype = evts{iE}; % could be 'nosepoke', 'reward', 'cue'
            cfg.trialdef.location = locs{iL}; % could be 'left', 'right', 'both'
            cfg.trialdef.block = block; % could be 'value', 'risk'
            cfg.trialdef.cue = cues(iC); % cell array with choice of elements {'c1','c3','c5','lo','hi'}
            cfg.ExpKeys = ExpKeys;
            
            [trl,evt] = ft_trialfun_lineartracktone2(cfg,'OutputFormat','Time','EvtOnly',1);
            
            x_interp = interp1(Range(x),Data(x),evt.ts);
            y_interp = interp1(Range(y),Data(y),evt.ts);
            
            h(iC) = plot(x_interp,y_interp,'.','Color',cl(iC),'MarkerSize',10);
            nTrials = nTrials + length(evt.ts);
        end
        
        if iPlot == 8
            l = legend(h,cues,'Orientation','horizontal','Location','Southwest','FontSize',5); legend boxoff;
        end
        t = title(sprintf('%s - %s - n = %d',evts{iE},locs{iL},nTrials)); set(t,'FontSize',8);
        
        set(gca,'YLim',[50 250]);
        
        iPlot = iPlot + 1;
        
    end
end

%% write output
if writeOutput
    print(gcf,'-dpng','-r300','-opengl','trial_validate.png');
end