%%
whatP = {'paramRec_gain1.mat','paramRec_gain2.mat','paramRec_drift.mat'}; % ALL_paramRec variables saved by MASTER_paramRecovery.m
true_param{1} = [0.95 0.975 1 1.025 1.05];
true_param{2} = [0.95 0.975 1 1.025 1.05];
true_param{3} = [-1 -0.5 0 0.5 1];

whichS = 'kalmanwrapped';

%% collect data
ALL_param = nan(length(whatP),length(gain_vals),15);
for iP = 1:length(whatP)
    
    load(whatP{iP}); % ALL_paramRec variable saved by MASTER_paramRecovery.m
    
    for iVal = 1:length(true_param{iP})
    
        for iS = 1:length(ALL_paramRec{iVal})
            
            this_param = ALL_paramRec{iVal}(iS).laser.(whichS).param{4};
            
            ALL_param(iP,iVal,iS) = nanmean(this_param(:,iP));
            
        end % of sessions
        
    end % of param values
    
end

%% plot
keep_idx = [5 7 8 9 10 12 14];
mean_width = 0.25;
for iP = 1:3
    
    subplot(2,3,iP);
    
    plot(true_param{iP},sq(ALL_param(iP,:,:)),'.r');
    hold on;
    plot(true_param{iP},sq(ALL_param(iP,:,keep_idx)),'.k');
    
    m = sq(nanmean(ALL_param(iP,:,:),3));
    
    d = median(diff(true_param{iP}));
    
    for iV = 1:length(true_param{iP})
        
        plot([true_param{iP}(iV)-mean_width*d true_param{iP}(iV)+mean_width*d],[m(iV) m(iV)],'LineWidth',2,'Color','k');
        
    end
    
    if iP == 3
        plot([-1.5 1.5],[-1.5 1.5],'k--');
        set(gca,'FontSize',18,'LineWidth',1,'TickDir','out','XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'XTick',-1.5:0.5:1.5,'YTick',-1.5:0.5:1.5); box off;
        xlabel('true \delta'); ylabel('recovered \delta');
    else
        plot([0.9 1.1],[0.9 1.1],'k--');
        set(gca,'FontSize',18,'LineWidth',1,'TickDir','out','XLim',[0.9 1.1],'YLim',[0.9 1.1],'XTick',0.9:0.05:1.1,'YTick',0.9:0.05:1.1); box off;
        xlabel('true \gamma'); ylabel('recovered \gamma');
    end
    
end