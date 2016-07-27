%% PLOT_Classify.m

%% parameters
rats = {'R117','R119','R131','R132'};
classifier = 'linearDiscriminant'; % {'linearDiscriminant','naiveBayes','knn9'}
what = 'min'; % {'min','avg'}
nMaxSessions = 70; % only used for initializing arrays to be filled

%%
available_rats = fieldnames(ALL_evt);
nRats = 0; nSessions = 0;

cellCount = 1:5:21;
%cellCount = 1:5:31;
DATA_err = nan(length(cellCount),nMaxSessions);
DATA_errCtrl = nan(length(cellCount),nMaxSessions);
DATA_errBest = nan(nMaxSessions,1);
DATA_ratID = nan(length(cellCount),nMaxSessions); % useful for later when exporting to R, etc..
DATA_nCells = nan(length(cellCount),nMaxSessions);

% first, collect the data
for iRat = 1:length(rats)
    
    this_rat = rats{iRat};
    
    if ~strmatch(this_rat,available_rats)
       warning('Rat %s not available -- skipping...',rats{iRat});
       continue;
    end
    
    available_sessions = fieldnames(ALL_evt.(this_rat));
    
    for iSession = 1:length(available_sessions)
    
        this_session = available_sessions{iSession};
        this_data = ALL_evt.(this_rat).(this_session).classify;
        
        iClass = strmatch(classifier,this_data.classifier_list,'exact');
        nSessions = nSessions + 1;
        
        [~,keep_idx] = intersect(this_data.nCells,cellCount);
        %
        switch what
            case 'avg'
                DATA_err(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = this_data.cdata(iClass,keep_idx);
                DATA_errCtrl(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = this_data.cdataC(iClass,keep_idx);
                DATA_errBest(nSessions) = this_data.bestSingleNeuron(iClass);
            case 'min'
                DATA_err(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = this_data.cdata_min(iClass,keep_idx);
                DATA_errCtrl(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = this_data.cdataC_min(iClass,keep_idx);
        end
        
        % also store some useful things for doing stats later
        DATA_ratID(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = iRat;
        DATA_nCells(1:length(this_data.cdata(iClass,keep_idx)),nSessions) = cellCount(keep_idx);
    end % sessions
    
end % rats

%% plot (average)

% compute SEM
nSessions_cell = sum(~isnan(DATA_err),2);
SEM_session = nanstd(DATA_err,[],2)./sqrt(4);
SEM_sessionCtrl = nanstd(DATA_errCtrl,[],2)./sqrt(4);

% find out how many sessions we have for each cell count - need this to set
% axis limits etc..
%cc_noCell = cellCount(min(find(nSessions_cell == 0))); % first cell count with no cells

eb_jitter = 0.01; % so that errorbars don't overlap

clear h;
h(1) = plot(cellCount+eb_jitter,nanmean(DATA_err,2),'k','LineWidth',2);
hold on;
plot(cellCount+eb_jitter,nanmean(DATA_err,2),'k.','MarkerSize',20);
ebh = errorbar(cellCount+eb_jitter,nanmean(DATA_err,2),SEM_session); set(ebh,'Color','k');

h(2) = plot(cellCount-eb_jitter,nanmean(DATA_errCtrl,2),'Color',[0.5 0.5 0.5],'LineWidth',2);
hold on;
plot(cellCount-eb_jitter,nanmean(DATA_errCtrl,2),'.','Color',[0.5 0.5 0.5],'MarkerSize',20);
ebh = errorbar(cellCount-eb_jitter,nanmean(DATA_errCtrl,2),SEM_sessionCtrl); set(ebh,'Color',[0.5 0.5 0.5]);

%axis([0 cc_noCell-4 0 0.55]);
axis([0 cellCount(end)+1 0 0.55]);

if strcmp(what,'avg') % also plot best single neuron
   
    plot(cc_noCell,nanmean(DATA_errBest),'k.','MarkerSize',20);
    axis([0 cc_noCell+1 0 0.55]);
    
end

plot(cellCount,0.5*ones(size(cellCount)),'k:');

set(gca,'FontSize',24,'LineWidth',1); box off;
legend(h,{'actual events','control (+2s)'},'Location','Southeast'); legend boxoff;

xlabel('number of cells'); ylabel('cross-validation error');
set(gca,'YTick',[0:0.1:0.5],'YTickLabel',{'0 (100% correct)','0.1','0.2','0.3','0.4','0.5 (chance)'}, ...
    'XTick',cellCount,'XTickLabel',cellCount);

%% stats - quick anova in matlab - save stats_out variable for analysis in R
stats_out.classif = cat(1,DATA_err(:),DATA_errCtrl(:));
stats_out.nCells = cat(1,DATA_nCells(:),DATA_nCells(:));
stats_out.actualCtrl = cat(1,zeros(size(DATA_nCells(:))),ones(size(DATA_nCells(:))));
stats_out.ratID = cat(1,DATA_ratID(:),DATA_ratID(:));

[P,T,STATS,TERMS] = anovan(stats_out.classif,{stats_out.nCells stats_out.actualCtrl},'varnames',{'nCells','offset'},'model','full');

%% plot (rats separately)

for iRat = 1:4

    keep_idx = DATA_ratID(1,:) == iRat;
    subplot(2,2,iRat)
    
    % compute SEM
    nSessions_cell = sum(~isnan(DATA_err(:,keep_idx)),2);
    SEM_session = nanstd(DATA_err(:,keep_idx),[],2)./sqrt(sum(keep_idx));
    SEM_sessionCtrl = nanstd(DATA_errCtrl(:,keep_idx),[],2)./sqrt(sum(keep_idx));
    
    % find out how many sessions we have for each cell count - need this to set
    % axis limits etc..
    %cc_noCell = cellCount(min(find(nSessions_cell == 0))); % first cell count with no cells
    
    eb_jitter = 0.01; % so that errorbars don't overlap
    
    clear h;
    h(1) = plot(cellCount+eb_jitter,nanmean(DATA_err(:,keep_idx),2),'k','LineWidth',2);
    hold on;
    plot(cellCount+eb_jitter,nanmean(DATA_err(:,keep_idx),2),'k.','MarkerSize',20);
    ebh = errorbar(cellCount+eb_jitter,nanmean(DATA_err(:,keep_idx),2),SEM_session); set(ebh,'Color','k');
    
    h(2) = plot(cellCount-eb_jitter,nanmean(DATA_errCtrl(:,keep_idx),2),'Color',[0.5 0.5 0.5],'LineWidth',2);
    hold on;
    plot(cellCount-eb_jitter,nanmean(DATA_errCtrl(:,keep_idx),2),'.','Color',[0.5 0.5 0.5],'MarkerSize',20);
    ebh = errorbar(cellCount-eb_jitter,nanmean(DATA_errCtrl(:,keep_idx),2),SEM_sessionCtrl); set(ebh,'Color',[0.5 0.5 0.5]);
    
    %axis([0 cc_noCell-4 0 0.55]);
    axis([0 cellCount(end)+1 0 0.55]);
    
    if strcmp(what,'avg') % also plot best single neuron
        
        plot(cc_noCell,nanmean(DATA_errBest),'k.','MarkerSize',20);
        axis([0 cc_noCell+1 0 0.55]);
        
    end
    
    plot(cellCount,0.5*ones(size(cellCount)),'k:');
    
    set(gca,'FontSize',18,'LineWidth',1); box off;
    legend(h,{'actual','control'},'Location','Southeast'); legend boxoff;
    
    xlabel('number of cells'); ylabel('cross-validation error');
    set(gca,'YTick',[0:0.1:0.5],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'}, ...
        'XTick',cellCount,'XTickLabel',cellCount);
    
    title(rats{iRat});

end
%% single sessions...
this_data = ALL_evt.R117.R117_2007_06_21.classify;
iClass = 1;

for iC = 1%:3

    figure(iC);
    
    subplot(121);
    h(1) = plot(this_data.nCells,this_data.cdata(iC,:),'b','LineWidth',2);
    hold on;
    plot(this_data.nCells,this_data.cdata(iC,:),'.b','MarkerSize',20);
    
    h(2) = plot(this_data.nCells,this_data.cdataC(iC,:),'r','LineWidth',2);
    plot(this_data.nCells,this_data.cdataC(iC,:),'.r','MarkerSize',20);
    
    h(3) = plot(this_data.nCells,this_data.cdataS(iC,:),'g:','LineWidth',2);
    
    plot(this_data.nCells,0.5*ones(size(this_data.nCells)),'k:');
    axis([0 max(this_data.nCells)+3 0 0.55]);
    
    plot(max(this_data.nCells)+2,this_data.bestSingleNeuron(iC),'.k','MarkerSize',20);
    
    title([classifier_list{iC} ' - avg']);

    % min
    subplot(122);
    h(1) = plot(this_data.nCells,this_data.cdata_min(iC,:),'b','LineWidth',2);
    hold on;
    plot(this_data.nCells,this_data.cdata_min(iC,:),'.b','MarkerSize',20);
    
    h(2) = plot(this_data.nCells,this_data.cdataC_min(iC,:),'r','LineWidth',2);
    plot(this_data.nCells,this_data.cdataC_min(iC,:),'.r','MarkerSize',20);
    
    h(3) = plot(this_data.nCells,this_data.cdataS_min(iC,:),'g:','LineWidth',2);
    
    plot(this_data.nCells,0.5*ones(size(this_data.nCells)),'k:');
    axis([0 max(this_data.nCells)+3 0 0.55]);
    
    plot(max(this_data.nCells)+2,this_data.bestSingleNeuron(iC),'.k','MarkerSize',20);

    title([classifier_list{iC} ' - min']);
    
end