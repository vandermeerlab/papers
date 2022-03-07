%% set up path
clear all;

MASTER_root = 'D:\My_Documents\Dropbox\projects\Alyssa';
cd(MASTER_root);
MASTER_path; % reset and then set up path for this project

%% which data files to use?

what = {'isidro_samelap_nonmatched_20160705.mat','isidro_nextlap_nonmatched_20160705.mat','isidro_loocv_nonmatched_20160705.mat'};
whatL = {'same-lap','next-lap','leave-one-out'}; % labels for plotting
    
%% config
cfg_plot = [];
cfg_plot.lr = {'left','right'};
cfg_plot.fs = 16;
cfg_plot.xltc = [0.5 111.5];
cfg_plot.yltc = [0 100];
cfg_plot.ylerr = [0 40];
cfg_plot.whichTC = 2;
cfg_plot.fnadd = '_all';
cfg_plot.nBins = 111;
cfg_plot.pPassed = 0.8; % parameter combo needs to exceed this threshold for included time bins (based on minimum activity)
cfg_plot.pnan = 0.2; % parameter combo needs to stay below this threshold for excluded time bins (based on decodability)
cfg_plot.minCells = 20; % both L and R must have this number of cells included
cfg_plot.writeOutput = 0;
cfg_plot.nMinSessions = 6;
cfg_plot.cax = [0 40];
cfg_plot.space2 = [4 4];% TC, Q idx into ALL_decErr for comparison to "best"
cfg_plot.nsem = 4;

%% collect data
for iW = 1:length(what)
    
    load(what{iW});
    
    % init vars for this file
    nSessions = length(ALL_decErr.fd);
    
    AVG.tc{1} = []; % tuning curves; 1 is left, 2 right as defined in cfg_plot.lr
    AVG.tc{2} = [];
    AVG.tc_sess{1} = []; % also track session each TC came from
    AVG.tc_sess{2} = []; % also track session each TC came from
    
    AVG.pcp{1} = nan(nSessions,cfg_plot.nBins); % histogram of place cell peaks
    AVG.pcp{2} = nan(nSessions,cfg_plot.nBins);
    
    AVG.des{1} = nan(nSessions,cfg_plot.nBins); AVG.des2{1} = nan(nSessions,cfg_plot.nBins); % dec err by space
    AVG.des{2} = nan(nSessions,cfg_plot.nBins); AVG.des2{2} = nan(nSessions,cfg_plot.nBins);
    
    AVG.occ{1} = nan(nSessions,cfg_plot.nBins); % occ
    AVG.occ{2} = nan(nSessions,cfg_plot.nBins);
    
    AVG.dpar{1} = nan(nSessions,length(ALL_decErr.cfg.TCsmooth),length(ALL_decErr.cfg.QsmoothSD));
    AVG.dpar{2} = nan(nSessions,length(ALL_decErr.cfg.TCsmooth),length(ALL_decErr.cfg.QsmoothSD));
    
    AVG.cmat{1} = nan(nSessions,cfg_plot.nBins,cfg_plot.nBins);
    AVG.cmat{2} = nan(nSessions,cfg_plot.nBins,cfg_plot.nBins);
    
    % loop over sessions, collect data
    for iS = 1:nSessions
        
        % check if should include
        if length(ALL_decErr.left.S{iS}.t) < cfg_plot.minCells | length(ALL_decErr.right.S{iS}.t) < cfg_plot.minCells
            fprintf('Session %d excluded (minCells)...\n',iS);
            continue;
        end
        
        %%
        for iLR = 1:length(cfg_plot.lr)
            
            % first, exclude anything that fails criteria
            this_err = sq(ALL_decErr.(cfg_plot.lr{iLR}).meanErr(iS,:,:));
            this_Ppassed = sq(ALL_decErr.(cfg_plot.lr{iLR}).Ppassed(iS,:,:));
            this_Pnan = sq(ALL_decErr.(cfg_plot.lr{iLR}).Pnan(iS,:,:));
            this_err(this_Ppassed < cfg_plot.pPassed) = NaN;
            this_err(this_Pnan > cfg_plot.pnan) = NaN;
            
            % decoding error by space - first find best overall
            [~,min_idx] = min(this_err(:));
            [min_x,min_y] = ind2sub(size(this_err),min_idx);
            
            this_spaceZ = sq(ALL_decErr.(cfg_plot.lr{iLR}).spaceErr(iS,min_x,min_y,:));
            AVG.des{iLR}(iS,:) = this_spaceZ;
            
            this_spaceZ = sq(ALL_decErr.(cfg_plot.lr{iLR}).spaceErr(iS,cfg_plot.space2(1),cfg_plot.space2(2),:));
            AVG.des2{iLR}(iS,:) = this_spaceZ;
            
            AVG.dpar{iLR}(iS,:,:) = this_err;
            
            this_confmat = ALL_decErr.(cfg_plot.lr{iLR}).confMat{iS,min_x,min_y}.thr';
            AVG.cmat{iLR}(iS,:,:) = this_confmat;
            
        end % of left/right
        
    end % of sessions
    
    % obtain and keep grand error matrix
    clear temp;
    temp(1,:,:,:) = AVG.dpar{1};
    temp(2,:,:,:) = AVG.dpar{2};
    grand_av{iW} = sq(nanmean(temp,1)); % within-session L/R average
    %grand_av = sq(nanmean(cat(1,AVG.dpar{1},AVG.dpar{2})));
    
    % only keep parameters with some minimum number of sessions
    clear nSess_temp;
    for iI = 1:size(grand_av{iW},2)
        for iJ = 1:size(grand_av{iW},3)
            nSess_temp(iI,iJ) = sum(~isnan(sq(grand_av{iW}(:,iI,iJ))));
            if nSess_temp(iI,iJ) < cfg_plot.nMinSessions
                grand_av{iW}(:,iI,iJ) = NaN;
            end
        end
    end
    
    this_grand_av = sq(nanmean(grand_av{iW}));
    this_grand_sem = sq(nanstd(grand_av{iW}))./sqrt(cfg_plot.nsem);
    
    %%% mean error %%%
    
    figure(1);
    
    subplot(2,3,iW)
    nan_imagesc_ec(this_grand_av');
    set(gca,'XTick',1:size(this_grand_av',2),'XTickLabel',ALL_decErr.cfg.TCsmooth);
    set(gca,'YTick',1:size(this_grand_av',1),'YTickLabel',ALL_decErr.cfg.QsmoothSD,'TickDir','out');
    xlabel('TC smoothing (bins)'); ylabel('Q smoothing SD (s)');
    
    hold on; box off;
    
    [~,min_idx] = min(this_grand_av(:));
    [min_x,min_y] = ind2sub(size(this_grand_av),min_idx);
    
    hold on;
    plot(min_x,min_y,'*','Color',[1 1 1],'MarkerSize',10);
    
    title(sprintf('%s err %.2f +/- %.2f',whatL{iW},this_grand_av(min_x,min_y),this_grand_sem(min_x,min_y)));
    caxis(cfg_plot.cax);
    colorbar;
    
    % confusion
    clear temp;
    temp(1,:,:,:) = AVG.cmat{1};
    temp(2,:,:,:) = AVG.cmat{2};
    this_cmat = sq(nanmean(temp,1)); % av L/R
    
    % only keep parameters with some minimum number of sessions
    clear nSess_temp;
    for iI = 1:size(grand_av{iW},2)
        for iJ = 1:size(grand_av{iW},3)
            nSess_temp(iI,iJ) = sum(~isnan(sq(this_cmat(:,iI,iJ))));
            if nSess_temp(iI,iJ) < cfg_plot.nMinSessions
                this_cmat(:,iI,iJ) = NaN;
            end
        end
    end
    
    this_cmat = sq(nanmean(this_cmat)); % av sessions
    
    subplot(2,3,iW+3);
    imagesc(1:111,1:111,this_cmat);
    colorbar; caxis([0 0.5]);
    set(gca,'XTick',[],'YTick',[]);
    xlabel('actual'); ylabel('decoded'); box off;
    
    %%% space error %%%
    figure(2);

    clear temp;
    temp(1,:,:) = AVG.des{1};
    temp(2,:,:) = AVG.des{2};
    this_spaceErr = sq(nanmean(temp,1)); % av L/R
    this_spaceErrZ = (this_spaceErr-repmat(nanmean(this_spaceErr,2),[1 size(this_spaceErr,2)]))./repmat(nanstd(this_spaceErr,[],2),[1 size(this_spaceErr,2)]);

    temp(1,:,:) = AVG.des2{1};
    temp(2,:,:) = AVG.des2{2};
    this_spaceErr2 = sq(nanmean(temp,1)); % av L/R
    this_spaceErrZ2 = (this_spaceErr2-repmat(nanmean(this_spaceErr2,2),[1 size(this_spaceErr2,2)]))./repmat(nanstd(this_spaceErr2,[],2),[1 size(this_spaceErr2,2)]);
    
    subplot(2,3,iW);
    plot(nanmean(this_spaceErr),'k'); hold on;
    %plot(nanmean(this_spaceErr2),'b');
    set(gca,'LineWidth',1,'FontSize',cfg_plot.fs,'XLim',[0 112],'YLim',[0 40],'XTick',[]);
    box off; xlabel('location'); ylabel('decoding error (cm)');
    title(whatL{iW});
    
    subplot(2,3,iW+3);
    plot(nanmean(this_spaceErrZ),'k'); hold on;
    %plot(nanmean(this_spaceErrZ2),'b');
    set(gca,'LineWidth',1,'FontSize',cfg_plot.fs,'XLim',[0 112],'YLim',[-2 2],'XTick',[]);
    box off; xlabel('location'); ylabel('decoding error (Z)');
    
end

%%
if cfg_plot.writeOutput
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300','decQuantify_multi1.png');
    print(gcf,'-depsc','decQuantify_multi1.eps');
end
%% compare for specific parameters
figure(5);

%what = {[0 0],[0.1 0.002],[0.25 0.005],[0.5 0.01],[1 0.025] [5 0.2]}; % works for all sessions
what2 = {[0 0],[0.5 0.005],[1 0.01],[1 0.05] [3 0.1]}; % works for all sessions
whatL = {'same-lap','next-lap','leave-one-out'};
 
cfg_plot2 = [];
cfg_plot2.col = 'rgbcmk';
cfg_plot2.x = 1:3;
cfg_plot2.xoff = -0.1:0.05:0.15;
cfg_plot2.ms = 20;
cfg_plot2.fs = 10;

clear whatP this_data this_dataSD h
for iW = 1:length(what2)
    
    this_data = []; this_dataSD = [];
    
    % get data for these params
    whatP{iW} = sprintf('sTC %.2f, sQ %.3f',what2{iW}(1),what2{iW}(2));
    for iD = 1:3
        
        this_x = find(ALL_decErr.cfg.TCsmooth == what2{iW}(1));
        this_y = find(ALL_decErr.cfg.QsmoothSD == what2{iW}(2));
        
        this_data(iD) = nanmean(grand_av{iD}(:,this_x,this_y));
        this_nSess = sum(~isnan((grand_av{iD}(:,this_x,this_y))));
        this_dataSD(iD) = nanstd(grand_av{iD}(:,this_x,this_y),[],1)/sqrt(this_nSess);
        
    end
    
    % plot these params
    subplot(223);
    h(iW) = plot(cfg_plot2.x+cfg_plot2.xoff(iW),this_data,'Color',cfg_plot2.col(iW),'Marker','.','MarkerSize',cfg_plot2.ms);
    hold on;
    errorbar(cfg_plot2.x+cfg_plot2.xoff(iW),this_data,this_dataSD,'Color',cfg_plot2.col(iW),'LineStyle','none');
    set(gca,'XLim',[0.5 3.5],'YLim',[0 55],'XTick',1:3,'XTickLabel',whatL,'FontSize',cfg_plot2.fs,'TickDir','out');
    box off;
    
end

ylabel('decoding error (cm)');
legend(h,whatP); legend boxoff;

%%
if cfg_plot.writeOutput
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300','decQuantify_multi2.png');
    print(gcf,'-depsc','decQuantify_multi2.eps');
end

%% normalized version -- should be able to use grand_av
clear grand_avN; figure;
for iW = 1:length(what)
    grand_avN{iW} = nan(size(grand_av{iW}));
    for iS = 2:size(grand_av{iW},1)
        grand_avN{iW}(iS,:,:) = sq(grand_av{iW}(iS,:,:))./sq(grand_av{1}(iS,:,:));
    end
end

for iW = 1:length(what)
    this_grand_av = sq(nanmean(grand_avN{iW},1));
    
    subplot(2,3,iW)
    nan_imagesc_ec(this_grand_av');
    set(gca,'XTick',1:size(this_grand_av',2),'XTickLabel',ALL_decErr.cfg.TCsmooth);
    set(gca,'YTick',1:size(this_grand_av',1),'YTickLabel',ALL_decErr.cfg.QsmoothSD,'TickDir','out');
    xlabel('TC smoothing (bins)'); ylabel('Q smoothing SD (s)');
    
    hold on; box off;
    
    [~,min_idx] = min(this_grand_av(:));
    [min_x,min_y] = ind2sub(size(this_grand_av),min_idx);
    
    hold on;
    plot(min_x,min_y,'*','Color',[1 1 1],'MarkerSize',10);
    
    title(sprintf('%s err %.2f',whatL{iW},this_grand_av(min_x,min_y)));
    %caxis([0 1]);
    colorbar;
    
end

%% compare for specific parameters
%figure;

%what = {[0 0],[0.1 0.002],[0.25 0.005],[0.5 0.01],[1 0.025] [5 0.2]}; % works for all sessions
what2 = {[0 0],[0.5 0.005],[1 0.01],[1 0.05] [3 0.1]}; % works for all sessions
whatL = {'same-lap','next-lap','leave-one-out'};
 
cfg_plot2 = [];
cfg_plot2.col = 'rgbcmk';
cfg_plot2.x = 1:3;
cfg_plot2.xoff = -0.1:0.05:0.15;
cfg_plot2.ms = 20;
cfg_plot2.fs = 10;

clear whatP this_data this_dataSD h
for iW = 1:length(what2)
    
    this_data = []; this_dataSD = [];
    
    % get data for these params
    whatP{iW} = sprintf('sTC %.2f, sQ %.3f',what2{iW}(1),what2{iW}(2));
    for iD = 1:3
        
        this_x = find(ALL_decErr.cfg.TCsmooth == what2{iW}(1));
        this_y = find(ALL_decErr.cfg.QsmoothSD == what2{iW}(2));
        
        this_data(iD) = nanmean(grand_avN{iD}(:,this_x,this_y));
        this_nSess = sum(~isnan((grand_avN{iD}(:,this_x,this_y))));
        this_dataSD(iD) = nanstd(grand_avN{iD}(:,this_x,this_y),[],1)/sqrt(this_nSess);
        
    end
    
    % plot these params
    subplot(224);
    h(iW) = plot(cfg_plot2.x+cfg_plot2.xoff(iW),this_data,'Color',cfg_plot2.col(iW),'Marker','.','MarkerSize',cfg_plot2.ms);
    hold on;
    errorbar(cfg_plot2.x+cfg_plot2.xoff(iW),this_data,this_dataSD,'Color',cfg_plot2.col(iW),'LineStyle','none');
    set(gca,'XLim',[0.5 3.5],'YLim',[0.75 5],'XTick',1:3,'XTickLabel',whatL,'FontSize',cfg_plot2.fs,'TickDir','out');
    box off;
    
end

ylabel('decoding error (normalized)');
legend(h,whatP,'Location','northwest'); legend boxoff;

%%
if cfg_plot.writeOutput
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpng','-r300','decQuantify_multi3.png');
    print(gcf,'-depsc','decQuantify_multi3.eps');
end