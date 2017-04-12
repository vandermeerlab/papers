%%
clear all;
cd('D:\My_Documents\Dropbox\projects\DecodingQuantified\results');

% output files to use -- from GENERATE_loocv_byLap.m
what = {'ALL_decErr_LOOCV_byLap_aTnbT_0_0.mat', ...
  'ALL_decErr_LOOCV_byLap_aTnbT_0005_05.mat', ...
    'ALL_decErr_LOOCV_byLap_aTnbT_001_1.mat', ...
    'ALL_decErr_LOOCV_byLap_aTnbT_005_1.mat', ...
    'ALL_decErr_LOOCV_byLap_aTnbT_01_3.mat',...
    };

cfg_plot = [];
cfg_plot.xbin = 1:10;
cfg_plot.ylim1 = [0 50];
cfg_plot.ylim2 = [0 1.2];
cfg_plot.fs = 18; % font
cfg_plot.ms = 20; % markersize
cfg_plot.lw = 1;
cfg_plot.cols = 'rgbcmk';
cfg_plot.N = 24; % N to use for errorbar SEM; note this gets overwritten based on actual number of sessions used
cfg_plot.offset = -0.15:0.05:0.15;
cfg_plot.maxPnan = 0.2; % maximum proportion of decoding time bins that can be NaN for session to be included
cfg_plot.minPassed = 0.8; % minimum proportion of nActiveNeurons criterion met for session to be included
cfg_plot.minCells = 20;

nP = length(what);

%%
for iP = 1:nP
    
    % load data
    load(what{iP});
    
    keep_maskL = ones(size(ALL_decErr.left.meanErr)); % track how many sessions survive selection
    keep_maskR = ones(size(ALL_decErr.right.meanErr));
    
    this_meanErrL = ALL_decErr.left.meanErr;
    this_meanErrR = ALL_decErr.right.meanErr;
    
    % mask
    toss_idx = find(ALL_decErr.left.Pnan > cfg_plot.maxPnan | ALL_decErr.left.Ppassed < cfg_plot.minPassed);
    keep_maskL(toss_idx) = 0;
    
    toss_idx = find(ALL_decErr.right.Pnan > cfg_plot.maxPnan | ALL_decErr.right.Ppassed < cfg_plot.minPassed);
    keep_maskR(toss_idx) = 0;
    
    % mask - nCells
    nSessions = size(ALL_decErr.left.meanErr,1);
    keep_idx = ones(nSessions,1);
    for iS = 1:nSessions
        if length(ALL_decErr.left.S{iS}.t) < cfg_plot.minCells | length(ALL_decErr.right.S{iS}.t) < cfg_plot.minCells
            keep_maskL(iS,:) = 0;
            keep_maskR(iS,:) = 0;
        end
    end
    
    this_meanErrL(~keep_maskL) = NaN;
    this_meanErrR(~keep_maskR) = NaN;
    
    fprintf('nSessions for param %d:\n',iP);
    cfg_plot.N = sum(keep_maskL | keep_maskR);
    disp(cfg_plot.N);
    
    clear temp;
    temp(1,:,:) = this_meanErrL; temp(2,:,:) = this_meanErrR;
    both_meanErr = sq(nanmean(temp));
    
    both_meanErrN = both_meanErr./repmat(both_meanErr(:,1),[1 size(both_meanErr,2)]);
    
    %%
    subplot(121);
    this_mean = nanmean(both_meanErr); this_sem = nanstd(both_meanErr)./sqrt(cfg_plot.N);
    
    h1(iP) = plot(cfg_plot.xbin+cfg_plot.offset(iP),this_mean(1:length(cfg_plot.xbin)),'Color', ...
        cfg_plot.cols(iP),'MarkerSize',cfg_plot.ms,'Marker','.');
    hold on;
    errorbar(cfg_plot.xbin+cfg_plot.offset(iP),this_mean(1:length(cfg_plot.xbin)),...
        this_sem(1:length(cfg_plot.xbin)),'Color',cfg_plot.cols(iP),'Marker','none');
    
    set(gca,'LineWidth',cfg_plot.lw,'XLim',[cfg_plot.xbin(1)-0.5 cfg_plot.xbin(end)+0.5], ...
        'FontSize',cfg_plot.fs,'TickDir','out','YLim',cfg_plot.ylim1);
    box off;
    xlabel('number of decoding laps'); ylabel('decoding error (cm)');
    
    h1t{iP} = sprintf('sTC %.2f, sQ %.3f',ALL_decErr.cfg.TCsmooth,ALL_decErr.cfg.QsmoothSD);
    
    %
    subplot(122);
    this_mean = nanmean(both_meanErrN); this_sem = nanstd(both_meanErrN)./sqrt(cfg_plot.N);
    
    h2(iP) = plot(cfg_plot.xbin+cfg_plot.offset(iP),this_mean(1:length(cfg_plot.xbin)),'Color', ...
        cfg_plot.cols(iP),'MarkerSize',cfg_plot.ms,'Marker','.');
    hold on;
    errorbar(cfg_plot.xbin+cfg_plot.offset(iP),this_mean(1:length(cfg_plot.xbin)),...
        this_sem(1:length(cfg_plot.xbin)),'Color',cfg_plot.cols(iP),'Marker','none');
    
    set(gca,'LineWidth',cfg_plot.lw,'XLim',[cfg_plot.xbin(1)-0.5 cfg_plot.xbin(end)+0.5], ...
        'FontSize',cfg_plot.fs,'TickDir','out','YLim',cfg_plot.ylim2);
    box off;
    xlabel('number of decoding laps'); ylabel('decoding error (normalized)');
    
    h2t{iP} = sprintf('sTC %.2f, sQ %.3f',ALL_decErr.cfg.TCsmooth,ALL_decErr.cfg.QsmoothSD);
    
end

%% legends
%subplot(121);
%legend(h1,h1t,'Location','southwest'); legend boxoff;

subplot(122);
legend(h2,h2t,'Location','southwest'); legend boxoff;

%%
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','-r300','decAcc_byLap_allTrials_multi.png');
print(gcf,'-depsc','decAcc_byLap_alLTrials_multi.eps');