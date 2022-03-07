function [N_Fig] = Naris_plot_mod(cfg_in, naris)
%% Naris_plot: plots the PSDs from each of the four naris phases.




cfg_def.linewidth = 1.5;
cfg_def.fontsize = 16;
cfg_def.Naris_phases = {'pre', 'ipsi', 'contra', 'post'};
cfg = ProcessConfig2(cfg_def, cfg_in);
%% plot the PSD
figHandles = get(0,'Children');
if isempty(figHandles) ==0 && sum(figHandles == 20) > 0
    close(20)
end
N_Fig = figure(20);
hold on
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
set(N_Fig, 'Position', [200, 200, 900 700])

% text(cfg.gamma(1)+3, -18, 'Gamma', 'Fontsize', 14)
% maximize
c_ord = linspecer(4);
for iPhase = 1:length(cfg.Naris_phases)
    for iCh = cfg.chan
        psd_norm = 10*log10(naris.(cfg.Naris_phases{iPhase}).data{iCh}.Data);
        plot(naris.(cfg.Naris_phases{iPhase}).data{iCh}.Frequencies,psd_norm,'color',c_ord(iPhase,:),'LineWidth',4);
        if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R054')
            set(gca,'XLim',[0 100],'XTick',0:10:100,'YLim',[-20 10],'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100},'YTick',[]); grid on;
            text(1, 45, num2str(naris.(cfg.Naris_phases{iPhase}).labels{iCh}), 'fontsize', cfg.fontsize)
        else
            set(gca,'XLim',[0 100],'XTick',0:10:100,'YLim',[-135 -115],'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100},'YTick',[]); grid on;
            text(1, 45, naris.(cfg.Naris_phases{iPhase}).labels{iCh}, 'fontsize', cfg.fontsize)
        end
        %         vline([45 60 75 85], {'g','g','b', 'b'}, {'Gamma50', ' ', 'Gamma80', ' '})
    end
end
y_lim = get(gca, 'ylim');
y_lim = y_lim(1)+0.01;
rectangle('position', [cfg.gamma(1,1), y_lim, (cfg.gamma(1,2)-cfg.gamma(1,1)), 1.4], 'FaceColor',[0.7 0.7 0.7], 'edgecolor', [1 1 1])
rectangle('position', [cfg.gamma(2,1), y_lim, (cfg.gamma(2,2)-cfg.gamma(2,1)), 1.4], 'FaceColor',[0.8 0.8 0.8], 'edgecolor', [1 1 1])

% legend(cfg.Naris_exp)
xlabel('Frequency (Hz)', 'fontsize', 16)%, 'fontweight', 'bold');
ylabel('Power', 'fontsize', 16)%, 'fontweight', 'bold')
box on
SetFigure([], gcf)

%% save the plots. and the data.

if ispc
    %save data.
    mkdir([cd '\' date])
    saveas(N_Fig, [cd '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.fig'])
    print(N_Fig, '-dpng','-r300',[cd '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
    % save it in the all Naris Folder
    mkdir([cfg.data_path '\'  date])
    saveas(N_Fig, [cfg.data_path  '\Naris_' cfg.fname '_comp.fig'])
    print(N_Fig, '-dpng','-r300',[cfg.data_path '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
else
    mkdir([cd '/' date])
    saveas(N_Fig, [cd  '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.fig'])
    print(N_Fig, '-dpng','-r300',[cd '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
    % save it in the all Naris Folder
    mkdir([cfg.data_path '/'  date])
    saveas(N_Fig, [cfg.data_path '/Naris_' cfg.fname '_comp.fig'])
    print(N_Fig, '-dpng','-r300',[cfg.data_path '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
end