function [N_fig] = Naris_plot(data, cfg )
%UNTITLED8 Summary of this function goes here
N_fig = figure(10);
maximize
for iCh = cfg.chan
    subtightplot(2,2,cfg.sub_num);
    
    psd_norm = 10*log10(data.psd{iCh}.Data);
    plot(data.psd{iCh}.Frequencies,psd_norm,'b','LineWidth',2);
    set(gca,'XLim',[0 100],'XTick',0:10:100,'YLim',[-20 10],'XTickLabel',{},'YTick',[]); grid on;
    vline([45 60 75 85], {'g','g','b', 'b'}, {'Gamma50', ' ', 'Gamma80', ' '})
    text(1, 45, num2str(data.labels(iCh)), 'fontsize', 16)
   
end

