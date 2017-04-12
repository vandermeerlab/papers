function Band_Plot_Labels
% used to add frequency band labels to PSDs
if exist('hline', 'file') ==0 || exist('vline', 'file') ==0
    warning('Band_Plot_Labels requires vline and hline.  Can be acquired from http://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline')
    axis tight; axis off;
    return
else
    y_ticks = [1.5 4 7 9 15 25 45 60 70 85 140 180]; % LFP ranges for delta-high Gamma
    vline(0, 'k')
    hline(y_ticks,{'y' 'y' 'g' 'g' 'r' 'r' 'b' 'b' 'c' 'c' 'm' 'm'})
    set(gca, 'XTick', [0])
    set(gca, 'YTick', median(reshape(y_ticks, 2, length(y_ticks)/2)));
    str ={'d' 't' 'b' 'g50' 'g80' 'U'};
    set(gca, 'YTickLabel', str, 'fontname', 'symbol')
end
end