function [naris, cfg] = Naris_fast(cfg_in)
%Naris_Fast  Loads the data for s single set of channels (normally a
%tetrode) and then pltos all four session PSDs in a 2x2 plot as well as a
%overlapping PSD plot.  Will save different file extensions for different
%hanning window sizes.
%
%Inputs:
%   cfg [struct]: can be blank if you want to use the defaults.
%% default parameters
cfg_def.chan = 1; % used to pick the first channel in the tetrode. left over from an older version
cfg_def.hann_win_fac = 4;
cfg_def.hann_win = 1024*cfg_def.hann_win_fac;
cfg_def.df = 10;
cfg_def.whitefilter = 'on';
cfg_def.gamma= [40 55; 70 85];
cfg  = ProcessConfig2(cfg_def, cfg_in);
cfg.data_path = '/Users/jericcarmichael/Documents/Nairs_data';
%%
cfg.fname = mkfile;
[cfg] = Naris_cfgs(cfg);
cfg
warning off

if strcmp(cfg.fname(1:4), 'R053')
    cfg.Naris_exp = {'pre', 'right', 'left', 'post'};
    cfg.method = 'raw';
    cfg.psd = 'normal';
elseif strcmp(cfg.fname(1:4), 'R060')
    cfg.Naris_exp = {'pre', 'right', 'left', 'post'};
    cfg.method = 'ratio';
    cfg.psd = 'white';
else
    cfg.Naris_exp = {'pre', 'ipsi', 'contra', 'post'};
    cfg.method = 'ratio';
    cfg.psd = 'white';
end

if ~isfield(cfg, 'tetrodes')
    if ~exist([cfg.fname '-pre_data.mat'])
        data = AMPX_loadData([cfg.fname '-pre.dat'],(1:64), cfg.df);
        save([cfg.fname '-pre_data.mat'], 'data');
    else
        load([cfg.fname '-pre_data.mat']);
    end
    LoadExpKeys
    % detect the channel with the highest gamma power.  only do this for the pre to keep it consistent across phases.
    cfg.ch = 1:64;
    ch_idx = cfg.ch(ExpKeys.BadChannels);
    cfg.ch(ch_idx) = [];
    cfg.tetrodes(1) = AMPX_detect_best_chan(cfg, data, ExpKeys);
    
    
    %     clear data
end

% cfg.notch = 0;
%%
for exp = 1:4
    cfg.naris_type = cfg.Naris_exp{exp};
    
    data = Naris_psd(cfg);
    disp([data.hdr.Filename(20:end-4) ' Data has been loaded'])
    %% plot the data
    cfg.sub_num = exp;
    
    n_fig = Naris_plot(data, cfg);
    
    naris.(cfg.Naris_exp{exp}).data = data.psd;
    for ich = 1:length(data.labels)
        naris.(cfg.Naris_exp{exp}).labels{ich} = data.labels;
        naris.(cfg.Naris_exp{exp}).cfg = cfg;
    end
end

%%
ah = axes('position',[0,0,1,1],'visible','off');
line([0,1],[0.5,.5],'parent',ah,'linewidth',1.5, 'color', 'k');
line([0.5,.5],[0,1],'parent',ah,'linewidth',1.5, 'color', 'k');
text(.255, .98, cfg.Naris_exp{1}, 'Parent', ah, 'fontsize', 28);
text(.7, .98, cfg.Naris_exp{2}, 'Parent', ah, 'fontsize', 28)
text(.255, .025, cfg.Naris_exp{3}, 'Parent', ah, 'fontsize', 28)
text(.705, .025, cfg.Naris_exp{4}, 'Parent', ah, 'fontsize', 28)
% if isfield(data.hdr, 'filter') && strcmp(data.hdr.filter, 'notch')
%     text(.05, .98, [cfg.fname 'w/notch'], 'Parent', ah, 'fontsize', 18)
% else
text(.05, .98, cfg.fname, 'Parent', ah, 'fontsize', 18)
% end
%% save the plots. and the data.
if isfield(cfg, 'whitefilter');
    save(['White_Naris_data_' num2str(cfg.hann_win_fac) '_new'], 'naris', '-v7.3')
else
    save(['Naris_data_' num2str(cfg.hann_win_fac)], 'naris', '-v7.3')
end
%save data.

if isfield(cfg, 'whitefilter'); cfg.fname = [cfg.fname '_white']; end

if ispc
    mkdir([cd '\' date])
    saveas(gcf, [cd '\' date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.fig'])
    print(gcf, '-dpng','-r300',[cd '\' date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.png'])
    
    % save it in the all Naris Folder
    mkdir(['G:\Naris\figures\'  date])
    saveas(gcf, ['G:\Naris\figures\'  date '\Naris_' cfg.fname '.fig'])
    print(gcf, '-dpng','-r300',['G:\Naris\figures\'  date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.png'])
    
else
    mkdir([cd '/' date])
    saveas(gcf, [cd '/' date '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.fig'])
    print(gcf, '-dpng','-r300',[cd '/' date '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.png'])
    % save it in the all Naris Folder
    mkdir([cfg.data_path '/figures/'  date])
    saveas(gcf, [cfg.data_path '/figures/' date '/Naris_' cfg.fname '.fig'])
    print(gcf, '-dpng','-r300',[cfg.data_path '/figures/'  date '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '.png'])
end


%%
figHandles = get(0,'Children');
if isempty(figHandles) ==0 && sum(figHandles == 20) > 0
    close(20)
end
F = figure(20);
hold on
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
set(F, 'Position', [200, 200, 900 700])

% text(cfg.gamma(1)+3, -18, 'Gamma', 'Fontsize', 14)

% maximize
c_ord = linspecer(4);
for exp = 1:4
    for iCh = cfg.chan
        psd_norm = 10*log10(naris.(cfg.Naris_exp{exp}).data{iCh}.Data);
        plot(naris.(cfg.Naris_exp{exp}).data{iCh}.Frequencies,psd_norm,'color',c_ord(exp,:),'LineWidth',4);
        set(gca,'XLim',[0 100],'XTick',0:10:100,'YLim',[-20 10],'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100},'YTick',[0 30]); grid on;
        %         vline([45 60 75 85], {'g','g','b', 'b'}, {'Gamma50', ' ', 'Gamma80', ' '})
        text(1, 45, num2str(data.labels(iCh)), 'fontsize', 16)
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
    saveas(F, [cd '\' date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.fig'])
    print(F, '-dpng','-r300',[cd '\' date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
    % save it in the all Naris Folder
    mkdir([cfg.data_path '\'  date])
    saveas(F, [cfg.data_path '\' date '\Naris_' cfg.fname '_comp.fig'])
    print(F, '-dpng','-r300',[cfg.data_path '\  date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
else
    mkdir([cd '/' date])
    saveas(F, [cd '/' date '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.fig'])
    print(F, '-dpng','-r300',[cd '/' date '\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
    % save it in the all Naris Folder
    mkdir([cfg.data_path '/'  date])
    saveas(F, [cfg.data_path '/'  date '/Naris_' cfg.fname '_comp.fig'])
    print(F, '-dpng','-r300',[cfg.data_path '/' date '/Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])
    
end


%% save the best ones:
%
% saveas(F, ['G:\Naris\paper_figs\' 'Naris_' cfg.fname '_comp.fig'])
% print(F, '-dpng','-r300',['G:\Naris\paper_figs\Naris_' num2str(cfg.hann_win_fac) '_' cfg.fname '_comp.png'])


