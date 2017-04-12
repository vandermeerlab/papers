function Naris_spec_fig(cfg_in, naris)
%%        :
%
%
%
%          Inputs:
%           -
%           -
%           -
%          Outputs:
%           -
%           -
%           -
%
% EC - 2017-01-08

cfg_def = [];
cfg_def.win = 512;
cfg_def.noverlap = cfg_def.win/4;
cfg_def.linewidth = 3;
cfg = ProcessConfig2(cfg_def, cfg_in);

%% append the data
all_naris.hdr.Fs = naris.pre.data_ampx.hdr.Fs;
[x,y] = size(naris.pre.data_ampx.channels{1});
if x>y
    all_naris.channels = {[naris.pre.data_ampx.channels{1}; naris.ipsi.data_ampx.channels{1};naris.contra.data_ampx.channels{1};naris.post.data_ampx.channels{1}]};
else
    all_naris.channels = {[naris.pre.data_ampx.channels{1}, naris.ipsi.data_ampx.channels{1},naris.contra.data_ampx.channels{1},naris.post.data_ampx.channels{1}]};
end
all_naris.tvec = [naris.pre.data_ampx.tvec', naris.ipsi.data_ampx.tvec',naris.contra.data_ampx.tvec',naris.post.data_ampx.tvec'];
all_naris.type = 'AMPX';
all_naris.labels = naris.pre.data_ampx.labels;
%% make a spectrogram of the combined data
[~,Fp,Tp,Pp] = spectrogram(naris.pre.data_ampx.channels{1},rectwin(cfg.win),cfg.noverlap,1:120,naris.pre.data_ampx.hdr.Fs);
[~,Fi,Ti,Pi] = spectrogram(naris.ipsi.data_ampx.channels{1},rectwin(cfg.win),cfg.noverlap,1:120,naris.ipsi.data_ampx.hdr.Fs);
[~,Fc,Tc,Pc] = spectrogram(naris.contra.data_ampx.channels{1},rectwin(cfg.win),cfg.noverlap,1:120,naris.contra.data_ampx.hdr.Fs);
 
[~,F,T,P] = spectrogram(all_naris.channels{1},rectwin(cfg.win),cfg.noverlap,1:120,all_naris.hdr.Fs);
% x_cor = corrcoef(P');
% imagesc(F,F,x_cor)
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');
if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R054')
    caxis([-50 80]);
else
    caxis([-150 -80])
end
h = vline([Tp(end),(Tp(end)+Ti(end)),(Tp(end)+Ti(end)+Tc(end))], {'k', 'k', 'k'} );
set(h(:), 'linewidth', cfg.linewidth)
% hold on
% if isfield(cfg, 'art')
%     for iArt = 1:length(cfg.art.tstart)
%         plot(cfg.art.tstart(iArt), 119, '*')
%     end
% end
SetFigure([],gcf)
set(gcf, 'position', [600 50 1600 780]);
saveas(gcf, 'All_phase_Spectrogram.fig')
saveas(gcf, 'All_phase_Spectrogram.png')
saveas(gcf, 'All_phase_Spectrogram', 'epsc')

% %% try the same thing but with FT
%
% data_ft = AMPX_makeft(all_naris);
% data_ft.label = data_ft.label{1};
%
% %%
% t_temp = data_ft.time{1} - data_ft.time{1}(1);
% data_ft.time{1} = t_temp;
% cfg = [];
% cfg.t = round(median([1 (data_ft.time{1}(end)-data_ft.time{1}(1))]));
% cfg.hdr = data_ft.hdr;
% cfg.toi = [data_ft.time{1}(1)/2000 data_ft.time{1}(end)/2000];
% trl = ft_maketrl(cfg)
%
% cfg = [];
% cfg.trl = trl;
% data_trl = ft_redefinetrial(cfg, data_ft)
%
% cfg              = []; % start with empty cfg
% cfg.output       = 'pow';
% cfg.trials     = 'all';
% cfg.method       = 'mtmconvol';
% cfg.taper        = 'hanning';
% cfg.foi          = 1:120; % frequencies of interest
% cfg.t_ftimwin    = ones(size(cfg.foi)).*0.5;  % window size: fixed at 0.5s
% % cfg.toi          = [data_ft.time{1}(1) data_ft.time{1}(end)];
%
% TFR = ft_freqanalysis(cfg, data_ft);
%
% figure
% cfg = [];
% ft_singleplotTFR(cfg, TFR);
end