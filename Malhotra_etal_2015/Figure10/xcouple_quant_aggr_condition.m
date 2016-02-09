% added zscore of power

clear all;
rng('shuffle');
%addpath('C:\Users\mvdm\Dropbox\projects\Sushant\shared');

fd = {'C:\data\R014-2012-09-14_promoted', ...
    'C:\data\R014-2012-09-16_promoted', ...
    'C:\data\R016-2012-10-05_promoted', ...
    'C:\data\R016-2012-10-08_promoted', ...
    'C:\data\R018-2012-12-08_promoted', ...
    'C:\data\R018-2012-12-14_promoted'};

%cd('D:\My_Documents\My Dropbox\projects\Sushant\2014-09-19') % athena
%load fd; % list of files to process, obtained from remove_cueviolation.m

cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
load fd_isidro

% only certain rat
%fd = fd(end-2:end); % R014
%fd = fd(1); % R020
%fd = fd(2:8); % R018
%fd = fd(9:13); % R016

%%
%twin = [0 1.25];
twin = [-1 3];
f_lg = [50 65];
f_hg = [70 85];

f_x = [8 12];

doXCorr = 0;

% first, build dataset
data_all = []; event_all = [];
for iFD = 1:length(fd)
    
   cd(fd{iFD});
   fprintf('\n***MAIN: Entering folder %s (%d/%d)\n\n',pwd,iFD,length(fd));
   
    run(FindFile('*keys.m'));
    fc = ExpKeys.goodGamma(1);
    data = ft_read_neuralynx_interp(fc);
    
    fprintf('\n***MAIN: Finished loading data\n\n');
 
    data.label{1} = 'vStr';
    
    %% get hilberted power -- filters could be improved!
    data2 = data;
    
    % low-gamma
    data2f = ft_filterLFP(data2,f_lg,'ford',4);
    data.trial{1}(2,:) = data2f.trial{1};
    data.label{2} = 'lg_pow';
    data.trial{1}(2,:) = abs(hilbert(data.trial{1}(2,:))).^2;
    data.trial{1}(2,:) = zscore(data.trial{1}(2,:));
    
    % high-gamma
    data2f = ft_filterLFP(data2,f_hg,'ford',4);
    data.trial{1}(3,:) = data2f.trial{1};
    data.label{3} = 'hg_pow';
    data.trial{1}(3,:) = abs(hilbert(data.trial{1}(3,:))).^2;
    data.trial{1}(3,:) = zscore(data.trial{1}(3,:));
    
    % cross-frequency
    d = fdesign.bandpass('N,F3dB1,F3dB2',4,f_x(1),f_x(2));
    Hd = design(d,'butter');
    
    rmpath('C:\Users\mvdm\Documents\GitHub\fieldtrip\external\signal'); % use correct filtfilt()
    rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal');
    rmpath('C:\Users\mvdm\GitHub\fieldtrip\external\signal'); % kronos
    
    data2f = ft_filterLFP(data2,f_x,'B',Hd.sosMatrix,'A',Hd.scaleValues);
    data.trial{1}(4,:) = data2f.trial{1};
    data.label{4} = 'x_phase';
    data.trial{1}(4,:) = angle(hilbert(data.trial{1}(4,:)));
    
    %% get all nosepokes
    
    all_cues = {'c1','c3','c5','lo','hi'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
    
    for iCue = 1:length(all_cues)
        
        fprintf('\n***MAIN: Obtaining cue set (%d/%d)\n\n',iCue,length(all_cues));
        
        cfg = [];
        cfg.trialdef.cue = all_cues(iCue); 
        cfg.trialdef.block = 'both';
        cfg.trialdef.eventtype = 'nosepoke';
        cfg.trialdef.location = 'both';
        cfg.trialfun = 'ft_trialfun_lineartracktone2';
        cfg.trialdef.hdr = data.hdr;
        cfg.trialdef.pre = 2; cfg.trialdef.post = 4;
        
        trl = ft_trialfun_lineartracktone2(cfg); cfg.trl = trl;
        this_data = ft_redefinetrial(cfg,data);
        
        event = repmat(iCue,[length(this_data.trial) 1]);
        
        if isempty(data_all)
            data_all = this_data;
            event_all = event;
        else
            fprintf('\n***MAIN: Appending data\n\n',iCue,length(all_cues));
            data_all = ft_appenddata([],data_all,this_data);
            event_all = cat(1,event_all,event);
        end
        
    end
    
end

%% analysis
dpi = pi/18;
pbin_edges = -pi:dpi:pi; pbin_centers = pbin_edges(1:end-1)+dpi/2;
maxlag = 2000; % samples for xcorr
nShuf = 1000;

temp = data_all;

% restrict data
cfg = []; cfg.latency = twin;
temp = ft_selectdata(cfg,temp);
trial_len = length(temp.trial{1}(1,:));

shuf_xf_lg = {}; shuf_xf_hg = {}; obs_xf_lg = {}; obs_xf_hg = {};
   obs_xc = {}; shuf_xc = {};

% do analysis for each trial type separately
trial_types = [1 3];
for iTT = 1:length(trial_types)
    
    cur_tt_idx = find(event_all == trial_types(iTT));
    
    % first get observed PAC value
    
    trial_mat = cell2mat(temp.trial(cur_tt_idx)); % trials concatenated
    
    lg_idx = strmatch('lg_pow',temp.label);
    obs_lg_pwr = trial_mat(lg_idx,:); % note, ft_selectdata switches labels!!
    
    hg_idx = strmatch('hg_pow',temp.label);
    obs_hg_pwr = trial_mat(hg_idx,:); % note, ft_selectdata switches labels!!
    
    xf_idx = strmatch('x_phase',temp.label);
    obs_xf_phase = trial_mat(xf_idx,:); % note, ft_selectdata switches labels!!
    
    obs_pac_lg = abs(mean(obs_lg_pwr.*exp(1i*obs_xf_phase)));
    obs_pac_hg = abs(mean(obs_hg_pwr.*exp(1i*obs_xf_phase)));
    
    % get observed phase-frequency relationships
    %xf1 = averageXbyYbin(tiedrank(obs_lg_pwr),obs_xf_phase,pbin_edges);
    %xf2 = averageXbyYbin(tiedrank(obs_hg_pwr),obs_xf_phase,pbin_edges);
    
    obs_xf_lg{iTT} = averageXbyYbin(obs_lg_pwr,obs_xf_phase,pbin_edges);
    obs_xf_hg{iTT} = averageXbyYbin(obs_hg_pwr,obs_xf_phase,pbin_edges);
    
    % improved xcorr that does not go across trials
    temp2 = temp; temp2.trial = temp2.trial(cur_tt_idx);
    if doXCorr
        for iT = length(temp2.trial):-1:1
            [obs_xc{iTT}(iT,:),lags] = xcov(tiedrank(temp2.trial{iT}(hg_idx,:)),tiedrank(temp2.trial{iT}(lg_idx,:)),maxlag,'coeff');
        end
        obs_xc{iTT} = nanmean(obs_xc{iTT});
        shuf_xc{iTT} = zeros(nShuf,length(lags));
    end
    
    % shuffle
    shuf_pac_lg = zeros(nShuf,1); shuf_pac_hg = zeros(nShuf,1);
    %shuf_xf_lg = zeros(nShuf,length(pbin_edges)-1); shuf_xf_hg = zeros(nShuf,length(pbin_edges)-1);
    
    for iShuf = 1:nShuf
        shuf_data = temp2;
        
        % random amount to shift data
        perm_shift1 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
        perm_shift2 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
        
        % for each trial, shift the data
        for iT = 1:length(shuf_data.trial)
            
            piece1 = shuf_data.trial{iT}(lg_idx,perm_shift1+1:end); piece2 = shuf_data.trial{iT}(lg_idx,1:perm_shift1);
            shuf_data.trial{iT}(lg_idx,1:trial_len-perm_shift1) = piece1;
            shuf_data.trial{iT}(lg_idx,end-perm_shift1+1:end) = piece2;
            
            piece1 = shuf_data.trial{iT}(hg_idx,perm_shift2+1:end); piece2 = shuf_data.trial{iT}(hg_idx,1:perm_shift2);
            shuf_data.trial{iT}(hg_idx,1:trial_len-perm_shift2) = piece1;
            shuf_data.trial{iT}(hg_idx,end-perm_shift2+1:end) = piece2;
            
            
        end
        
        % get PAC for this shuffle
        trial_mat = cell2mat(shuf_data.trial); % trials concatenated
        
        shuf_lg_pwr = trial_mat(lg_idx,:);
        shuf_hg_pwr = trial_mat(hg_idx,:);
        
        shuf_xf_lg{iTT}(iShuf,:) = averageXbyYbin(shuf_lg_pwr,obs_xf_phase,pbin_edges);
        shuf_xf_hg{iTT}(iShuf,:) = averageXbyYbin(shuf_hg_pwr,obs_xf_phase,pbin_edges);
        
        shuf_pac_lg(iShuf) = abs(mean(shuf_lg_pwr.*exp(1i*obs_xf_phase)));
        shuf_pac_hg(iShuf) = abs(mean(shuf_hg_pwr.*exp(1i*obs_xf_phase)));
        
        
        % get xcorr for this shuffle
        temp_shuf = [];
        if doXCorr
            for iT = length(shuf_data.trial):-1:1
                [temp_shuf(iT,:),~] = xcov(tiedrank(shuf_data.trial{iT}(hg_idx,:)),tiedrank(shuf_data.trial{iT}(lg_idx,:)),maxlag,'coeff');
            end
            shuf_xc{iTT}(iShuf,:) = nanmean(temp_shuf);
            
        end
        
        % compute some PAC stats
        pacz_lg{iTT} = (obs_pac_lg-mean(shuf_pac_lg))/std(shuf_pac_lg);
        pacz_hg{iTT} = (obs_pac_hg-mean(shuf_pac_hg))/std(shuf_pac_hg);
        
        pacm_lg{iTT} = 1-(sum(shuf_pac_lg < obs_pac_lg)/nShuf);
        pacm_hg{iTT} = 1-(sum(shuf_pac_hg < obs_pac_hg)/nShuf);
        
    end % shuffles
    
end % loop over trial types

%% plot
for iTT = 1:length(trial_types)

    
figure(iTT);

if doXCorr
    %lags = lags.*1/data.fsample;
    
    subplot(221);
    plot(lags,obs_xc{iTT},'LineWidth',2)
    hold on
    plot(lags,nanmean(shuf_xc{iTT}),'k');
    plot(lags,nanmean(shuf_xc{iTT})+nanstd(shuf_xc{iTT}),'k:');
    plot(lags,nanmean(shuf_xc{iTT})-nanstd(shuf_xc{iTT}),'k:');
    
    plot([0 0],get(gca,'YLim'),'--','Color',[0.7 0.7 0.7]);
    

    set(gca,'FontSize',18,'YLim',[-0.05 0.15]);
    title('xcorr low-high gamma');
    xlabel('time (s)'); ylabel('power rank correlation');
end

subplot(222);
plot(pbin_centers,obs_xf_lg{iTT},'LineWidth',2)

hold on
plot(pbin_centers,nanmean(shuf_xf_lg{iTT}),'k');
plot(pbin_centers,nanmean(shuf_xf_lg{iTT})+nanstd(shuf_xf_lg{iTT}),'k:');
plot(pbin_centers,nanmean(shuf_xf_lg{iTT})-nanstd(shuf_xf_lg{iTT}),'k:');

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('lg pacz %.2f (p = %.3f)',pacz_lg{iTT},1-normcdf(pacz_lg{iTT},0,1)));


subplot(224);
plot(pbin_centers,obs_xf_hg{iTT},'LineWidth',2)

hold on
plot(pbin_centers,nanmean(shuf_xf_hg{iTT}),'k');
plot(pbin_centers,nanmean(shuf_xf_hg{iTT})+nanstd(shuf_xf_hg{iTT}),'k:');
plot(pbin_centers,nanmean(shuf_xf_hg{iTT})-nanstd(shuf_xf_hg{iTT}),'k:');

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('hg pacz %.2f (p = %.3f)',pacz_hg{iTT},1-normcdf(pacz_hg{iTT},0,1)));

subplot(223);
title(sprintf('ALL(%d) %d-%d Hz, %d sess %d trls twin %.2f-%.2f',trial_types(iTT),f_x(1),f_x(2),length(fd),sum(event_all == trial_types(iTT)),twin(1),twin(2)));
axis off;

end

%% together
figure(3);
cols = 'cm';
clear h;

%
if doXCorr
for iC = 1:2
    
    subplot(221);
    h(iC) = plot(lags.*1/data.fsample,obs_xc{iC},'LineWidth',2,'Color',cols(iC));
    hold on
    
end

plot([0 0],get(gca,'YLim'),'--','Color',[0.7 0.7 0.7]);

set(gca,'FontSize',18,'YLim',[-0.04 0.06]);
title('xcorr low-high gamma');
xlabel('time (s)'); ylabel('power rank correlation');
legend(h,{'1p','5p'}); legend boxoff;
end
%
clear h;
for iC = 1:2
    subplot(222);
    h(iC) = plot(pbin_centers,obs_xf_lg{iC},'LineWidth',2,'Color',cols(iC));
    hold on;
    
end

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('low-gamma'));
legend(h,{'1p','5p'}); legend boxoff;

%
clear h;
for iC = 1:2
    subplot(224);
    h(iC) = plot(pbin_centers,obs_xf_hg{iC},'LineWidth',2,'Color',cols(iC));
    hold on;
    
end

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('high-gamma'));
legend(h,{'1p','5p'}); legend boxoff;

%% write figures
figure(1)
fn = 'ALL_xcouple_quant_1p-8-12Hz';
saveas(gcf,[fn '.fig']);
print(gcf,'-dill',[fn '.ai']);
print(gcf,'-dpng','-r300',[fn '.png']);

figure(2)
fn = 'ALL_xcouple_quant_5p-8-12Hz';
saveas(gcf,[fn '.fig']);
print(gcf,'-dill',[fn '.ai']);
print(gcf,'-dpng','-r300',[fn '.png']);

figure(3)
fn = 'ALL_xcouple_quant_1vs5p-8-12Hz';
saveas(gcf,[fn '.fig']);
print(gcf,'-dill',[fn '.ai']);
print(gcf,'-dpng','-r300',[fn '.png']);
