clear all;
rng('shuffle');
%addpath('C:\Users\mvdm\Dropbox\projects\Sushant\shared');

fd = {'C:\data\R014\R014-2012-09-14_promoted', ...
    'C:\data\R014\R014-2012-09-16_promoted', ...
    'C:\data\R016\R016-2012-10-05_promoted', ...
    'C:\data\R016\R016-2012-10-08_promoted', ...
    'C:\data\R018\R018-2012-12-08_promoted', ...
    'C:\data\R018\R018-2012-12-14_promoted'};

%cd('D:\My_Documents\My Dropbox\projects\Sushant\2014-09-19') % athena
cd('D:\My_Documents\Dropbox\projects\Sushant\2014-09-19'); % equinox
load fd; % list of files to process, obtained from remove_cueviolation.m

% only certain rat
%fd = fd(end-2:end); % R014
%files = fd(1); % R020
%fd = fd(2:8); % R018
%fd = fd(9:13); % R016

%%
twin = [-1 3];
f_lg = [50 65];
f_hg = [70 85];

f_x = [3 5];

doXCorr = 0;

% first, build dataset
data_all = []; event_all = [];
for iFD = 1:length(fd)
    
   cd(fd{iFD});

    run(FindFile('*keys.m'));
    fc = ExpKeys.goodGamma(1);
    data = ft_read_neuralynx_interp(fc); % this doesn't give NaNs...
 
    data.label{1} = 'vStr';
    
    %% get hilberted power -- filters could be improved!
    data2 = data;
    
    % low-gamma
    [data2f,nan_idx] = ft_filterLFP(data2,f_lg,'ford',4,'nan_mode','zero');
    data.trial{1}(2,:) = data2f.trial{1};
    data.label{2} = 'lg_pow';
    data.trial{1}(2,:) = abs(hilbert(data.trial{1}(2,:))).^2;
    data.trial{1}(2,:) = zscore(data.trial{1}(2,:));
    data.trial{1}(2,nan_idx) = nan;
    
    % high-gamma
    data2f = ft_filterLFP(data2,f_hg,'ford',4,'nan_mode','zero');
    data.trial{1}(3,:) = data2f.trial{1};
    data.label{3} = 'hg_pow';
    data.trial{1}(3,:) = abs(hilbert(data.trial{1}(3,:))).^2;
    data.trial{1}(3,:) = zscore(data.trial{1}(3,:));
    data.trial{1}(3,nan_idx) = nan;
    
    % cross-frequency
    d = fdesign.bandpass('N,F3dB1,F3dB2',4,f_x(1),f_x(2),2e3);
    Hd = design(d,'butter');
    
    rmpath('C:\Users\mvdm\Documents\GitHub\fieldtrip\external\signal'); % use correct filtfilt()
    rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal');
    rmpath('C:\Users\mvdm\GitHub\fieldtrip\external\signal');
    
    data2f = ft_filterLFP(data2,f_x,'B',Hd.sosMatrix,'A',Hd.scaleValues,'nan_mode','zero');
    data.trial{1}(4,:) = data2f.trial{1};
    data.label{4} = 'x_phase';
    data.trial{1}(4,:) = angle(hilbert(data.trial{1}(4,:)));
    data.trial{1}(4,nan_idx) = nan;
    
    %% make some trials
    
    cfg.twin = [-2 4];
    
    if ~isempty(ExpKeys.Prerecord)
        pre_t = ExpKeys.Prerecord(1)-cfg.twin(1):diff(cfg.twin):ExpKeys.Prerecord(2)-cfg.twin(2);
    else
        pre_t = [];
    end
    
    post_t = ExpKeys.Postrecord(1)-cfg.twin(1):diff(cfg.twin):ExpKeys.Postrecord(2)-cfg.twin(2);
    
    cfg = [];
    cfg.t = cat(2,pre_t,post_t);
    
    % select specific t's with low speed
    load(FindFile('*vt.mat'))
    s = GetLinSpd(x,y);
    t_idx = nearest_idx3(cfg.t,Range(s));
    sData = Data(s); sData = smooth(sData,15);
    keep_idx = sData(t_idx) < 5;
    cfg.t = cfg.t(keep_idx);
    
    cfg.mode = 'nlx';
    cfg.hdr = data.hdr;
    
    trl = ft_maketrl(cfg);
    
    cfg = []; cfg.trl = trl;
    this_data = ft_redefinetrial(cfg,data);
    
    mf = @(x) ~any(isnan(x(:))); % could relax this to a few NaNs
    keep_idx = cellfun(mf,this_data.trial);
    fprintf('\n*This session has %d out of %d trials without NaNs...\n\n',sum(keep_idx),length(keep_idx));
    
    if isempty(data_all)
        data_all = this_data;
        %event_all = event;
    else
        data_all = ft_appenddata([],data_all,this_data);
        %event_all = cat(1,event_all,event);
    end
    
end

%% remove NaN trials
mf = @(x) ~any(isnan(x(:))); % could relax this to a few NaNs
keep_idx = cellfun(mf,data_all.trial);
fprintf('Keeping %d out of %d trials...\n',sum(keep_idx),length(keep_idx));

cfg = []; cfg.trials = find(keep_idx);
data_all = ft_selectdata(cfg,data_all);

% randomize trials to enable faster xcorr computation later
%data_all.trial = data_all.trial(randperm(length(data_all.trial)));

%% check phase distribution
figure;
clear ht
for iTrial = 1:length(data_all.trial)
ht(iTrial,:) = hist(data_all.trial{iTrial}(4,:),100);
end
plot(nanmean(ht))

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

% first get observed PAC value
trial_mat = cell2mat(temp.trial); % trials concatenated

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

obs_xf_lg = averageXbyYbin(obs_lg_pwr,obs_xf_phase,pbin_edges);
obs_xf_hg = averageXbyYbin(obs_hg_pwr,obs_xf_phase,pbin_edges);

% improved xcorr that does not go across trials
obs_xc = [];
if doXCorr
    for iT = length(temp.trial):-1:1
        [obs_xc(iT,:),lags] = xcov(tiedrank(temp.trial{iT}(hg_idx,:)),tiedrank(temp.trial{iT}(lg_idx,:)),maxlag,'coeff');
    end
    obs_xc = nanmean(obs_xc);
    shuf_xc = zeros(nShuf,length(lags));
end

% shuffle
shuf_pac_lg = zeros(nShuf,1); shuf_pac_hg = zeros(nShuf,1);
shuf_xf_lg = zeros(nShuf,length(pbin_edges)-1); shuf_xf_hg = zeros(nShuf,length(pbin_edges)-1);

for iShuf = 1:nShuf
    shuf_data = temp;
    
    % random amount to shift data
    perm_shift1 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
    perm_shift2 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
    
    % for each trial, shift the data
    for iT = 1:length(shuf_data.trial)
        
        %shuf_data.trial{iT}(lg_idx,:) = circshift(shuf_data.trial{iT}(lg_idx,:),[0 perm_shift1]);
        %shuf_data.trial{iT}(hg_idx,:) = circshift(shuf_data.trial{iT}(hg_idx,:),[0 perm_shift2]);
        
        %shuf_data.trial{iT}(lg_idx,:) = [shuf_data.trial{iT}(lg_idx,perm_shift1+1:end) shuf_data.trial{iT}(lg_idx,1:perm_shift1)];
        %shuf_data.trial{iT}(hg_idx,:) = [shuf_data.trial{iT}(hg_idx,perm_shift2+1:end) shuf_data.trial{iT}(hg_idx,1:perm_shift2)];
        
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
    
    shuf_xf_lg(iShuf,:) = averageXbyYbin(shuf_lg_pwr,obs_xf_phase,pbin_edges);
    shuf_xf_hg(iShuf,:) = averageXbyYbin(shuf_hg_pwr,obs_xf_phase,pbin_edges);
    
    shuf_pac_lg(iShuf) = abs(mean(shuf_lg_pwr.*exp(1i*obs_xf_phase)));
    shuf_pac_hg(iShuf) = abs(mean(shuf_hg_pwr.*exp(1i*obs_xf_phase)));
    
    % get xcorr for this shuffle
    temp_shuf = [];
    if doXCorr
        for iT = length(shuf_data.trial):-1:1
            [temp_shuf(iT,:),~] = xcov(tiedrank(shuf_data.trial{iT}(hg_idx,:)),tiedrank(shuf_data.trial{iT}(lg_idx,:)),maxlag,'coeff');
        end
        shuf_xc(iShuf,:) = nanmean(temp_shuf);
        
    end

    
end

% compute some PAC stats
pacz_lg = (obs_pac_lg-mean(shuf_pac_lg))/std(shuf_pac_lg);
pacz_hg = (obs_pac_hg-mean(shuf_pac_hg))/std(shuf_pac_hg);

pacm_lg = 1-(sum(shuf_pac_lg < obs_pac_lg)/nShuf);
pacm_hg = 1-(sum(shuf_pac_hg < obs_pac_hg)/nShuf);

%% plot
figure;

if doXCorr
    lags = lags.*1/data.fsample;
    
    subplot(221);
    plot(lags,obs_xc,'LineWidth',2)
    hold on
    plot(lags,nanmean(shuf_xc),'k');
    plot(lags,nanmean(shuf_xc)+nanstd(shuf_xc),'k:');
    plot(lags,nanmean(shuf_xc)-nanstd(shuf_xc),'k:');
    
    plot([0 0],get(gca,'YLim'),'--','Color',[0.7 0.7 0.7]);
    

    %set(gca,'FontSize',18,'YLim',[-0.05 0.15]);
    title('xcorr low-high gamma');
    xlabel('time (s)'); ylabel('power rank correlation');
end

subplot(222);
plot(pbin_centers,obs_xf_lg,'LineWidth',2)

hold on
plot(pbin_centers,nanmean(shuf_xf_lg),'k');
plot(pbin_centers,nanmean(shuf_xf_lg)+nanstd(shuf_xf_lg),'k:');
plot(pbin_centers,nanmean(shuf_xf_lg)-nanstd(shuf_xf_lg),'k:');

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('lg pacz %.2f (p = %.3f)',pacz_lg,1-normcdf(pacz_lg,0,1)));


subplot(224);
plot(pbin_centers,obs_xf_hg,'LineWidth',2)

hold on
plot(pbin_centers,nanmean(shuf_xf_hg),'k');
plot(pbin_centers,nanmean(shuf_xf_hg)+nanstd(shuf_xf_hg),'k:');
plot(pbin_centers,nanmean(shuf_xf_hg)-nanstd(shuf_xf_hg),'k:');

set(gca,'XLim',[-pi pi],'XTick',[-pi:pi/2:pi],'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'},'FontSize',18,'YTick',[]);
xlabel('phase (rad)'); ylabel('power (a.u.)');
title(sprintf('hg pacz %.2f (p = %.3f)',pacz_hg,1-normcdf(pacz_hg,0,1)));

subplot(223);
title(sprintf('ALL %d-%d Hz, %d sess %d trls',f_x(1),f_x(2),length(fd),length(data_all.trial)));
axis off;