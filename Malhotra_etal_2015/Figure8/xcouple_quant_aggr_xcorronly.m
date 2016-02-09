clear all;
rng('shuffle');
%addpath('C:\Users\mvdm\Dropbox\projects\Sushant\shared');

fd = {'C:\data\R014-2012-09-14_promoted', ...
    'C:\data\R014-2012-09-16_promoted', ...
    'C:\data\R016-2012-10-05_promoted', ...
    'C:\data\R016-2012-10-08_promoted', ...
    'C:\data\R018-2012-12-08_promoted', ...
    'C:\data\R018-2012-12-14_promoted'};

%cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
%load fd_isidro
%cd('D:\My_Documents\My Dropbox\projects\Sushant\2014-09-19') % athena
cd('D:\My_Documents\Dropbox\projects\Sushant\2014-09-19'); % equinox
load fd; % list of files to process, obtained from remove_cueviolation.m

% only certain rat
%fd = fd(end-2:end); % R014
%fd = fd(1); % R020
%fd = fd(2:8); % R018
%fd = fd(9:13); % R016

%%
twin = [0 1.25];
%twin = [-1 3];
f_lg = [50 65];
f_hg = [70 85];

f_x = [3 5];

doXCorr = 1;

% first, build dataset
data_all = []; event_all = [];
for iFD = 1:length(fd)
    
   cd(fd{iFD});

    run(FindFile('*keys.m'));
    fc = ExpKeys.goodGamma(1);
    data = ft_read_neuralynx_interp(fc);
 
    cur_label = data.label{1};
    
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
    d = fdesign.bandpass('N,F3dB1,F3dB2',4,f_x(1),f_x(2),2e3);
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

% first get observed PAC value
trial_mat = cell2mat(temp.trial); % trials concatenated

lg_idx = strmatch('lg_pow',temp.label);
obs_lg_pwr = trial_mat(lg_idx,:); % note, ft_selectdata switches labels!!

hg_idx = strmatch('hg_pow',temp.label);
obs_hg_pwr = trial_mat(hg_idx,:); % note, ft_selectdata switches labels!!

% improved xcorr that does not go across trials
if doXCorr
    for iT = length(temp.trial):-1:1
        [obs_xc(iT,:),lags] = xcov(tiedrank(temp.trial{iT}(hg_idx,:)),tiedrank(temp.trial{iT}(lg_idx,:)),maxlag,'coeff');
    end
    obs_xc = nanmean(obs_xc);
    shuf_xc = zeros(nShuf,length(lags));
else
    obs_xc = [];
end

% shuffle
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
    
    % get xcorr for this shuffle
    temp_shuf = [];
    if doXCorr
        for iT = length(shuf_data.trial):-1:1
            [temp_shuf(iT,:),~] = xcov(tiedrank(shuf_data.trial{iT}(hg_idx,:)),tiedrank(shuf_data.trial{iT}(lg_idx,:)),maxlag,'coeff');
        end
        shuf_xc(iShuf,:) = nanmean(temp_shuf);
        
    end
    
end


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
    

    set(gca,'FontSize',18,'YLim',[-0.05 0.15]);
    title('xcorr low-high gamma');
    xlabel('time (s)'); ylabel('power rank correlation');
end
