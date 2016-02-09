clear all;
rng('shuffle');
%addpath('C:\Users\mvdm\Dropbox\projects\Sushant\shared');

cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-17');
load fd_isidro
%cd('D:\My_Documents\My Dropbox\projects\Sushant\2014-09-19') % athena
%load fd; % list of files to process, obtained from remove_cueviolation.m

% load all promoted sessions
%cd('C:\Users\mvdm\Dropbox\projects\Sushant\2014-10-20');
%load all_files

% only certain rat
%files = fd(end-2:end); % R014
%files = fd(1); % R020
%files = fd(2:8); % R018
%files = fd(9:13); % R016
%% for files that passed behavior
FDtoRat(1) = 20;
FDtoRat(2:8) = 18;
FDtoRat(9:13) = 16;
FDtoRat(14:16) = 14;

%% for all
%FDtoRat(1:8) = 16;
%FDtoRat(9:11) = 14;
%FDtoRat(12:18) = 18;
%FDtoRat(19:25) = 20;

%% build data set -- all trials together!
data_np_all = []; data_cue_all = []; data_rest_all = [];
data_np_allID = []; data_cue_allID = []; data_rest_allID = [];
for iFD = 1:length(fd)
    
    cd(fd{iFD});

    run(FindFile('*keys.m'));
    fc = ExpKeys.goodGamma(1);
    data = ft_read_neuralynx_interp(fc);
 
    data.label{1} = 'vStr';

    % nosepokes   
    cfg = [];
    cfg.trialdef.cue = {'c1','c3','c5','lo','hi'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
    cfg.trialdef.block = 'both';
    cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
    cfg.trialdef.location = 'both';
    cfg.trialfun = 'ft_trialfun_lineartracktone2';
    cfg.trialdef.hdr = data.hdr;
    cfg.trialdef.pre = 0; cfg.trialdef.post = 1.25;
    
    [trl, event] = ft_trialfun_lineartracktone2(cfg); 
    
    cfg = []; cfg.trl = trl;
    this_data = ft_redefinetrial(cfg,data);
    
    % remove NaN trials?
    mf = @(x) ~any(isnan(x(:)));
    keep_idx = cellfun(mf,this_data.trial);
    fprintf('\n*This session has %d out of %d trials without NaNs...\n\n',sum(keep_idx),length(keep_idx));
    
    cfg = []; cfg.trials = find(keep_idx);
    this_data = ft_selectdata(cfg,this_data);
    
    if isempty(data_np_all)
        data_np_all = this_data;
        data_np_allID = repmat(FDtoRat(iFD),[length(this_data.trial) 1]);
    else
        data_np_all = ft_appenddata([],data_np_all,this_data);
        data_np_allID = cat(1,data_np_allID,repmat(FDtoRat(iFD),[length(this_data.trial) 1]));
    end
    
    % cue
    cfg.trialdef.cue = {'c1','c3','c5','lo','hi'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
    cfg.trialdef.block = 'both';
    cfg.trialdef.location = 'both';
    cfg.trialfun = 'ft_trialfun_lineartracktone2';
    cfg.trialdef.hdr = data.hdr;
    cfg.trialdef.eventtype = 'cue';
    cfg.trialdef.pre = 0.5; cfg.trialdef.post = 0.75;
    
    [trl, event] = ft_trialfun_lineartracktone2(cfg); 
    
    cfg = []; cfg.trl = trl;
    this_data = ft_redefinetrial(cfg,data);
    
    % remove NaN trials?
    mf = @(x) ~any(isnan(x(:)));
    keep_idx = cellfun(mf,this_data.trial);
    fprintf('\n*This session has %d out of %d trials without NaNs...\n\n',sum(keep_idx),length(keep_idx));
    
    cfg = []; cfg.trials = find(keep_idx);
    this_data = ft_selectdata(cfg,this_data);
    
    % build data set
    if isempty(data_cue_all)
        data_cue_all = this_data;
        data_cue_allID = repmat(FDtoRat(iFD),[length(this_data.trial) 1]);
    else
        data_cue_all = ft_appenddata([],data_cue_all,this_data);
        data_cue_allID = cat(1,data_cue_allID,repmat(FDtoRat(iFD),[length(this_data.trial) 1]));
    end
    
    % rest
    cfg.twin = [-1 1.5]; % original: [-0.5 0.75]
    
    if ~isempty(ExpKeys.Prerecord)
        pre_t = ExpKeys.Prerecord(1)+5-cfg.twin(1):diff(cfg.twin)*2:ExpKeys.Prerecord(2)-cfg.twin(2);
    else
        pre_t = [];
    end
    
    post_t = ExpKeys.Postrecord(1)+5-cfg.twin(1):diff(cfg.twin)*2:ExpKeys.Postrecord(2)-cfg.twin(2);
    
    cfg = [];
    cfg.t = cat(2,pre_t,post_t);
    cfg.twin = [-1 1.5];
    cfg.mode = 'nlx';
    cfg.hdr = data.hdr;
    
    trl = ft_maketrl(cfg);
    
    cfg = []; cfg.trl = trl;
    this_data = ft_redefinetrial(cfg,data);
    
    % remove NaN trials?
    mf = @(x) ~any(isnan(x(:)));
    keep_idx = cellfun(mf,this_data.trial);
    fprintf('\n*This session has %d out of %d trials without NaNs...\n\n',sum(keep_idx),length(keep_idx));
    
    cfg = []; cfg.trials = find(keep_idx);
    this_data = ft_selectdata(cfg,this_data);
    
    if isempty(data_rest_all)
        data_rest_all = this_data;
        data_rest_allID = repmat(FDtoRat(iFD),[length(this_data.trial) 1]);
    else
        data_rest_all = ft_appenddata([],data_rest_all,this_data);
        data_rest_allID = cat(1,data_rest_allID,repmat(FDtoRat(iFD),[length(this_data.trial) 1]));
    end
    
end

%% plot PSDs
rats = [14 16 18 20];
what = {'data_rest_all','data_cue_all','data_np_all'};
cols = 'rgbc';
for iRat = 1:length(rats)
    
   cfg = [];
   cfg.method = 'mtmfft';
   cfg.output = 'pow';
   cfg.keeptrials = 'no';
   cfg.taper = 'hanning';
   cfg.tapsmofrq = 2;
   cfg.foi = 0:0.5:150;
   
   for iW = 1:length(what)
   
       trial_idx = eval(['find(' what{iW} 'ID == ' num2str(rats(iRat)) ');']); 
       data_temp = eval(what{iW});
       
       cfg.trials = trial_idx;
    
       F = ft_freqanalysis(cfg,data_temp);
    
       figure(iW); hold on; grid on;
       h(iW,iRat) = plot(F.freq,10*log10(F.powspctrm),cols(iRat),'LineWidth',2);
    
       set(gca,'FontSize',18,'LineWidth',1,'XTick',0:25:150,'YTick',0:10:40,'YLim',[-5 40],'XLim',[0 150]);
       title(what{iW});
       xlabel('Frequency (Hz)'); ylabel('power (10*log10(uV))');
       
   end
    
end
legend(h(3,:),cellfun(@num2str,num2cell(rats),'UniformOutput',0));
