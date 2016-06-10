%%
cfg = [];
cfg.fsample = 2000;
cfg.trllen = 0.4;
cfg.numtrl = 100;

cfg.s1.freq = 10;
cfg.s1.ampl = 1;
cfg.s1.phase = 0;

cfg.noise.ampl = 0.1;

data = ft_freqsimulation(cfg);
data.label{1} = 's1'; data.label{2} = 's2';

cfg_select = [];
cfg_select.channel = {'s1','s2'};
data = ft_selectdata(cfg_select,data);

%%
clear ALL_icoh ALL_icpsd;

phis = -pi:pi/36:pi;
%phis = -180:0:180;
for iPhi = 1:length(phis)
    
    f1 = 10; phi1 = phis(iPhi);
    f2 = 10; phi2 = 0;
    
    for iT = 1:cfg.numtrl
        
        data.trial{iT}(1,:) = sin(2*pi*f1*data.time{1} + phi1) + cfg.noise.ampl*randn(size(data.time{1}));
        data.trial{iT}(2,:) = sin(2*pi*f2*data.time{1} + phi2) + cfg.noise.ampl*randn(size(data.time{1}));
        
    end
    
    %%
    
    cfg_freqanalysis              = [];
    cfg_freqanalysis.output       = 'fourier';
    cfg_freqanalysis.method       = 'mtmfft';
    cfg_freqanalysis.taper        = 'dpss';
    cfg_freqanalysis.pad          = 10;
    cfg_freqanalysis.tapsmofrq    = 8; % def 8
    cfg_freqanalysis.foi          = 1:1:100;
    cfg_freqanalysis.keeptrials   = 'yes';
    cfg_freqanalysis.channel      = {'s1','s2'};
    cfg_freqanalysis.channelcmb   = 'all';
    
    freq = ft_freqanalysis(cfg_freqanalysis, data);
    
    
    %% compute imaginary coherence: NOTE should add variability estimate
    cfg_conn       = []; cfg_conn.method    = 'coh'; cfg_conn.complex = 'angle';
    cohi           = ft_connectivityanalysis(cfg_conn, freq);
    
    ALL_icoh(iPhi) = cohi.cohspctrm(1,2,nearest_idx3(f2,cohi.freq));
    %plot(cohi.freq,sq(cohi.cohspctrm(1,2,:)));
    
    %% alternate
    clear these_icpsd
    for iT = cfg.numtrl:-1:1
        [this_icpsd,F] = cpsd(data.trial{iT}(1,:),data.trial{iT}(2,:),length(data.trial{iT}(1,:)),[],[],cfg.fsample);
        this_icpsd = angle(this_icpsd);
        these_icpsd(iT) = this_icpsd(nearest_idx3(f2,F));
    end
   ALL_icpsd(iPhi) = circ_mean(these_icpsd');
end

%% 
subplot(221)
plot(phis,ALL_icoh)
title('ft angle')

subplot(222)
plot(phis,ALL_icpsd)
title('icpsd')