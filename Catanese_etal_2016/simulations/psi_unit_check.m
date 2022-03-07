%%
cfg = [];
cfg.fsample = 2000;
cfg.trllen = 1;
cfg.numtrl = 100;

cfg.s1.freq = 10;
cfg.s1.ampl = 1;
cfg.s1.phase = 0;

cfg.noise.ampl = 0.01;

data = ft_freqsimulation(cfg);
data.label{1} = 's1'; data.label{2} = 's2';

cfg_select = [];
cfg_select.channel = {'s1','s2'};
data = ft_selectdata(cfg_select,data);

%%
clear ALL_icoh ALL_icpsd;

t_shift = 0.004; % in seconds
t_shift_samples = round(t_shift/(1/cfg.fsample));

fvec = 10:10:100;

for iF = 1:length(fvec)
    
    this_f = fvec(iF);
    
    f1 = this_f-5;
    f2 = this_f+5;
    
    tvec = data.time{1};
    chrp = chirp(tvec,f1,tvec(end),f2);
    
    for iT = 1:cfg.numtrl
        
        data.trial{iT}(1,:) = chrp + cfg.noise.ampl*randn(size(data.time{1}));
        data.trial{iT}(2,:) = circshift(chrp,[1 t_shift_samples]) + cfg.noise.ampl*randn(size(data.time{1}));
        
    end

    %%
    
    cfg_freqanalysis              = [];
    cfg_freqanalysis.output       = 'fourier';
    cfg_freqanalysis.method       = 'mtmfft';
    cfg_freqanalysis.taper        = 'hanning'; % hanning, dpss
    cfg_freqanalysis.pad          = 10;
    cfg_freqanalysis.tapsmofrq    = 4; % def 8, for dpss only
    cfg_freqanalysis.foi          = 1:1:100;
    cfg_freqanalysis.keeptrials   = 'yes';
    cfg_freqanalysis.channel      = {'s1','s2'};
    cfg_freqanalysis.channelcmb   = 'all';
    
    freq = ft_freqanalysis(cfg_freqanalysis, data);
    
    
    %% compute imaginary coherence: NOTE should add variability estimate
    cfg_psi      = []; cfg_psi.method    = 'psi'; cfg_psi.bandwidth = 4;
    psi           = ft_connectivityanalysis(cfg_psi, freq);
    
    cfg_coh     = []; cfg_coh.method    = 'coh'; cfg_coh.complex = 'angle';
    coh           = ft_connectivityanalysis(cfg_coh, freq);
   
    phase_diff = sq(coh.cohspctrm(1,2,:)).*(360/(2*pi));
    phase_slope = sq(psi.psispctrm(1,2,:)).*(360/(2*pi));
    
    f_idx = nearest_idx3(this_f,coh.freq);
    ALL_phase_diff(iF) = phase_diff(f_idx);
    ALL_phase_slope(iF) = phase_slope(f_idx);
end

%%
subplot(221);
plot(fvec,ALL_phase_diff);
title('phase diffs');
hold on;
plot(fvec,(fvec.*t_shift)*360,'r');

subplot(222);
plot(fvec,ALL_phase_slope/(2*pi));
title('phase slopes');


%% 
subplot(221)
plot(psi.freq,sq(psi.psispctrm(1,2,:)).*(360/(2*pi)))
title('psi raw')

subplot(222)
plot(coh.freq,sq(coh.cohspctrm(1,2,:)).*(360/(2*pi)));
title('icoh raw')

