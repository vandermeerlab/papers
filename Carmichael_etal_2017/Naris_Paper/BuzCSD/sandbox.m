%%
% restoredefaultpath;
% addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
% addpath('D:\My_Documents\GitHub\fieldtrip');
% ft_defaults;

% rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal');

addpath('/Users/ecarmichael/Dropbox/Matlab/BuzCSD');

%%
cfg = []; cfg.method = 'superimposed';
cfg.fsample = 2000; cfg.trllen = 0.4;
cfg.s1.freq = 50; cfg.s1.ampl = 100; % units will be uV
cfg.s2.freq = 50; cfg.s2.ampl = 100; cfg.s2.phase = degtorad(10);
cfg.s3.freq = 50; cfg.s3.ampl = 100; cfg.s3.phase = degtorad(1);

d_orig = ft_freqsimulation(cfg);
%% create wavelet and replicate into 64 channels
d = d_orig;

ch_label_to_replicate = 's1';
nReplicates = 64;
cfg_wv.sd = 6; % SD of gaussian
cfg_wv.pad = 200; % in samples -- hack (should be in time)
cfg_wv.na = 1; % noise amplitude

ch_idx = strmatch(ch_label_to_replicate,d.label);

for iR = 1:nReplicates
    d.trial{1}(iR,:) = d_orig.trial{1}(ch_idx,:);
    
    d.trial{1}(iR,:) = d.trial{1}(iR,:).*gausswin(length(d.trial{1}(iR,:)),cfg_wv.sd)';
    
    d.trial{1}(iR,:) = d.trial{1}(iR,:) *(1+iR/64);
    
    d.trial{1}(iR,1:cfg_wv.pad) = 0; % hack!
    d.trial{1}(iR,end-cfg_wv.pad+1:end) = 0;
    
    d.trial{1}(iR,:) = d.trial{1}(iR,:) + cfg_wv.na*randn(size(d.trial{1}(iR,:))); % add noise
end
% invert a couple to create source
%invrt = [7 8]; % note these will fall off the final output (edge
% electrodes) but nonzero CSD still shows up
invrt = [4];
d_og = d.trial{1}(invrt,:) ;
d.trial{1}(invrt,:) = -d.trial{1}(invrt,:);


figure(1);
subplot(421);
plot(d.time{1},d.trial{1}(1:8,:));

% create sites X time matrix (input to csd5pt.m)
sites = 1:8;
xcpots = d.trial{1}(sites,:);
xcpots(5,:) = NaN;
subplot(421)
plot(d.time{1}, xcpots)
% CSD
elec_sep_mm = sqrt(0.2^2+0.2^2);
[csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);

subplot(422)
imagesc(d.time{1},CSDelecinds,csd);
colorbar; % note units now in uV/mm^2... but numbers here seem very small?


% phase offset %%%%%%%%%%%%

invrt = 4;
d.trial{1}(invrt,:) =  d_orig.trial{1}(3,:);

d.trial{1}(invrt,:) = d.trial{1}(invrt,:).*gausswin(length(d.trial{1}(invrt,:)),cfg_wv.sd)';

d.trial{1}(invrt,1:cfg_wv.pad) = 0; % hack!
d.trial{1}(invrt,end-cfg_wv.pad+1:end) = 0;

d.trial{1}(invrt,:) = d.trial{1}(invrt,:) *(1+invrt/64);

d.trial{1}(invrt,:) = d.trial{1}(invrt,:) + cfg_wv.na*randn(size(d.trial{1}(invrt,:))); % add noise

disp(['Channel ' num2str(invrt) ' Phase off by ' num2str(cfg.s2.phase)])

subplot(423)
xcpots_phase = d.trial{1}(sites,:);
plot(d.time{1}, xcpots_phase)

sites = 1:8;
xcpots = d.trial{1}(sites,:);

subplot(423)
plot(d.time{1}, xcpots)
title(['Channel ' num2str(invrt) ' Phase off by ' num2str(rad2deg(cfg.s2.phase))]);

% CSD
elec_sep_mm = sqrt(0.2^2+0.2^2);
[csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);

subplot(424)
imagesc(d.time{1},CSDelecinds,csd);
colorbar; % note units now in uV/mm^2... but numbers here seem very small?


% phae off by 10
d.trial{1}(invrt,:) =  d_orig.trial{1}(4,:);

d.trial{1}(invrt,:) = d.trial{1}(invrt,:).*gausswin(length(d.trial{1}(invrt,:)),cfg_wv.sd)';

d.trial{1}(invrt,1:cfg_wv.pad) = 0; % hack!
d.trial{1}(invrt,end-cfg_wv.pad+1:end) = 0;

d.trial{1}(invrt,:) = d.trial{1}(invrt,:) + cfg_wv.na*randn(size(d.trial{1}(invrt,:))); % add noise

disp(['Channel ' num2str(invrt) ' Phase off by ' num2str(cfg.s3.phase)])

subplot(425)
xcpots_phase = d.trial{1}(sites,:);
plot(d.time{1}, xcpots_phase)
title(['Channel ' num2str(invrt) ' Phase off by ' num2str(rad2deg(cfg.s3.phase))]);

sites = 1:8;
xcpots = d.trial{1}(sites,:);

% subplot(426)
% plot(d.time{1}, xcpots)
% CSD
elec_sep_mm = sqrt(0.2^2+0.2^2);
[csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);

subplot(426)
imagesc(d.time{1},CSDelecinds,csd);
colorbar; % note units now in uV/mm^2... but numbers here seem very small?


%% check all ranges
val_out_csd_min = cell(1,8);
val_out_csd_max = cell(1,8);
val_out_rad = cell(1,8);
val_out_csd_min_abs = cell(1,8);

for iAmp = 1:8
    for irad = 0:5:360;
        cfg = []; cfg.method = 'superimposed';
        cfg.fsample = 2000; cfg.trllen = .4;
        cfg.s1.freq = 50; cfg.s1.ampl = 100; cfg.s1.phase = degtorad(irad);
        cfg.s2.freq = 50; cfg.s2.ampl = 100; cfg.s2.phase = degtorad(90);
        cfg.s3.freq = 50    ; cfg.s3.ampl = 100; cfg.s3.phase = degtorad(5);
        d_orig = ft_freqsimulation(cfg);
        invrt = 4;
        %% get the channels and convolve
        d.trial{1}(invrt,:) =  d_orig.trial{1}(2,:);
        
        d.trial{1}(invrt,:) = d.trial{1}(invrt,:).*gausswin(length(d.trial{1}(invrt,:)),cfg_wv.sd)';
        
        d.trial{1}(invrt,1:cfg_wv.pad) = 0; % hack!
        d.trial{1}(invrt,end-cfg_wv.pad+1:end) = 0;
        
        d.trial{1}(invrt,:) = d.trial{1}(invrt,:) + cfg_wv.na*randn(size(d.trial{1}(invrt,:))); % add noise
        
        %     disp(['Channel ' num2str(invrt) ' Phase off by ' num2str(cfg.s1.phase)])
        d.trial{1}(invrt,:) = d.trial{1}(invrt,:) *(1+iAmp);
        
        xcpots_phase = d.trial{1}(sites,:);
        
        sites = 1:8;
        xcpots = d.trial{1}(sites,:);
        
        % CSD
        elec_sep_mm = sqrt(0.2^2+0.2^2);
        [csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);
        val_out_csd_min{iAmp} = [val_out_csd_min{iAmp}, min(min(csd))];
        val_out_csd_max{iAmp} = [val_out_csd_max{iAmp}, max(max(csd))];
        val_out_csd_min_abs{iAmp} = [val_out_csd_min_abs{iAmp}, min(min(abs(csd)))];
        val_out_rad{iAmp} = [val_out_rad{iAmp}, irad];
    end
end

figure
c_ord = linspecer(10);
sub_place = [1 2 4 5 7 8 10 11];
for iAmp = 1:8
    subplot(4,3,sub_place(iAmp))
    plot(val_out_rad{iAmp}, val_out_csd_min{iAmp}(1,:), val_out_rad{iAmp}, val_out_csd_max{iAmp}, '-*', 'color', c_ord(iAmp,:), 'linewidth', 2)
    ylabel(['amp * ' num2str(1+iAmp)]);
end

% get the difference between the max and thin 
subplot(4,3,[3 6])
hold on
for iAmp = 1:8
    plot(val_out_rad{iAmp}, (val_out_csd_max{iAmp}-val_out_csd_min{iAmp}),'*', 'color', c_ord(iAmp,:))
end

subplot(4,3,[9 12])
hold on
for iAmp = 1:8
    plot(val_out_rad{iAmp}, (val_out_csd_max{iAmp}-val_out_csd_min_abs{iAmp}),'*', 'color', c_ord(iAmp,:))
end


%% match the data best fit
% create 8 channels with offsets similar to the cpsd phase differences and
% the power differences.  
cfg_wv.sd = 6; % SD of gaussian
cfg_wv.pad = 200; % in samples -- hack (should be in time)
cfg_wv.na = 1; % noise amplitude
amp = 100; 
cfg = []; cfg.method = 'superimposed';
cfg.fsample = 2000; cfg.trllen = 0.4;
cfg.s1.freq = 50; cfg.s1.ampl = amp*1;  cfg.s1.phase = degtorad(10)*-1;
cfg.s2.freq = 50; cfg.s2.ampl = amp*1.1; cfg.s2.phase = degtorad(9)*-1;
cfg.s3.freq = 50; cfg.s3.ampl = amp*1.15; cfg.s3.phase = degtorad(8)*-1;
d_1_3 = ft_freqsimulation(cfg);
 % set 4-6
 cfg = []; cfg.method = 'superimposed';
cfg.fsample = 2000; cfg.trllen = 0.4;
cfg.s1.freq = 50; cfg.s1.ampl = amp*1.18;  cfg.s1.phase = degtorad(6)*-1;
cfg.s2.freq = 50; cfg.s2.ampl = amp*2; cfg.s2.phase = degtorad(4)*-1;
cfg.s3.freq = 50; cfg.s3.ampl = amp*3; cfg.s3.phase = degtorad(3.5)*-1;
d_4_6 = ft_freqsimulation(cfg);
 % set 4-6
 cfg = []; cfg.method = 'superimposed';
cfg.fsample = 2000; cfg.trllen = 0.4;
cfg.s1.freq = 50; cfg.s1.ampl = amp*4.2;  cfg.s1.phase = degtorad(2)*-1;
cfg.s2.freq = 50; cfg.s2.ampl = amp*6; cfg.s2.phase = degtorad(0)*-1;
d_7_8 = ft_freqsimulation(cfg);

% restructure
d = d_1_3;

d_orig.trial{1}(1,:) = d_1_3.trial{1}(1,:);
d_orig.trial{1}(2,:) = d_1_3.trial{1}(2,:);
d_orig.trial{1}(3,:) = d_1_3.trial{1}(3,:);
d_orig.trial{1}(4,:) = d_4_6.trial{1}(1,:);
d_orig.trial{1}(5,:) = d_4_6.trial{1}(2,:);
d_orig.trial{1}(6,:) = d_4_6.trial{1}(3,:);
d_orig.trial{1}(7,:) = d_7_8.trial{1}(1,:);
d_orig.trial{1}(8,:) = d_7_8.trial{1}(2,:);

%%
for iR = 1:8
    d.trial{1}(iR,:) = d_orig.trial{1}(iR,:);
    
    d.trial{1}(iR,:) = d.trial{1}(iR,:).*gausswin(length(d.trial{1}(iR,:)),cfg_wv.sd)';
    
%     d.trial{1}(iR,:) = d.trial{1}(iR,:) *(1+iR/64);
    
    d.trial{1}(iR,1:cfg_wv.pad) = 0; % hack!
    d.trial{1}(iR,end-cfg_wv.pad+1:end) = 0;
    
    d.trial{1}(iR,:) = d.trial{1}(iR,:) + cfg_wv.na*randn(size(d.trial{1}(iR,:))); % add noise
end

figure(2);
subplot(421);
plot(d.time{1},d.trial{1}(1:8,:));

% create sites X time matrix (input to csd5pt.m)
sites = 1:8;
xcpots = d.trial{1}(sites,:);

subplot(421)
plot(d.time{1}, xcpots)
% CSD
elec_sep_mm = sqrt(0.2^2+0.2^2);
[csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm);
label = num2cell(1:8);
legend('show')

subplot(422)
imagesc(d.time{1},CSDelecinds,csd);
colorbar; % note units now in uV/mm^2... but numbers here seem very small?
