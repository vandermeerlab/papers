function out = hdfit_crossval_func(cfg_in)

cfg_def = [];
cfg_def.datapath = 'D:\My_Documents\Dropbox\projects\HDfit\data';
cfg_def.session = 3; % which session to load
cfg_def.mode = 'kalmanwrapped'; % 'kalmanwrapped', 'smooth', 'kalman'
cfg_def.debug = 0;
cfg_def.target_session = 'std'; % 'std', 'laser'
cfg_def.gainbin_centers = 0.9:0.025:1.1;
cfg_def.driftbin_centers = -1.5:0.25:1.5;
cfg_def.simparams = []; % if specified, use simulated data
cfg_def.models = 1:4; % models to run

cfg_master = ProcessConfig(cfg_def,cfg_in);


global param; % to access internals of crossval()
global param_count;
global param_hist;

%% cfg
cfg_master.subsample_factor = 7;
cfg_master.smoothwin = 11; % samples of 1/60


%% data loading
cfg_load = [];
cfg_load.fd = cfg_master.datapath;
cfg_load.ATIshift = 2;

data = WillDataLoader(cfg_load,cfg_master.session);


%% preprocess HD
cfg_pp = [];
cfg_pp.smoothwin = cfg_master.smoothwin;
cfg_pp.subsample_factor = cfg_master.subsample_factor;
cfg_pp.debug = cfg_master.debug;
cfg_pp.mode = cfg_master.mode; % 'smooth', 'kalman'

data = WillPreProcData(cfg_pp,data);

%% sanity check: should be able to get HD back from AHV
fnames = fieldnames(data);
for iF = 1:length(fnames)
    
    this_data = data.(fnames{iF});
    
    cfg_test = [];
    cfg_test.hd0 = this_data.hd_ss.data(1);
    hd_ss_test = AHVtoHDfast(this_data.ahv_ss,[1 1 0 cfg_test.hd0 0]);
    
    if cfg_master.debug
        figure;
        plot(this_data.hd_ss,'.k'); hold on;
        plot(this_data.ahv_ss.tvec,hd_ss_test,'.r');
        title(fnames{iF});
    end
    
end

%% estimate firing rates (spike density functions) from binned source data
cfg_sdf = [];
cfg_sdf.subsample_factor = cfg_master.subsample_factor;
data = GetSDFfromHDdata(cfg_sdf,data);


%% compute TCs (real)
cfg_tc = [];
cfg_tc.bin_edges = 0:3.6:360;
cfg_tc.bin_centers = 1/2+cfg_tc.bin_edges(1:end-1);
cfg_tc.occ_dt = 1/60;
cfg_tc.nseg = 8;

for iF = 1:length(fnames)
    
    tc = TuningCurvesSDF(cfg_tc,data.(fnames{iF}).sdf,data.(fnames{iF}).hd);
    tc.xbin = cfg_tc.bin_centers;
    
    data.(fnames{iF}).tc = tc;
    
    % plot overall and short time segment TCs
    if cfg_master.debug
        figure;
        subplot(221);
        plot(tc.xbin,tc.tc'); xlim([0 360]);
        title(fnames{iF});
        cmap = colormap(jet);
    end
    
    nCells = size(data.(fnames{iF}).obs_fr,1);

    for iC = 1:nCells
        
        if cfg_master.debug
            subplot(2,2,1+iC);
        end
        
        seg = linspace(data.(fnames{iF}).sdf.tvec(1),data.(fnames{iF}).sdf.tvec(end),cfg_tc.nseg+1);
        
        for iSeg = 1:cfg_tc.nseg
            this_sdf = restrict(data.(fnames{iF}).sdf,seg(iSeg),seg(iSeg+1));
            this_hd = restrict(data.(fnames{iF}).hd,seg(iSeg),seg(iSeg+1));
            this_tc = TuningCurvesSDF(cfg_tc,this_sdf,this_hd);
            
            if cfg_master.debug
                plot(tc.xbin,this_tc.tc(iC,:)','Color',cmap(round(iSeg*(64/cfg_tc.nseg)),:)); xlim([0 360]);
                hold on;
            end
            
            data.(fnames{iF}).tc_seg{iSeg} = this_tc.tc;
            
            % if laser session, compute shifts relative to std session
            if strcmp(fnames{iF},'laser')
                [tc_corr,tc_shift] = TCcorr(data.std.tc.tc,data.(fnames{iF}).tc_seg{iSeg});
                data.(fnames{iF}).tc_segcorr(iSeg) = tc_corr;
                data.(fnames{iF}).tc_segshift(iSeg) = tc_shift;
            end
            
        end    
        
    end
end

%% generate synthetic firing rates based on TCs from specified session
% with known gain, drift parameters
% used for parameter recovery analysis
% NOTE assumes we know internal HD at start of session!
if ~isempty(cfg_master.simparams)
    
    cfg_sim = [];
    cfg_sim.params = cfg_master.simparams; % gain_left, gain_right, drift
    cfg_sim.sourceTC = 'std'; % which TCs to use
    for iF = 1:length(fnames)
        
        cfg_sim.hd0 = data.(fnames{iF}).hd.data(1);
        cfg_sim.hd0_idx = 0; % don't shift tuning curves, take as is
        
        sim_internal_hd = AHVtoHDfast(data.(fnames{iF}).ahv,cat(2,cfg_sim.params,cfg_sim.hd0,cfg_sim.hd0_idx));
        sim_internal_hd = tsd(data.(fnames{iF}).ahv.tvec,sim_internal_hd);
        
        data.(fnames{iF}).sim_sdf = GenerateSDFfromTC([],sim_internal_hd,data.(cfg_sim.sourceTC).tc);
        
        data.(fnames{iF}).sim_sdf_ss = data.(fnames{iF}).sim_sdf;
        data.(fnames{iF}).sim_sdf_ss.tvec = data.(fnames{iF}).sim_sdf_ss.tvec(1:cfg_master.subsample_factor:end);
        data.(fnames{iF}).sim_sdf_ss.data = data.(fnames{iF}).sim_sdf_ss.data(:,1:cfg_master.subsample_factor:end);
        
    end
end
%% tuning curves from simulated SDFs
if ~isempty(cfg_master.simparams)
    for iF = 1:length(fnames)
        
        tc = TuningCurvesSDF(cfg_tc,data.(fnames{iF}).sim_sdf,data.(fnames{iF}).hd);
        tc.xbin = cfg_tc.bin_centers;
        
        data.(fnames{iF}).tc_sim = tc;
        
        % plot overall and short time segment TCs
        if cfg_master.debug
            figure;
            subplot(221);
            plot(tc.xbin,tc.tc'); xlim([0 360]);
            title(fnames{iF});
            cmap = colormap(jet);
        end
        
        nCells = size(data.(fnames{iF}).obs_fr,1);
        
        for iC = 1:nCells
            
            if cfg_master.debug
                subplot(2,2,1+iC);
            end
            
            seg = linspace(data.(fnames{iF}).sim_sdf.tvec(1),data.(fnames{iF}).sim_sdf.tvec(end),cfg_tc.nseg+1);
            
            for iSeg = 1:cfg_tc.nseg
                this_sdf = restrict(data.(fnames{iF}).sim_sdf,seg(iSeg),seg(iSeg+1));
                this_hd = restrict(data.(fnames{iF}).hd,seg(iSeg),seg(iSeg+1));
                this_tc = TuningCurvesSDF(cfg_tc,this_sdf,this_hd);
                
                if cfg_master.debug
                    plot(tc.xbin,this_tc.tc(iC,:)','Color',cmap(round(iSeg*(64/cfg_tc.nseg)),:)); xlim([0 360]);
                    hold on;
                end
                
                data.(fnames{iF}).tc_seg_sim{iSeg} = this_tc.tc;
                
                % if laser session, compute shifts relative to std session
                if strcmp(fnames{iF},'laser')
                    [tc_corr,tc_shift] = TCcorr(data.std.tc.tc,data.(fnames{iF}).tc_seg_sim{iSeg});
                    data.(fnames{iF}).tc_segcorr_sim(iSeg) = tc_corr;
                    data.(fnames{iF}).tc_segshift_sim(iSeg) = tc_shift;
                end
            end
            
        end
    end
end

%% direct approach -- compute goodness of fit for various parameter combinations
cfg_param.target_session = cfg_master.target_session; % fit data from this session
cfg_param.reference_session = 'std'; % get TCs from here

if isempty(cfg_master.simparams)
    cfg_param.data_to_fit = 'sdf_ss'; % sdf_ss (actual spiking data) or sim_sdf_ss (simulated)
else
    cfg_param.data_to_fit = 'sim_sdf_ss'; % sdf_ss (actual spiking data) or sim_sdf_ss (simulated)
end

clear model;
model(1).gain_l = 1;
model(1).gain_r = 1;
model(1).drift = 0;
model(1).label = 'M0: fixed gain and drift';

model(2).gain_l = 1;
model(2).gain_r = 1;
model(2).drift = cfg_master.driftbin_centers;
model(2).label = 'M1: fixed gain';

model(3).gain_l = cfg_master.gainbin_centers;
model(3).gain_r = cfg_master.gainbin_centers;
model(3).drift = 0;
model(3).label = 'M2: fixed drift';

model(4).gain_l = cfg_master.gainbin_centers;
model(4).gain_r = cfg_master.gainbin_centers;
model(4).drift = cfg_master.driftbin_centers;
model(4).label = 'M3: full model';

clear err all_param all_wraperr;
for iM = cfg_master.models

    cfg_param.gain_l = model(iM).gain_l;
    cfg_param.gain_r = model(iM).gain_r;
    cfg_param.drift = model(iM).drift;
    
    ref_hd0 = data.(cfg_param.target_session).hd_ss.data(1);
    cfg_param.hd0 = ref_hd0-30:10:ref_hd0+30;
    %cfg_param.hd0 = ref_hd0;
    
    cfg_param.histBins{1} = model(4).gain_l; cfg_param.histBins{2} = model(4).gain_r;
    cfg_param.histBins{3} = model(4).drift; cfg_param.histBins{4} = cfg_param.hd0;
    
    % obtain correctly shifted tuning curves - matching first segment of target
    % session (should check that looks OK)
    switch cfg_param.data_to_fit
        case 'sdf_ss'
            cfg_param.tc.tc = circshift(data.(cfg_param.reference_session).tc.tc,[0 data.(cfg_param.target_session).tc_segshift(1)]);
            cfg_param.tc.xbin = data.(cfg_param.reference_session).tc.xbin;
        case 'sim_sdf_ss'
            cfg_param.tc.tc = circshift(data.(cfg_param.reference_session).tc.tc,[0 data.(cfg_param.target_session).tc_segshift_sim(1)]);
            cfg_param.tc.xbin = data.(cfg_param.reference_session).tc.xbin;
    end
    
    %plot(cfg_param.tc.xbin,cfg_param.tc.tc')
    %figure;
    %plot(cfg_param.tc.xbin,data.(cfg_param.target_session).tc_seg_sim{1}')
    
    % generate param combos (particles) to test
    [p,q,v,r,s] = ndgrid(cfg_param.gain_l,cfg_param.gain_r,cfg_param.drift,cfg_param.hd0,0);
    particle_vals = cat(2,p(:),q(:),v(:),r(:),s(:));
    nP = size(particle_vals,1);
    
    %%% crossvalidation starts here %%%
    
    % have: time points Xtrain with corresponding SDF values Ytrain
    % need: Yfit for excluded points Xtest
    
    this_ahv = data.(cfg_param.target_session).ahv_ss;
    this_sdf = data.(cfg_param.target_session).(cfg_param.data_to_fit);

    cfg_param.nKeep = 10; % number of particles to keep for optimization after initial sweep
    cfg_param.addRandom = 0; % if 1, also choose random particles to keep for optimization
    predfun = @(xtrain,xtest)HDerrfun_xval_2D(xtrain,xtest,particle_vals,this_ahv,this_sdf,cfg_param.tc,cfg_param);
    
    %xtrain = data.(cfg_param.target_session).ahv_ss.tvec(1:end-500)';
    %xtest = data.(cfg_param.target_session).ahv_ss.tvec(end-499:end)';
    %err = predfun(xtrain(~isnan(xtrain)),xtest(~isnan(xtest)));
    
    opts = statset; opts.UseParallel = 'true';
    
    param = []; param_count = 1; param_hist = [];
    err(iM,:) = crossval(predfun,this_ahv.tvec','kfold',10,'mcreps',1,'options',opts);
    
    fprintf('%s %s, err %.2f +/- %.2f\n',cfg_param.target_session,model(iM).label,nanmean(err(iM,:)),nanstd(err(iM,:)));
    all_param{iM} = param;
    all_param_hist{iM} = param_hist;
    
    % now compute expected HD PFD changes (error) for gain and drift
    % estimates
    cfg_hd = []; cfg_hd.hd0 = ref_hd0;
    [~,true_hd_unwrapped] = AHVtoHD(cfg_hd,this_ahv);
    clear this_err_drift this_err_gain this_err_both;
    for iFold = 1:size(param,1) % param struct holds 
        
        cfg_hd = []; cfg_hd.hd0 = param(iFold,4);
        cfg_hd.gain = [param(iFold,1) param(iFold,2)];
        [~,this_hd_unwrapped] = AHVtoHD(cfg_hd,this_ahv);
        this_err_gain(iFold) = true_hd_unwrapped(end)-this_hd_unwrapped(end);
        
        cfg_hd = []; cfg_hd.hd0 = param(iFold,4);
        cfg_hd.drift = param(iFold,3);
        [~,this_hd_unwrapped] = AHVtoHD(cfg_hd,this_ahv);
        this_err_drift(iFold) = true_hd_unwrapped(end)-this_hd_unwrapped(end);
        
        cfg_hd = []; cfg_hd.hd0 = param(iFold,4);
        cfg_hd.gain = [param(iFold,1) param(iFold,2)];
        [~,this_hd_unwrapped] = AHVtoHD(cfg_hd,this_ahv);
        this_err_gain(iFold) = true_hd_unwrapped(end)-this_hd_unwrapped(end);
        
        cfg_hd = []; cfg_hd.hd0 = param(iFold,4);
        cfg_hd.gain = [param(iFold,1) param(iFold,2)];
        cfg_hd.drift = param(iFold,3);
        [~,this_hd_unwrapped] = AHVtoHD(cfg_hd,this_ahv);
        this_err_both(iFold) = true_hd_unwrapped(end)-this_hd_unwrapped(end);
        
    end
    all_wraperr{iM} = cat(1,this_err_gain,this_err_drift,this_err_both);
    all_wraperr{iM} = nanmean(all_wraperr{iM},2); % avg over folds
    
    % compute expected tuning curves for first cell, model 4
    if iM == 4
        % first, get internal HD
        this_param = nanmean(param); % note will not work for angle (hd0)
        cfg_hd = []; cfg_hd.hd0 = this_param(4);
        cfg_hd.gain = [this_param(1) this_param(2)];
        cfg_hd.drift = this_param(3);
        internal_hd = AHVtoHD(cfg_hd,this_ahv);
        
        % get sdf
        this_sim_sdf = GenerateSDFfromTC([],internal_hd,data.(cfg_param.target_session).tc);
        
        % now get TCs
        tc = TuningCurvesSDF(cfg_tc,this_sim_sdf,data.(cfg_param.target_session).hd_ss);
        tc.xbin = cfg_tc.bin_centers;
        out.tc_fit = tc;
        
        % segmented
        for iC = 1:nCells

            seg = linspace(data.(fnames{iF}).sdf.tvec(1),data.(fnames{iF}).sdf.tvec(end),cfg_tc.nseg+1);
            
            for iSeg = 1:cfg_tc.nseg
                this_sdf = restrict(this_sim_sdf,seg(iSeg),seg(iSeg+1));
                this_hd = restrict(data.(cfg_param.target_session).hd_ss,seg(iSeg),seg(iSeg+1));
                this_tc = TuningCurvesSDF(cfg_tc,this_sdf,this_hd);
                
                out.tc_fit_seg{iSeg} = this_tc.tc;
            end
            
        end
        
        
        % keep track of originals
        out.tc_actual = data.(cfg_param.target_session).tc;
        out.tc_actual_seg = data.(cfg_param.target_session).tc_seg;
        
    end
    
end % of loop over models

out.err = err;
out.wraperr = all_wraperr;
out.param = all_param;
out.param_hist = all_param_hist; % this has two layers of cells -- different models, and dimensions
out.model = model;
out.target_fn = data.(cfg_param.target_session).fn;
out.histBinCenters = cfg_param.histBins;