function data = WillPreProcData(cfg_in,data)
% function data = WillPreProcData(cfg_in,data)
%
% preprocess HD data

cfg_def = [];
cfg_def.dt = 1/60;
cfg_def.smoothwin = 11;
cfg_def.subsample_factor = 7;
cfg_def.debug = 0;
cfg_def.mode = 'kalman'; % 'smooth', 'kalman', 'kalmanwrapped'

cfg = ProcessConfig(cfg_def,cfg_in);


fnames = fieldnames(data);
for iF = 1:length(fnames)
   
    this_data = data.(fnames{iF});
    
    %% need to unwrap HD to interpolate missing data
    [~,hd_unwrapped] = HDtoAHV(this_data.obs_hd);
    
    tvec_interp = (0:max(this_data.bin_idx)-1)*cfg.dt;
    
    switch cfg.mode
        case 'smooth'
            hd_unwrapped_interp = interp1(hd_unwrapped.tvec,hd_unwrapped.data,tvec_interp,'linear');
            hd_unwrapped_interp = medfilt1(hd_unwrapped_interp,cfg.smoothwin);
            hd_unwrapped_interp = smooth(hd_unwrapped_interp,cfg.smoothwin); % why does this take so long?
            
            hd_unwrapped_interp = tsd(tvec_interp,hd_unwrapped_interp);
            hd = wrapHD(hd_unwrapped_interp.data)';
            
        case 'kalman'
            hd_raw = nan(size(tvec_interp));
            [~,keep_idx,~] = intersect(tvec_interp,hd_unwrapped.tvec);
            hd_raw(keep_idx) = hd_unwrapped.data;
            
            hd_est = generic_kalman([],hd_raw);
            hd_estR = generic_kalman([],hd_raw(end:-1:1));
            
            hd_unwrapped_interp = nanmean(cat(1,hd_est,hd_estR(end:-1:1)))';
            
            hd_unwrapped_interp = tsd(tvec_interp,hd_unwrapped_interp);
            hd = wrapHD(hd_unwrapped_interp.data)';
            
        case 'kalmanwrapped'
            hd_raw = nan(size(tvec_interp));
            [~,keep_idx,~] = intersect(tvec_interp,hd_unwrapped.tvec);
            hd_raw(keep_idx) = this_data.obs_hd.data;
            x_obs = (hd_raw-180)/(180/pi); 
            dt = 1/60;
            
            % initial state
            v0 = diffang_rad(nanmedian(x_obs(24:37)),nanmedian(x_obs(1:10)))*2; % starting velocity
            if isnan(v0), v0 = 0; end
            x0 = [nanmedian(x_obs(1:10)); v0; 0]; % assumes 60 samples/s
            
            % state transition matrix
            A = [1 dt dt^2/2;
                0  1 dt;
                0  0  1];
            
            % state covariance
            Cx = [0.25 0 0;
                0 0.5 0;
                0 0 1];
            
            cy = 5;
            
            z = wkf(x_obs,A,Cx,cy,x0); % forward pass

            v0 = diffang_rad(x_obs(end),nanmedian(x_obs(end-36:end-24)))*2; % starting velocity
            if isnan(v0), v0 = 0; end
            
            x0 = [nanmedian(x_obs(end:-1:end-9)); v0; 0]; % nanmedian to avoid NaNs...
            zRev = wkf(x_obs(end:-1:1),A,Cx,cy,x0); % reverse pass
            
            zRev = zRev(1,end:-1:1);
            
            clear hd;
            for iS = length(z):-1:1
                hd(iS) = circ_mean([z(1,iS),zRev(1,iS)]');
            end
            hd = hd*(180/pi)+180;
            
        otherwise
            error('Unknown mode.');
            
            
    end
    
    hd = tsd(tvec_interp,hd);
    ahv = HDtoAHV(hd);
    
    if cfg.debug
        %% could consider kalman filter instead of smoothing + linear interpolation... but this looks ok for now
        figure;
        plot(hd,'.r'); hold on;
        plot(this_data.obs_hd,'.k');
    end
    
    data.(fnames{iF}).hd = hd;
 
    %% subsample
    dt_ss = cfg.dt*cfg.subsample_factor;
    hd_ss = hd;
    hd_ss.tvec = hd_ss.tvec(1:cfg.subsample_factor:end);
    hd_ss.data = hd_ss.data(1:cfg.subsample_factor:end);
    
    %% get AHV
    ahv_ss = HDtoAHV(hd_ss);
    
    data.(fnames{iF}).ahv_ss = ahv_ss;
    data.(fnames{iF}).ahv = ahv;
    data.(fnames{iF}).hd_ss = hd_ss;
    
end