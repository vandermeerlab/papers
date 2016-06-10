function data = GetSDFfromHDdata(cfg_in,data)
% function data = GetSDFfromHDdata(cfg_in,data)
%
%

cfg_def = [];
cfg_def.dt = 1/60;
cfg_def.subsample_factor = 7;

cfg = ProcessConfig(cfg_def,cfg_in);

fnames = fieldnames(data);

for iF = 1:length(fnames)
    
    this_data = data.(fnames{iF});

    nCells = size(this_data.obs_fr,1);
    tvec_orig = this_data.obs_hd.tvec;
    tvec_interp = this_data.hd.tvec;
    
    % first, need to get firing rates on full timebase (tvec_interp)
    [~,~,interp_idx] = intersect(tvec_orig,tvec_interp); % compare original tvec with interpolated one
    clear obs_fr_interp;
    for iC = 1:nCells
        obs_fr_interp(iC,:) = nan(size(this_data.hd.tvec));
        obs_fr_interp(iC,interp_idx) = this_data.obs_fr(iC,:);
    end
    
    % convolve spikes
    k = gausskernel(61,cfg.subsample_factor);
    for iC = 1:nCells
        % first, replace nans with zeros so can smooth
        nan_idx = find(isnan(obs_fr_interp(iC,:)));
        obs_fr_interp(iC,nan_idx) = 0;
        obs_fr_interp(iC,:) = conv(obs_fr_interp(iC,:),k,'same')./cfg.dt;
        obs_fr_interp(iC,nan_idx) = NaN;
    end
    
    % get SDFs
    sdf_data = []; sdf_data_ss = [];
    for iC = 1:nCells
       sdf_data = cat(1,sdf_data,obs_fr_interp(iC,:));
       sdf_data_ss = cat(1,sdf_data_ss,obs_fr_interp(iC,1:cfg.subsample_factor:end));
    end
    
    sdf = tsd(tvec_interp,sdf_data);
    sdf_ss = tsd(tvec_interp(1:cfg.subsample_factor:end),sdf_data_ss);
 
    % update
    data.(fnames{iF}).sdf = sdf;
    data.(fnames{iF}).sdf_ss = sdf_ss;
    
end