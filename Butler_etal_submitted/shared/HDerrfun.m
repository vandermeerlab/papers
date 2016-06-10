function err = HDerrfun(params,obs_ahv,obs_sdf,tc)
%

internal_hd = AHVtoHDfast(obs_ahv,params);
internal_hd = tsd(obs_ahv.tvec,internal_hd);

cfg = []; cfg.mode = 'interp';
predicted_sdf = GenerateSDFfromTC(cfg,internal_hd,tc);

err = (predicted_sdf.data-obs_sdf.data).^2;
err = nanmean(err(:));


