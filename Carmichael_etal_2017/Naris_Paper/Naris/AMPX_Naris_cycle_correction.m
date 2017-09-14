
session_type = 'pre';
for iSess  = 1:length(Session_list)
    data_out = all_data.(strrep(Session_list{iSess}, '-', '_'));
    %% smoothed plane fitting
    cfg = [];
    cfg.name = strrep(Session_list{iSess}, '-', '_');
    cfg.session_type = session_type;
    cfg.plot = 0;
    %low gamma
    disp('low gamma')
    cfg.type = 'low';
    data_out.lg.power.smooth.stats = AMPX_get_plane_fitting(cfg,data_out.lg.power);
    plane_stats.low_gamma.smooth.rsq =data_out.lg.power.smooth.stats.rsq;
    plane_stats.low_gamma.smooth.rsq2 =data_out.lg.power.smooth.stats.rsq2;
    
    % high gamma
    disp('high gamma')
    cfg.type = 'high';
    data_out.hg.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.hg.power);
    plane_stats.high_gamma.smooth.rsq =data_out.hg.power.smooth.stats.rsq;
    plane_stats.high_gamma.smooth.rsq2 =data_out.hg.power.smooth.stats.rsq2;
    %
    % low control
    disp('low control')
    cfg.type = 'low_ctrl';
    data_out.lg_ran.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.lg_ran.power);
    plane_stats.low_control.smooth.rsq =data_out.lg_ran.power.smooth.stats.rsq;
    plane_stats.low_control.smooth.rsq2 =data_out.lg_ran.power.smooth.stats.rsq2;
    
    % high control
    disp('high control')
    cfg.type = 'high_ctrl';
    data_out.hg_ran.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.hg_ran.power);
    plane_stats.high_control.smooth.rsq =data_out.hg_ran.power.smooth.stats.rsq;
    plane_stats.high_control.smooth.rsq2 =data_out.hg_ran.power.smooth.stats.rsq2;
    
    all_data.(strrep(Session_list{iSess}, '-', '_')) = data_out;
end