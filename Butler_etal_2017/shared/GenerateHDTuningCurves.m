function tc_out = GenerateHDTuningCurves(cfg_in)
%


cfg_def = [];
cfg_def.binsize = 1;
cfg_def.pfd = 180;
cfg_def.sd = 30;
cfg_def.maxfr = 40;

cfg = ProcessConfig(cfg_def,cfg_in);

%
nCells = length(cfg.pfd);

bin_edges = 0:cfg.binsize:360;
bin_centers = cfg.binsize/2+bin_edges(1:end-1);

for iC = 1:nCells
    
    % first, get gaussian centered at 180
    this_tc = cfg.maxfr(iC)*exp((-(bin_centers-180).^2)/(2*cfg.sd(iC).^2));
    
    % then shift to correct PFD
    shift = round((cfg.pfd(iC) - 180)/cfg.binsize);
    tc(iC,:) = circshift(this_tc,[1 shift]);
    
end

tc_out.tc = tc;
tc_out.xbin = bin_centers;