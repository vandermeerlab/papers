function at_feeder = Julien_atFeederTimes()

pos = LoadPos([]);

load(FindFile('*Coord.mat'));
load(FindFile('*IndMarker.mat'));

feeder_x = Coord(1,IndMarker.Feeders(2));
feeder_y = Coord(2,IndMarker.Feeders(2));

feeder_dist = sqrt((getd(pos,'x')-feeder_x).^2 + (getd(pos,'y')-feeder_y).^2);
feeder_dist = tsd(pos.tvec,feeder_dist);

% now use TSDtoIV
cfg = [];
cfg.dcn = '<';
cfg.method = 'raw';
cfg.threshold = 20;
at_feeder = TSDtoIV(cfg,feeder_dist);

%cfg.plot = [];
%cfg_plot.display = 'tsd'; % tsd, iv
%PlotTSDfromIV(cfg_plot,at_feeder,feeder_dist);

