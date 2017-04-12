cfg_def = [];
cfg_def.datapath = 'D:\My_Documents\Dropbox\projects\HDfit\data';
cfg_def.session = 3; % which session to load
cfg_def.mode = 'kalmanwrapped'; % 'kalmanwrapped', 'smooth', 'kalman'
cfg_def.debug = 1;
cfg_def.target_session = 'std'; % 'std', 'laser'

cfg_master = ProcessConfig(cfg_def,[]);


%% cfg
cfg_master.subsample_factor = 1; % 7
%cfg_master.smoothwin = 11; % samples of 1/60


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

%% show ahv
figure;
h(1) = plot(this_data.ahv_ss.tvec,this_data.ahv_ss.data,'.k');
set(gca,'XLim',[0 10],'FontSize',24,'TickDir','out','LineWidth',1,'XTick',0:5:10,'YLim',[-400 400]); box off;
xlabel('time (s)'); ylabel('AHV (deg/s)');

%% show example HD profiles
this_data = data.std;

cfg_test = [];
cfg_test.hd0 = this_data.hd_ss.data(1);

true_hd = AHVtoHDfast(this_data.ahv_ss,[1 1 0 cfg_test.hd0 0]);
gain = AHVtoHDfast(this_data.ahv_ss,[1.5 1 0 cfg_test.hd0 0]);
drift = AHVtoHDfast(this_data.ahv_ss,[1 1 -10 cfg_test.hd0 0]);

figure;
h(1) = plot(this_data.ahv_ss.tvec,true_hd,'.k'); hold on;
h(2) = plot(this_data.ahv_ss.tvec,gain,'.r');
h(3) = plot(this_data.ahv_ss.tvec,drift,'.b');

set(gca,'XLim',[0 10],'FontSize',24,'TickDir','out','LineWidth',1,'YLim',[0 360],'YTick',0:90:360); box off;
legend(h,{'true HD','gain change','drift change'},'Location','Northwest'); legend boxoff;
xlabel('time (s)'); ylabel('HD (deg)');

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

iF = 1;

tc = TuningCurvesSDF(cfg_tc,data.(fnames{iF}).sdf,data.(fnames{iF}).hd);
tc.xbin = cfg_tc.bin_centers;

data.(fnames{iF}).tc = tc;

figure;
plot(tc.xbin,tc.tc'); 
set(gca,'FontSize',24,'TickDir','out','LineWidth',1,'YLim',[0 80],'XLim',[0 360],'XTick',0:90:360); box off;

cmap = colormap(jet);
