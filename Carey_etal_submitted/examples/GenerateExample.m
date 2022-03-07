%% setup
SET_GitHub_root = 'C:\Users\mvdm\Documents\GitHub'; % replace this with the location of your local GitHub root folder
SET_data_fd = 'D:\data\R042\R042-2013-08-18'; % replace this with the location of your local data folder

addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\shared')));

%% load the data
cd(SET_data_fd);
 
please = []; please.load_questionable_cells = 1;
S = LoadSpikes(please); % `please` variable overrides default LoadSpikes() options
 
pos = LoadPos([]); % empty input [] causes LoadPos() to use default options

LoadExpKeys; % annotation file containing some basic information about this data session
LoadMetadata; % loads experimenter-annotated file associated with each data session

%% arrange data
clear expCond;

expCond(1).label = 'left'; % this is a T-maze, we are interested in 'left' and 'right' trials
expCond(2).label = 'right'; % these are just labels we can make up here to keep track of which condition means what

expCond(1).t = metadata.taskvars.trial_iv_L; % previously stored trial start and end times for left trials
expCond(2).t = metadata.taskvars.trial_iv_R; 

expCond(1).coord = metadata.coord.coordL; % previously user input idealized linear track
expCond(2).coord = metadata.coord.coordR; % note, this is in units of "camera pixels", not cm

expCond(1).S = S;
expCond(2).S = S;

%% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
nCond = length(expCond);
for iCond = 1:nCond
   
    cfg_linpos = []; cfg_linpos.Coord = expCond(iCond).coord;
    expCond(iCond).linpos = LinearizePos(cfg_linpos,pos);
   
end

%% find intervals where rat is running
spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; 
run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 10 pix/s

%% restrict (linearized) position data and spike data to desired intervals
for iCond = 1:nCond
   
    fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
end

%% get tuning curves, see lab wiki at:
% http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week12
for iCond = 1:nCond
    
    cfg_tc = [];
    expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos); 
    
end

%% detect place fields
for iCond = 1:nCond
    
    expCond(iCond).fields = DetectPlaceCells1D([],expCond(iCond).tc.tc);
    
end
    
%% load CSC with good SWRs
cfg = []; cfg.fc = ExpKeys.goodSWR(1);
lfp = LoadCSC(cfg); lfp = decimate_tsd([],lfp);

LoadCandidates;

%% decode with moving window
cfg = []; cfg.keepPosterior = 1;
out = Generate_DecSeq(cfg);

%% visualize with decoding
ttl = {'left','right'}; 
xl = [5042.3 5043]; % R042-2013-08-18 (uncomment lines below for different sessions)
%xl = [5809.9 5810.6]; % R050-2014-04-03
%xl = [4942.4 4943.1]; % R044-2013-12-21
%xl = [5838.6 5839.3]; % R064-2015-04-21
for iCond = 1:nCond
    % raster
    ax(iCond) = subplot(2,2,iCond);
    S_temp = SelectTS([],out.expCond(iCond).decS,out.expCond(iCond).fields.template_idx);
    cfg_plot = []; 
    cfg_plot.lfp(1) = lfp;
    cfg_plot.evt = out.expCond(iCond).seq_iv;
    cfg_plot.openNewFig = 0;
    MultiRaster(cfg_plot,S_temp);
    xlim(xl);
    set(gca,'YTick',[],'TickDir','out','XTick',[]);
    title(ttl{iCond}); xlabel(''); ylabel('');
    
    % decoded posterior
    ax(iCond+2) = subplot(2,2,iCond+2);
    imagesc(out.expCond(iCond).P2.tvec,1:111,out.expCond(iCond).P2.data);
    set(gca,'YDir','normal','YTick',[],'XTick',xl(1):0.2:xl(2),'TickDir','out');
    xlabel('time (s)');
    caxis([0 0.3]);
    xlim(xl);
end
linkaxes(ax,'x')

%% export
print(gcf,'-painters','-depsc','example.eps');
print(gcf,'-painters','-dpng','-r300','example.png');