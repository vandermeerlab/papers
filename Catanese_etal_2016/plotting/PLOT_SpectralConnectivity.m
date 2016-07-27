%% PLOT_SpectralConnectivity.m
% 
% Julien Catanese & Matthijs van der Meer

%% set paths
restoredefaultpath;
cd('D:\My_Documents\GitHub\fieldtrip');
ft_defaults;

rmpath('D:\My_Documents\GitHub\fieldtrip\external\signal\');

addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('D:\My_Documents\GitHub\vandermeerlab\code-matlab\tasks\Julien_linear_track')); % Detect events, CountCycles live here

%% load gamma events (use MASTER_CollectGammaEvents.m to obtain) -- puts ALL_evt variable in workspace
cd('D:\My_Documents\Dropbox\projects\Julien_multiLFP\2016-01-07');
%load(FindFile('gamma*conn.mat'));

%% parameters
PARAM_plotcolor = 'b';

%% define what to run
rats = {'R026','R032','R033','R039'};
available_rats = fieldnames(ALL_evt);

nSession = 0;
for iRat = 1:length(rats)
    
    this_rat = rats{iRat};
    
    if ~strmatch(this_rat,available_rats)
       warning('Rat %s not available -- skipping...',rats{iRat});
       continue;
    end
    
    available_sessions = fieldnames(ALL_evt.(this_rat));
    
    for iSession = 1:length(available_sessions)
    
        nSession = nSession + 1;
        
        this_session = available_sessions{iSession};
        this_session_data = ALL_evt.(this_rat).(this_session);
        
        fprintf('\n\nProcessing session %s...\n',this_session);
        
        % plot indiv session psi
        figure(1);
        subplot(3,3,nSession);
        
        plot(this_session_data.psi.freq,this_session_data.psi.psispctrm);
        hold on,  plot(this_session_data.psi.freq,this_session_data.psi.psispctrm + this_session_data.psi.psispctrmsem,':r');
        hold on,  plot(this_session_data.psi.freq,this_session_data.psi.psispctrm - this_session_data.psi.psispctrmsem,':r');
        hold on,  plot(this_session_data.psi.freq,0, '--k');
        xlim([30 105]); set(gca,'FontSize',14);
        
        title(sprintf('psiraw, %s',this_session),'FontSize',14,'Interpreter','none');
        
        psi_all_raw(nSession,:) = this_session_data.psi.psispctrm;
        psi_all(nSession,:) = this_session_data.psi.psispctrm./this_session_data.psi.psispctrmsem;
        
        % plot indiv session icoh
        figure(2);
        subplot(3,3,nSession);
        
        plot(this_session_data.coh.cohfreq,this_session_data.coh.cohspctrm);
        hold on,  plot(this_session_data.coh.cohfreq,0, '--k');
        xlim([30 105]); set(gca,'FontSize',14);
        
        title(sprintf('icoh, %s',this_session),'FontSize',14,'Interpreter','none');
        
        coh_all_raw(iSession,:) = this_session_data.coh.cohspctrm; % convert to deg
        
    end % over sessions
    
end % over rats

% mean psi
figure(1)
subplot(339)
plot(this_session_data.psi.freq,nanmean(psi_all),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  plot(this_session_data.psi.freq,0, '--k');
xlim([30 105]); set(gca,'FontSize',14,'LineWidth',1); box off;
title(sprintf('psi average',this_session),'FontSize',14,'Interpreter','none');

% mean icoh
figure(2)
subplot(339)
plot(this_session_data.coh.cohfreq,nanmean(coh_all_raw),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  plot(this_session_data.psi.freq,0, '--k');
xlim([30 105]); set(gca,'FontSize',14,'LineWidth',1); box off;
title(sprintf('icoh average',this_session),'FontSize',14,'Interpreter','none');

%% single session example (for Figure 5 -- used ALL data)
figure(3)
iRat = 3; iSession = 1;
this_rat = rats{iRat}; 
available_sessions = fieldnames(ALL_evt.(this_rat));
this_session = available_sessions{iSession};
this_session_data = ALL_evt.(this_rat).(this_session);

fvec = this_session_data.psi.freq;
psl = this_session_data.psi.psispctrm;
psl_sem = this_session_data.psi.psispctrmsem;
psi = this_session_data.psi.psispctrm./this_session_data.psi.psispctrmsem;
icoh = this_session_data.coh.cohspctrm;

lgx = [45 45 65 65]; lgc = [0.8 1 0.8];
hgx = [70 70 90 90]; hgc = [1 0.8 0.8];
fa = 1;

subplot(321);
plot(fvec,icoh,'g','LineWidth',2);
box off; hold on;
plot([35 105],[0 0],'k--','LineWidth',1);
set(gca,'LineWidth',1,'XLim',[40 100],'XTick',40:20:100,'YLim',[-20 0],'YTick',-20:5:0);
ylabel('phase lag \phi_{vStr-mPFC} (\circ)');
yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))

subplot(323);
plot(fvec,psl,'c','LineWidth',2);
box off; hold on;
plot(fvec,psl+psl_sem,'c--','LineWidth',1);
plot(fvec,psl-psl_sem,'c--','LineWidth',1);
plot([35 105],[0 0],'k--','LineWidth',1);
set(gca,'LineWidth',1,'XLim',[40 100],'XTick',40:20:100,'YLim',[-1 1],'YTick',-1:0.5:1);
ylabel('phase slope (\circ/Hz)');
yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))

subplot(325);
plot(fvec,psi,'k','LineWidth',2);
box off; hold on;
plot([35 105],[0 0],'k--','LineWidth',1);
set(gca,'LineWidth',1,'XLim',[40 100],'XTick',40:20:100,'YLim',[-15 15],'YTick',-15:5:15);
ylabel('phase slope index (\sigma)');
xlabel('Frequency (Hz)');
yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))

%% more examples (Figure 6) - TASK data
figure;
h = subplot(431);
iRat = 1; iSession = 2;
this_rat = rats{iRat}; 
available_sessions = fieldnames(ALL_evt.(this_rat));
this_session = available_sessions{iSession};
this_session_data = ALL_evt.(this_rat).(this_session);

FUNC_plotPSI(h,this_session_data);

h = subplot(434);
iRat = 2; iSession = 1;
this_rat = rats{iRat}; 
available_sessions = fieldnames(ALL_evt.(this_rat));
this_session = available_sessions{iSession};
this_session_data = ALL_evt.(this_rat).(this_session);

FUNC_plotPSI(h,this_session_data);

h = subplot(437);
iRat = 3; iSession = 2; % note different example
this_rat = rats{iRat}; 
available_sessions = fieldnames(ALL_evt.(this_rat));
this_session = available_sessions{iSession};
this_session_data = ALL_evt.(this_rat).(this_session);

FUNC_plotPSI(h,this_session_data);

h = subplot(4,3,10);
iRat = 4; iSession = 1;
this_rat = rats{iRat}; 
available_sessions = fieldnames(ALL_evt.(this_rat));
this_session = available_sessions{iSession};
this_session_data = ALL_evt.(this_rat).(this_session);

FUNC_plotPSI(h,this_session_data);

%% grand average PSI
PARAM_plotcolor = 'k';
figure;
subplot(221);
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20)),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))-nanstd(psi_all(:,1:end-20))./2,'Color',[0.7 0.7 0.7],'LineWidth',1);
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))+nanstd(psi_all(:,1:end-20))./2,'Color',[0.7 0.7 0.7],'LineWidth',1);

plot(this_session_data.psi.freq,0, '--k');
xlim([35 105]); set(gca,'FontSize',18,'LineWidth',1,'XTick',[40:20:100],'YTick',[-10:5:10],'YLim',[-10 10]); 
box off;
ylabel('mean phase slope index');

lgx = [45 45 65 65]; lgc = [0.8 1 0.8];
hgx = [70 70 90 90]; hgc = [1 0.8 0.8];
fa = 1;
yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))

%% task vs. rest average PSI
PARAM_plotcolor = 'r';
figure(10);
subplot(221);

% load task data
PARAM_plotcolor = 'r';
all_h(1) = plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20)),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))-nanstd(psi_all(:,1:end-20))./2,'Color',[1 0.7 0.7],'LineWidth',1);
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))+nanstd(psi_all(:,1:end-20))./2,'Color',[1 0.7 0.7],'LineWidth',1);

plot(this_session_data.psi.freq,0, '--k');
xlim([35 105]); set(gca,'FontSize',18,'LineWidth',1,'XTick',[40:20:100],'YTick',[-10:5:15]); 
box off;
ylabel('mean phase slope index');

% load rest next
figure(10);
subplot(221); hold on;
PARAM_plotcolor = 'b';
all_h(2) = plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20)),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))-nanstd(psi_all(:,1:end-20))./2,'Color',[0.7 0.7 1],'LineWidth',1);
plot(this_session_data.psi.freq(1:end-20),nanmean(psi_all(:,1:end-20))+nanstd(psi_all(:,1:end-20))./2,'Color',[0.7 0.7 1],'LineWidth',1);

plot(this_session_data.psi.freq,0, '--k');
xlim([35 105]); set(gca,'FontSize',18,'LineWidth',1,'XTick',[40:20:100],'YTick',[-10:5:15]); 
box off;

legend(all_h,{'task','rest'}); legend boxoff;

lgx = [45 45 65 65]; lgc = [0.8 1 0.8];
hgx = [70 70 90 90]; hgc = [1 0.8 0.8];
fa = 1;
yl = ylim; y = [yl yl(end:-1:1)];
p = patch(lgx,y,lgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
p = patch(hgx,y,hgc); set(p,'FaceAlpha',fa,'EdgeColor','none');
set(gca,'children',flipud(get(gca,'children')))

%% mean phase slope
figure(3)
subplot(221)
plot(this_session_data.psi.freq,nanmean(psi_all_raw),'Color',PARAM_plotcolor,'LineWidth',2);
hold on,  plot(this_session_data.psi.freq,0, '--k');
xlim([30 105]); set(gca,'FontSize',14,'LineWidth',1); box off;
title(sprintf('phase slope average',this_session),'FontSize',14,'Interpreter','none');

mf = this_session_data.psi.freq(this_session_data.psi.freq < 100);
m = nanmean(psi_all_raw(:,this_session_data.psi.freq < 100));

m = m./360; m = m * 1000;

[mmax,max_idx] = max(m);
[mmin,min_idx] = min(m);
fprintf('avg max slope %1.3e (at %dHz; %.2f ms), min slope %1.3e (at %dHz; %.2f ms)\n',...
    mmax,mf(max_idx),mmax, ...
    mmin,mf(min_idx),mmin);