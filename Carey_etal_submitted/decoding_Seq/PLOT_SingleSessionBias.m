function fh = PLOT_SingleSessionBias(cfg_in,data)
% function fh = PLOT_SingleSessionBias(cfg_in,data)
%
% plots behavior and SWR bias across sessions
% (MotivationalT data set)
%
% INPUTS:
%
% cfg_in: options that control display properties
%
% data: struct with .all, .pre, .task, .post fields containing data output
% by ALL_CollectDecSeq.m, see ALL_Plot_DecSeq for an example of how to load
% these
%
% OUTPUTS:
%


% TODO: break out actual plotter into function, needs to handle multiple
% data points (for all rats) and single data points (for indiv rats)

cfg_def = [];

cfg = ProcessConfig(cfg_def,cfg_in);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALL SEQUENCES, AVERAGED ACROSS SUBJECTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.output_fn = cat(2,cfg.outbasefn,'sessionProps_overall');

fh{1} = figure;
subplot(4,3,1);

% before anything else, draw background
for iC = 1:6 % we have 6 sessions
    
    h(iC) = rectangle('Position',[iC-0.5 0 1 1]);
    if mod(iC,2) % odd
        set(h(iC),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
    else % even
        set(h(iC),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
    end
    
end
hold on;
plot([0.5 6.5],[0.5 0.5],'LineStyle','--','LineWidth',1,'Color',[1 1 1]);

% now collect data
this_data = data.all.data.all.ALL_sig_seq;

sess_label = repmat(1:6,[1 4]);
rat_label = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));

if max(this_data.sess) == 19 % collector was run using 19 included sessions only (instead of all 24)
    fprintf('PLOT_SingleSessionBias: manual session IDs used.\n');
    this_data.sess = [2 2 3 3 4 4 5 5 6 6 10 10 11 11 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24];   
end

sess_idx = sess_label(this_data.sess);
rat_idx = rat_label(this_data.sess);

% need to swap water and food sessions for rats 3 and 4, which started with food not water
swp = [2 1 4 3 6 5];
for iRat = 3:4
   this_rat_idx = find(rat_idx == iRat);
   sess_idx(this_rat_idx) = swp(sess_idx(this_rat_idx));
end

% normalized (N) input data is redundant for left and right, so remove
sess_idx = sess_idx(1:2:end); rat_idx = rat_idx(1:2:end);
seq = this_data.countN(1:2:end);
behav = this_data.allTrialsN(1:2:end);
choice = this_data.choiceN(1:2:end);
type = this_data.type(1:2:end);

% get averages for each session, and plot each data point individually
for iSess = 1:6
   sess_out(iSess) = iSess;
   seq_out(iSess) = nanmean(seq(sess_idx == iSess));
   behav_out(iSess) = nanmean(behav(sess_idx == iSess));
   
   plot(iSess-cfg.offs,seq(sess_idx == iSess),'.','MarkerSize',5,'Color',[0 0.7 0]);
   plot(iSess+cfg.offs,behav(sess_idx == iSess),'.','MarkerSize',5,'Color',[0 0 0]);
end

plot([0.5 6.5],[0.5 0.5],'LineStyle','--','LineWidth',1,'Color',[1 1 1]);
hold on;
set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
box off;

yyaxis right; set(gca,'YColor',[0 0.7 0],'LineWidth',1,'FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'YLim',[0 1]);
ylabel('SWR sequence p(food)');


plot(sess_out+cfg.offs,behav_out,'k','LineStyle','-');
plot(sess_out+cfg.offs,behav_out,'.k','MarkerSize',20);

plot(sess_out-cfg.offs,seq_out,'Color',[0 0.7 0],'LineStyle','-');
plot(sess_out-cfg.offs,seq_out,'.','MarkerSize',20,'Color',[0 0.7 0]);

% make nice
set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
box off;

yyaxis left;
ylabel('choice p(food)');

if cfg.writeOutput
    maximize; drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc',[cfg.output_fn,'.eps']);
    cd(cfg.originalFolder)
end

% print some stats to go with this figure
fprintf('likelihood ratio test for model comparison with and without motiv state:\n');
tbl = table(categorical(rat_idx)',categorical(sess_idx)',behav,categorical(type),seq,choice,'VariableNames',{'Subject','Session','Behav','MotType','SeqContent','Choice'});
lme = fitglme(tbl,'SeqContent ~ MotType + (1|Subject)');
lme2 = fitglme(tbl,'SeqContent ~ 1 + (1|Subject)');
compare(lme2,lme)

fprintf('Pearson correlation between choice and sequence content bias:\n');
[r,p] = corrcoef(tbl.Choice,tbl.SeqContent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BREAK OUT BY SUBJECT AND PHASE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.output_fn = cat(2,cfg.outbasefn,'sessionProps_byPhase');

fh{2} = figure;

for iRat = 1:length(cfg.rats)
    for iW = 1:length(cfg.what)
        
        this_data = eval(cat(2,'data.',cfg.what{iW},'.data.',cfg.rats{iRat},'.ALL_sig_seq'));

        switch iRat % collector was run using 19 included sessions only (instead of all 24)
            case 2
                this_data.sess = [4 4 5 5];
            case {3, 4}
                this_data.sess = [1 1 2 2 3 3 4 4 5 5 6 6];
             
        end
        
        subplot(4,3,(iRat-1)*3+iW);
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        box off;
        
        yyaxis right; set(gca,'YColor',[0 0.7 0],'LineWidth',1,'FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'YLim',[0 1]);
      
        for iC = 1:6
            
            h(iC) = rectangle('Position',[iC-0.5 0 1 1]);
            if mod(iC,2) % odd
                if iRat >= 3
                    set(h(iC),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                else
                    set(h(iC),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                end
            else % even
                if iRat >= 3
                    set(h(iC),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                else
                    set(h(iC),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                end
            end
            
        end
        hold on;
        plot([0.5 6.5],[0.5 0.5],'LineStyle','--','LineWidth',1,'Color',[1 1 1]);
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        
        plot(this_data.sess(1:2:end),this_data.allTrialsN(1:2:end),'k','LineStyle','-');
        plot(this_data.sess(1:2:end),this_data.allTrialsN(1:2:end),'.k','MarkerSize',20);
        
        plot(this_data.sess(1:2:end),this_data.countN(1:2:end),'Color',[0 0.7 0],'LineStyle','-');
        hold on;
        plot(this_data.sess(1:2:end),this_data.countN(1:2:end),'.','MarkerSize',20,'Color',[0 0.7 0]);
               
        % make nice
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',cfg.fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        box off;
          
        if iRat == 1
            title(cfg.what{iW});
        end
        
    end
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'_supp.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'_supp.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'_supp.eps']);
    cd(cfg.originalFolder)
end