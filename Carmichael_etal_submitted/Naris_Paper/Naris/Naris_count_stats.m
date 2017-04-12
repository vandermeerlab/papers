function stats_out = Naris_count_stats(cfg_in, all_stats)
%%  Naris_count_stats : collects and compiles all occlusion count stats from the naris occlusion analysi pipeline
%   
% 
% 
%          Inputs: 
%           - cfg_in: [struct] contains configuration parameters
%           - all_stats: [struct] contains the low and high gamma counts
%           across all sessions for all rats.
%          Outputs: 
%           - stats: [struct] contains the compiled sount values for high
%           and low gamma across all subjects.  
% 
% EC - 2017-01-19
global PARAMS
cfg_def = [];

cfg = ProcessConfig2(cfg_def, cfg_in)
cfg.date = date; 
%% collect and compile all the count values across subjects/sessions
sess_list = fieldnames(all_stats);
phases = {'control', 'pre', 'ipsi', 'contra', 'post'};
for iPhase = 1:length(phases) % initialize the structrues
    low_gamma.(phases{iPhase}) = [];
    high_gamma.(phases{iPhase}) = [];
end

for iSess = 1:length(sess_list)
    for iPhase = 1:length(phases) % loop through the phases skipping the first one which is the contol. 
        low_gamma.(phases{iPhase})(iSess) = all_stats.(sess_list{iSess}).low.rate(iPhase);
        high_gamma.(phases{iPhase})(iSess) = all_stats.(sess_list{iSess}).high.rate(iPhase);
       
    end
end

%% Normalize to the "control"
low_gamma.ipsi = low_gamma.ipsi./low_gamma.control;
low_gamma.contra = low_gamma.contra./low_gamma.control;
low_gamma.control = low_gamma.control./low_gamma.control;

high_gamma.ipsi = high_gamma.ipsi./high_gamma.control;
high_gamma.contra = high_gamma.contra./high_gamma.control;
high_gamma.control = high_gamma.control./high_gamma.control;

%% Stats and output
ks = [];
for i = 1:4
    labels = {'ipsi', 'contra', 'control'};
    h = kstest(low_gamma.(labels{1}));
    if h ~=1
        disp('***************************************************************')
        disp(['KS test FAIL for ' low_gamma.(labels{1})])
        disp('***************************************************************')
        ks = [ks ; 1];
    end
    ks = [ks; 0];
end
%%

% ipsi = reshape(comp_data.ipsi,1,numel(comp_data.ipsi));
% control = reshape(comp_data.control,1,numel(comp_data.control));
% contra = reshape(comp_data.contra,1,numel(comp_data.contra));
if sum(ks)>=1
    [p_ip_con, h_ip_con] = signrank(low_gamma.ipsi, low_gamma.contra);
    [p_ip_ctr, h_ip_ctr] = signrank(low_gamma.ipsi, low_gamma.control);
    [p_con_ctr, h_con_ctr] = signrank(low_gamma.contra, low_gamma.control);
    
    [h_p_ip_con, h_h_ip_con] = signrank(high_gamma.ipsi, high_gamma.contra);
    [h_p_ip_ctr, h_h_ip_ctr] = signrank(high_gamma.ipsi, high_gamma.control);
    [h_p_con_ctr, h_h_con_ctr] = signrank(high_gamma.contra, high_gamma.control);
else
    fprintf('\n\nUsing T-Test\n\n')
    [h_ip_con, p_ip_con, ~, l_stats_ip_con] = ttest(low_gamma.ipsi, low_gamma.contra);
    [h_ip_ctr, p_ip_ctr, ~,l_stats_ip_ctr] = ttest(low_gamma.ipsi, low_gamma.control);
    [h_con_ctr, p_con_ctr, ~,l_stats_con_ctr] = ttest(low_gamma.contra, low_gamma.control);
    
    [h_h_ip_con, h_p_ip_con, ~, h_stats_ip_con] = ttest(high_gamma.ipsi, high_gamma.contra);
    [h_h_ip_ctr, h_p_ip_ctr, ~,h_stats_ip_ctr] = ttest(high_gamma.ipsi, high_gamma.control);
    [h_h_con_ctr, h_p_con_ctr, ~,h_stats_con_ctr] = ttest(high_gamma.contra, high_gamma.control);
end
%low SEM
low_pre_SEM = std(low_gamma.pre)/(sqrt(length(low_gamma.pre)));
low_ipsi_SEM = std(low_gamma.ipsi)/(sqrt(length(low_gamma.ipsi)));
low_contra_SEM = std(low_gamma.contra)/(sqrt(length(low_gamma.contra)));
low_post_SEM = std(low_gamma.post)/(sqrt(length(low_gamma.post)));
low_control_SEM = std(low_gamma.control)/(sqrt(length(low_gamma.control)));
% high SEM
high_pre_SEM = std(high_gamma.pre)/(sqrt(length(high_gamma.pre)));
high_ipsi_SEM = std(high_gamma.ipsi)/(sqrt(length(high_gamma.ipsi)));
high_contra_SEM = std(high_gamma.contra)/(sqrt(length(high_gamma.contra)));
high_post_SEM = std(high_gamma.post)/(sqrt(length(high_gamma.post)));
high_control_SEM = std(high_gamma.control)/(sqrt(length(high_gamma.control)));


%% print the stats 

% low gamma
fileID = fopen([PARAMS.stats_dir '\Naris_stats_count.txt'],'w+');
fprintf(fileID, ['\n_________________________________________\n'])
fprintf(fileID,datestr(date, 'yyyy-mm-dd-HH'))
fprintf(fileID, ['\nLow Gamma Event Count \n'])
fprintf(fileID,'Basic Stats\n')
fprintf(fileID,['Ipsilateral:     Mean ' num2str(mean(low_gamma.ipsi), '%4.4f') '   SD +/-' num2str(std(low_gamma.ipsi), '%4.4f') '   SEM +/-' num2str(std(low_gamma.ipsi)/sqrt(length(low_gamma.ipsi))) '\n'])
fprintf(fileID,['Contralateral:   Mean ' num2str(mean(low_gamma.contra),'%4.4f') '   SD +/-' num2str(std(low_gamma.contra), '%4.4f') '   SEM +/-' num2str(std(low_gamma.contra)/sqrt(length(low_gamma.contra))) '\n'])
fprintf(fileID,['Control:         Mean ' num2str(mean(low_gamma.control),'%4.4f') '   SD +/-' num2str(std(low_gamma.control), '%4.4f') '   SEM +/-' num2str(std(low_gamma.control)/sqrt(length(low_gamma.control))) '\n'])
if sum(ks) >=1
    fprintf(fileID,'\nWilcoxin Sign Rank test\n')
    fprintf(fileID,['Ipsilateral   vs. Contralateral:    P:' num2str(p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID,['Ipsilateral   vs. Control:          P:' num2str(p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID,['Contralateral vs. Control:          P:' num2str(p_con_ctr, '%4.4f') '\n' ])
else
    fprintf(fileID,'\nPaired T-Test\n')
    fprintf(fileID,['Ipsilateral   vs. Contralateral:   df(' num2str(l_stats_ip_con.df) ')   t:' num2str(l_stats_ip_con.tstat, '%4.4f') '  P:' num2str(p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID,['Ipsilateral   vs. Control:         df(' num2str(l_stats_ip_ctr.df) ')   t:' num2str(l_stats_ip_ctr.tstat, '%4.4f') '  P:' num2str(p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID,['Contralateral vs. Control:         df(' num2str(l_stats_con_ctr.df) ')   t:' num2str(l_stats_con_ctr.tstat, '%4.4f') '   P:' num2str(p_con_ctr, '%4.4f') '\n' ])
end

% high gamma
fprintf(fileID,'\n\nHigh Gamma\n')
fprintf(fileID,'Basic Stats\n')
fprintf(fileID,['Ipsilateral:     Mean ' num2str(mean(high_gamma.ipsi), '%4.4f') '   SD +/-' num2str(std(high_gamma.ipsi), '%4.4f') '   SEM +/-' num2str(std(high_gamma.ipsi)/sqrt(length(high_gamma.ipsi))) '\n'])
fprintf(fileID,['Contralateral:   Mean ' num2str(mean(high_gamma.contra),'%4.4f') '   SD +/-' num2str(std(high_gamma.contra), '%4.4f') '   SEM +/-' num2str(std(high_gamma.contra)/sqrt(length(high_gamma.contra))) '\n'])
fprintf(fileID,['Control:         Mean ' num2str(mean(high_gamma.control),'%4.4f') '   SD +/-' num2str(std(high_gamma.control), '%4.4f') '   SEM +/-' num2str(std(high_gamma.control)/sqrt(length(high_gamma.control))) '\n'])
if sum(ks) >=1
    fprintf(fileID,'\nWilcoxin Sign Rank test\n')
    fprintf(fileID,['Ipsilateral   vs. Contralateral:    P:' num2str(h_p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID,['Ipsilateral   vs. Control:          P:' num2str(h_p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID,['Contralateral vs. Control:          P:' num2str(h_p_con_ctr, '%4.4f') '\n' ])
else
    fprintf(fileID,'\nPaired T-Test\n')
    fprintf(fileID,['Ipsilateral   vs. Contralateral:   df(' num2str(h_stats_ip_con.df) ')   t:' num2str(h_stats_ip_con.tstat, '%4.4f') '  P:' num2str(h_p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID,['Ipsilateral   vs. Control:         df(' num2str(h_stats_ip_ctr.df) ')   t:' num2str(h_stats_ip_ctr.tstat, '%4.4f') '  P:' num2str(h_p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID,['Contralateral vs. Control:         df(' num2str(h_stats_con_ctr.df) ')   t:' num2str(h_stats_con_ctr.tstat, '%4.4f') '   P:' num2str(h_p_con_ctr, '%4.4f') '\n' ])
end

fclose(fileID);
% type 'G:\Naris\Naris_stats_count.txt'
% copyfile('G:\Naris\Naris_stats_count.txt','D:\Users\mvdmlab\My_Documents\GitHub\Papers_EC\vStr_Naris', 'f')

%% plot some stuff
F = figure(200);
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
set(F, 'Position', [200, -50, 900 700])

% y = [ mean(low_gamma.control) mean(high_gamma.control); mean(low_gamma.ipsi) mean(high_gamma.ipsi); mean(low_gamma.contra) mean(high_gamma.contra)];
% E = [low_control_SEM, high_control_SEM ; low_ipsi_SEM, high_ipsi_SEM; low_contra_SEM, high_contra_SEM];
% y = [1 2 ; 3 4 ; 5 6]
y = [ mean(low_gamma.control) mean(high_gamma.control); mean(low_gamma.ipsi) mean(high_gamma.ipsi); mean(low_gamma.contra), mean(high_gamma.contra)];
E = [low_control_SEM, high_control_SEM; low_ipsi_SEM high_ipsi_SEM; low_contra_SEM high_contra_SEM];
colors = linspecer(numel(y));

b = bar(y, .8);
set(b(1), 'facecolor', colors(2,:))
set(b(2), 'facecolor', colors(3,:))
ylim([0 2])
% set(b(1), 'facecolor', colors(2,:))
% set(b(2), 'facecolor', colors(3,:))
hold on
h = errorbar([.86 1.14; 1.86 2.14; 2.86 3.16], y, E);
for i  = 1:2
set(h(i),'color','k', 'LineStyle','none', 'linewidth', 1)
end
d_x = get(h(1), 'xData');
% d_x2 = get(h(2), 'xData');
% set(h(1), 'xData', d_x-0.15)
% set(h(2), 'xData', d_x+0.15)
ylim([0 2.7])
ylabel({'Gamma event count'; '(normalized to control)'}, 'fontsize', 16)
set(gca, 'fontsize', 14, 'Xticklabel', {'Control','Ipsi', 'Contra'},'FontName', 'helvetica','ytick', [0:.5:2.5], 'fontweight', 'bold')
set(gcf, 'color', 'w')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
%%
figure(200)
hold on
l_width = 1;
sig_font = 30; 
% add bars for low gamma
if h_ip_con == 1;
    ip_con = 1.7:0.001:1.85;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*1.85, ip_con(1,:), 'color',colors(2,:) ,'linewidth', l_width);
    plot((ip_con(2,:).*2.85), ip_con(1,:),'color',colors(2,:),'linewidth', l_width)
    plot(1.85:2.85, [1.85 1.85], 'color',colors(2,:), 'linewidth', l_width)
    if p_ip_con <0.005; text(2.25, 1.86, '**', 'FontSize', sig_font); elseif p_ip_con >0.0051 && p_ip_con <0.05; text(2.25, 1.86, '*', 'FontSize', sig_font) ;end
end

if h_ip_ctr == 1;
    ip_con = 1.45:0.001:1.6;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:)*.85, ip_con(1,:), 'color',colors(2,:),'linewidth', l_width);
    plot((ip_con(2,:).*1.85), ip_con(1,:), 'color',colors(2,:)','linewidth', l_width)
    plot(.85:1.85, [1.6 1.6], 'color',colors(2,:),'linewidth', l_width)
    if p_ip_ctr <0.005; text(1.25, 1.61, '**', 'FontSize', sig_font); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(2.25, 1.61, '*', 'FontSize', sig_font) ;end
end

if h_con_ctr == 1;
    ip_con = 1.1:0.001:1.4;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:)*.85, ip_con(1,:), 'color',colors(2,:));
    plot((ip_con(2,:).*2.85), ip_con(1,:), 'color',colors(2,:))
    plot(1.85:2:2.85, [1.4 1.4], 'color',colors(2,:))
    if p_con_ctr <0.005; text(1.85, 1.31, '**', 'FontSize', sig_font); elseif p_con_ctr >0.0051 && p_con_ctr <0.05; text(2.5, 1.31, '*', 'FontSize', sig_font) ;end
    
end

%% add bars for high gamma
if h_h_ip_con == 1; % ispsi vs contra
    ip_con = (1.7:0.001:1.85)+.4;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*2.15, ip_con(1,:), 'color',colors(3,:) ,'linewidth', l_width);
    plot((ip_con(2,:).*3.15), ip_con(1,:),'color',colors(3,:),'linewidth', l_width)
    plot(2.15:3.15, [1.85+.4 1.85+.4 ], 'color',colors(3,:), 'linewidth', l_width)
    if h_p_ip_con >0.0011 && h_p_ip_con <0.005; text(2.55, 1.85+.4, '**', 'FontSize', sig_font);elseif h_p_ip_con <0.001; text(2.55, 1.85+.4, '***', 'FontSize', sig_font); elseif h_p_ip_con >0.0051 && h_p_ip_con <0.05; text(2.55, 1.85+.4, '*', 'FontSize', sig_font) ;end
end

if h_h_ip_ctr == 1; % ipsi vs control
    ip_con = (1.45:0.001:1.6)+.4 ;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*1.15, ip_con(1,:), 'color',colors(3,:),'linewidth', l_width);
    plot((ip_con(2,:).*2.15), ip_con(1,:), 'color',colors(3,:)','linewidth', l_width)
    plot(1.15:2.15, [1.6+.4  1.6+.4 ], 'color',colors(3,:),'linewidth', l_width)
    if h_p_ip_ctr >0.0011 && h_p_ip_ctr <0.005; text(1.5, 1.61+.4, '**', 'FontSize', sig_font); elseif h_p_ip_ctr <0.001; text(1.5, 1.61+.4 , '***', 'FontSize', sig_font); elseif h_p_ip_ctr >0.0051 && h_p_ip_ctr <0.05; text(2.5, 1.61+.4 , '*', 'FontSize', sig_font) ;end
end

if h_h_con_ctr == 1;
    ip_con = 2.4:0.001:2.55;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:)*1.15, ip_con(1,:), 'color',colors(3,:),'linewidth', l_width);
    plot((ip_con(2,:).*3.15), ip_con(1,:), 'color',colors(3,:),'linewidth', l_width)
    plot(1.15:2:3.15, [2.55 2.55], 'color',colors(3,:),'linewidth', l_width)
    if h_p_con_ctr >0.0011 && h_p_con_ctr <0.005; text(2, 2.55, '**', 'FontSize', sig_font); elseif h_p_con_ctr <0.001; text(2, 2.55, '***', 'FontSize', sig_font);elseif h_p_con_ctr >0.0051 && h_p_con_ctr <0.05; text(2, 2.55, '*', 'FontSize',sig_font) ;end
    
end
% tx_x = 3.3;
% text(tx_x, -.155, '** {\itp}<0.005', 'FontSize', 14)
% text(tx_x, -.1, ' *  {\itp}<0.05', 'FontSize', 14)
% text(tx_x-0.15, -.2, 'Wilcoxon rank sum', 'FontSize', 14)

SetFigure([], gcf)
%% save figure
saveas(F, [PARAMS.naris_figures_dir '\Naris_gamma_count.fig'])
print(F, '-depsc','-r300',[PARAMS.naris_figures_dir 'Naris_gamma_count.eps'])
print(F, '-dpng','-r300',[PARAMS.naris_figures_dir 'Naris_gamma_count.png'])

%% collect stats
stats_out.low.count = low_gamma;
stats_out.high.count = high_gamma;

stats_out.low.stats.ip_ctr = l_stats_ip_ctr;
stats_out.low.stats.ip_con = l_stats_ip_con;
stats_out.low.stats.con_ctr = l_stats_con_ctr;

stats_out.high.stats.ip_ctr = h_stats_ip_ctr;
stats_out.high.stats.ip_con = h_stats_ip_con;
stats_out.high.stats.con_ctr = h_stats_con_ctr;
stats_out.cfg =cfg; 
