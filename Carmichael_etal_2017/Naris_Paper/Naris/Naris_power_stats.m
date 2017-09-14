function stats_out = Naris_power_stats(cfg_in, all_naris)


%% Initialize
global PARAMS
cfg_def.marker_ord = {'x', '+', '*'};
cfg_def.c_order = linspecer(4);
cfg_def.c_order_pow = [ 255 66 0; 102 25 255; 64 155 191]/255; % specific colours to match the histology markers in figure 1a
cfg_def.phases = {'pre', 'ipsi', 'contra', 'post'};
% cfg.gamma_freq = [70 85];
cfg_def.gamma_freq = [45 65];
% cfg.gamma_freq = [6 8];
cfg_def.notch = [59.5 61.5 ];
cfg_def.plots =1; % turn on if you want to see plots for each subject across sessions and phases.
cfg_def.normalize = 1; % turn this off if you do not want to use the normalized power for the figure 2 plot.
cfg_def.dev_order = 1; % the order value for the white noise filter gradient function See ft_preproc_derivative.m

cfg_def.whitefilter = 'yes'; 
cfg = ProcessConfig2(cfg_def, cfg_in);

if cfg.gamma_freq(1) <= 45
    cfg.band = 'low';
else
    cfg.band = 'high';
end
%% set up some basic variables
ids = fieldnames(all_naris); % gets the number of rats in the all_naris data

%% get some averages for the gamma bands in each case
for iphase = 1:4
    comp_data.(cfg.phases{iphase}) = [];
end
% if ~strcmp(ids{1}, 'R053')
%     ids = flip(ids); 
% end
% get the averages within the desired frequencies
for id = 1:length(ids) % make sure to go backwards so that R053 is the first one. 
    sess_list = fieldnames(all_naris.(ids{id})); 
    for iphase = 1:length(cfg.phases)
        if cfg.plots ==1
            figure(id)
            subplot(2,2,iphase)
        end
        for nSess = 1:4
            phases = fieldnames(all_naris.(ids{id}).(sess_list{nSess}));

            %find the notch values and replace them with NaNs
            freq = all_naris.(ids{id}).(sess_list{nSess}).(phases{iphase}).data{1}.Frequencies;
            notch_range = find(freq >cfg.notch(1) & freq < cfg.notch(2));
            gamma_range = find(freq >cfg.gamma_freq(1) & freq < cfg.gamma_freq(2));
            % replace the notch frequencies with NaNs and then get the average power within the gamma range in cfg.gamma_freq
            data = all_naris.(ids{id}).(sess_list{nSess}).(phases{iphase}).data{1}.Data;
            data(notch_range) = NaN; % replace notch with NaNs
            comp_data.(cfg.phases{iphase})(nSess, id) = nanmean(data(gamma_range));
            disp([(cfg.phases{iphase}) ' Sess# ' num2str(nSess) ' mean: ' num2str(nanmean(data(gamma_range)))])
            if cfg.plots ==1
                hold on
                plot(freq, 10*log10(data), 'Color', cfg.c_order(nSess,:))
                if strcmp(ids{id}, 'R053')
                    set(gca,'XLim',[0 100],'XTick',0:20:100,'YLim',[-20 10],'XTickLabel',{0 20 40 60 80 100},'YTick',[]); grid on;
                else
                    set(gca,'XLim',[0 100],'XTick',0:20:100,'YLim',[-135 -115],'XTickLabel',{0 20 40 60 80 100},'YTick',[]); grid on;
                end
                xlabel(phases{iphase})
            end
        end
        if cfg.plots ==1; h = legend('1', '2', '3', '4'); v= get(h, 'title'); set(v, 'string', 'Session num'); end
    end
            SetFigure([], gcf)
end
%% convert the comp_data to a format for plotting all the points.

% normalize to the pre recording
if isfield(cfg, 'whitefilter')==0
    comp_data.ipsi = comp_data.ipsi./comp_data.pre;
    comp_data.contra = comp_data.contra./comp_data.pre;
    comp_data.post = comp_data.post./comp_data.pre;
    comp_data.pre = comp_data.pre./comp_data.pre; % this needs to be last or everything gets divided by 1.
else
    
    %% or normalize to the average between the pre and the post.
    % if cfg.normalize == 2
    comp_data.ipsi = comp_data.ipsi./((comp_data.pre+comp_data.post)/2);
    comp_data.contra = comp_data.contra./((comp_data.pre+comp_data.post)/2);
    comp_data.control = ((comp_data.pre+comp_data.post)/2)./((comp_data.pre+comp_data.post)/2); % this needs to be last or everything gets divided by 1.
    % end
end
%% plot everything on to one figure

figHandles = get(0,'Children');
if sum(figHandles == 200) > 0
    close(200)
end
cfg.conditions = {'control', 'ipsi', 'contra'};

F = figure(200);
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
set(F, 'Position', [200, 200, 900 700])
if isfield(cfg, 'whitefilter')
    cfg.phases = cfg.conditions;
    note = ('w/diff');
    labels = {'Control', 'Ipsi', 'Contra'};
else
    labels = {'Pre', 'Ipsi', 'Contra', 'Post'};
    note = '';
end
hold on
for iphase = 1:length(cfg.phases)
    for id= 1:length(ids)
        if id ==1 % "*iphase-(1/10) is to give a slight horizontal offset to the data for easier viewing. 
            h1 = plot(ones(4,3)*iphase-(1/10),comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3) ;
        elseif id==2
            h2 = plot(ones(4,3)*iphase,comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3);
        elseif id ==3
            h3 = plot(ones(4,3)*iphase+(1/10),comp_data.(cfg.phases{iphase})(:,id),  cfg.marker_ord{id}, 'MarkerEdgeColor', cfg.c_order_pow(id,:), 'MarkerSize', 14, 'LineWidth', 3);
        end
        disp(num2str(comp_data.(cfg.phases{iphase})(:,id)))
    end
%         legend([h1(1), h2(1), h3(1)], 'R2', 'R4', 'R5', 'orientation', 'horizontal', 'fontsize', 16)

    legend([h1(1), h2(1), h3(1)],{'R5', 'R6', 'R7'}, 'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 16)
    line_of_best_fit(1,iphase) = mean(nanmean(comp_data.(cfg.phases{iphase}),1));
end

plot(1:length(cfg.phases), line_of_best_fit, 'k', 'LineWidth', 1)
% set the figure properties
box on
set(gcf, 'color', [1 1 1])
ylabh = get(gca,'yLabel');
% ylab_str('interpreter','latex','string','\fontsize{20}{0}\selectfont$Power$\fontsize{16}{0}\selectfont$(Normalized)$');
% ylabel('\fontsize{20}{0}\selectfont$Power$\n\fontsize{16}{0}\selectfont$(Normalized)$','Interpreter','LaTex')
ylabel('Power (normalized)')%, 'position', get(ylabh,'Position') - [.2 0 0] )
text(2, -3, 'Contidtion', 'FontSize', 20)
SetFigure([], gcf)
%% Stats and output
ipsi = reshape(comp_data.ipsi,1,numel(comp_data.ipsi));
control = reshape(comp_data.control,1,numel(comp_data.control));
contra = reshape(comp_data.contra,1,numel(comp_data.contra));

[h_ip_con, p_ip_con,~, x_ip_con] = ttest(ipsi, contra);
[h_ip_ctr, p_ip_ctr, ~, x_ip_ctr] = ttest(ipsi, control);
[h_con_ctr, p_con_ctr,~,x_con_ctr] = ttest(contra, control);

%% stats text out
fprintf('\n\n\nBasic Stats\n')
fileID = fopen([PARAMS.stats_dir '\Naris_stats_power_' cfg.band '.txt'],'w');
fprintf(fileID, ['\n_________________' cfg.band '________________________\n'])
fprintf(fileID,['Ipsilateral:     Mean ' num2str(mean(ipsi), '%4.4f') '   SEM +/-' num2str(std(ipsi)/sqrt(length(ipsi))) '\n'])
fprintf(fileID,['Contralateral:   Mean ' num2str(mean(contra),'%4.4f') '   SEM +/-' num2str(std(contra)/sqrt(length(contra))) '\n'])
fprintf(fileID,['Control:         Mean ' num2str(mean(control),'%4.4f') '   SEM +/-' num2str(std(control)/sqrt(length(control))) '\n'])

fprintf(fileID,'\nt-test\n')
fprintf(fileID,['Ipsilateral   vs. Contralateral:  df: ' num2str(x_ip_con.df) ' t: ' num2str(x_ip_con.tstat)  '   P:' num2str(p_ip_con, '%4.4f') '\n' ])
fprintf(fileID,['Ipsilateral   vs. Control:        df: ' num2str(x_ip_con.df) ' t: ' num2str(x_ip_ctr.tstat)  '   P:' num2str(p_ip_ctr, '%4.4f') '\n' ])
fprintf(fileID,['Contralateral vs. Control:        df: ' num2str(x_ip_con.df) ' t: ' num2str(x_con_ctr.tstat) '   P:' num2str(p_con_ctr, '%4.4f') '\n' ])

fclose(fileID);
% copyfile(['G:\Naris\Naris_stats_power_' cfg.band '.txt'],'D:\Users\mvdmlab\My_Documents\GitHub\Papers_EC\vStr_Naris', 'f')
%% update figure
l_width = 1;
sig_font = 30; 
max_val = max(max([ipsi, contra, control])); % used to determine the significance marker offset
figure(200)
hold on
if h_ip_con == 1;
    ip_con = max_val+0.1:0.001:max_val+.2;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*2, ip_con(1,:), '-k','linewidth', l_width);
    plot((ip_con(2,:).*3), ip_con(1,:), '-k','linewidth', l_width)
    plot(2:3, [ip_con(1,end) ip_con(1,end)], '-k','linewidth', l_width)
    if p_ip_con <0.005; text(2.45, ip_con(1,end)+0.01, '**', 'FontSize', sig_font); elseif p_ip_con >0.0051 && p_ip_con <0.05; text(2.5, ip_con(1,end)+0.01, '*', 'FontSize', sig_font) ;end
end

if h_ip_ctr == 1;
    ip_con = max_val+0.05:0.001:max_val+.15;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:), ip_con(1,:), '-k','linewidth', l_width);
    plot((ip_con(2,:).*2), ip_con(1,:), '-k','linewidth', l_width)
    plot(1:2, [ip_con(1,end) ip_con(1,end)], '-k','linewidth', l_width)
    if p_ip_ctr <0.005; text(1.45, ip_con(1,end)+0.01, '**', 'FontSize', sig_font); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(1.5, ip_con(1,end)+0.01, '*', 'FontSize', sig_font) ;end
end

if h_con_ctr == 1;
    ip_con = max_val+0.2:0.001:max_val+.3;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:), ip_con(1,:), '-k','linewidth', l_width);
    plot((ip_con(2,:).*3), ip_con(1,:), '-k','linewidth', l_width)
    plot(1:2:3, [ip_con(1,end) ip_con(1,end)], '-k','linewidth', l_width)
    if p_con_ctr <0.005; text(1.85, ip_con(1,end)+0.01, '**', 'FontSize', sig_font); elseif p_con_ctr >0.0051 && p_con_ctr <0.05; text(1.95, ip_con(1,end)+0.01, '*', 'FontSize', sig_font) ;end
    
end
% tx_x = 3.1;
% text(tx_x, .4, '** {\itp}<0.005', 'FontSize', 18)
% text(tx_x, .35, ' *  {\itp}<0.05', 'FontSize', 18)
% text(tx_x, .3, 'Wilcoxon Sign Rank Test', 'FontSize', 12)
set(gca,'XLim',[0.5 (length(cfg.phases)+0.5)],'XTick',1:1:length(cfg.phases),'XTickLabel',labels,'Ylim', [.25 max_val+.5],'YTick',[.25:.25:max_val+.5], 'FontSize', 20) %'Color', [1 1 1]) ; grid off;

SetFigure([], gcf)
%%

saveas(F, [PARAMS.naris_figures_dir '\Naris_gamma_power_' cfg.band '.fig'])
print(F, '-depsc','-r300',[PARAMS.naris_figures_dir '\Naris_gamma_power_' cfg.band '.eps'])
print(F, '-dpng','-r300',[PARAMS.naris_figures_dir '\Naris_gamma_power_' cfg.band '.png'])
% print(F, '-deps','-r300',['Naris_comp_high.eps'])

%% collect the stats
stats_out.cfg = cfg;

stats_out.power = comp_data;
stats_out.labels = ids;


stats_out.ip_ctr = x_ip_ctr;
stats_out.ip_con = x_ip_con;
stats_out.con_ctr = x_con_ctr;


