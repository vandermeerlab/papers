%% load ALL_out variable

%% session by session analysis

whatS = {'std','laser'};
whatF = {'kalmanwrapped'};

ALL_err = nan(15,2,4);
ALL_param = nan(15,2,3); ALL_paramSD = nan(15,2,3);
ALL_paramerr = nan(15,2,3);

for iFD = [5 7:10 12 14]; % sessions to analyze -- obtained from criterion in last cell
%for iFD = 1:15
    
    for iS = 1:length(whatS)
    
        for iF = 1:length(whatF)
            
            figure(iFD+iF-1);
        
            this_data = ALL_out(iFD).(whatS{iS}).(whatF{iF});
            
            % normalize
            m0_mean_err = nanmean(this_data.err(1,:));
            m4_mean_err = nanmean(this_data.err(4,:));
            this_data.err = this_data.err./repmat(m0_mean_err,size(this_data.err));
            
            this_mean_err = nanmean(this_data.err,2);
            this_std_err = nanstd(this_data.err,[],2);
            
            % keep track of errors
            ALL_err(iFD,iS,:) = this_mean_err;
            
            for iM = 1:4
                this_gl_mean(iM) = nanmean(this_data.param{iM}(:,1)-1);
                this_gl_std(iM) = nanstd(this_data.param{iM}(:,1)-1);
                
                this_gr_mean(iM) = nanmean(this_data.param{iM}(:,2)-1);
                this_gr_std(iM) = nanstd(this_data.param{iM}(:,2)-1);
                
                this_d_mean(iM) = nanmean(this_data.param{iM}(:,3));
                this_d_std(iM) = nanstd(this_data.param{iM}(:,3));
            end
            
            % keep track of joint fit params
            ALL_param(iFD,iS,:) = cat(1,this_gl_mean(4),this_gr_mean(4),this_d_mean(4));
            ALL_paramSD(iFD,iS,:) = cat(1,this_gl_std(4),this_gr_std(4),this_d_std(4));
            ALL_paramerr(iFD,iS,:) = this_data.wraperr{4};
            
            subplot(2,2,iS);
            
            model_x = 1:4;
            
            bh = barwitherr(cat(1,this_gl_std,this_gr_std,this_d_std)',cat(1,this_gl_mean,this_gr_mean,this_d_mean)');
            set(gca,'XTickLabel',{'M0','M1','M2','M3'}); box off;
            legend('GainL','GainR','drift','Location','Southwest'); legend boxoff;
            if iS == 1
                title(this_data.target_fn);
            else
                title(whatF{iF});
            end
            
            subplot(2,2,iS+2);
            
            model_x = 1:4;
        
            bh = barwitherr(this_std_err',this_mean_err');
            set(gca,'XTickLabel',{'M0','M1','M2','M3'}); box off;
            title(sprintf('M0 err %.2f, M4 err %.2f',m0_mean_err,m4_mean_err));
            
        end
           
    end
        
end

%% plot relative errors
figure(100);
for iS = 1:length(whatS)
   
    subplot(2,2,iS);
    plot(sq(nanmean(ALL_err(:,iS,:))),'.k','MarkerSize',20);
    hold on;
    plot(sq(ALL_err(:,iS,:))','.k');
    
    set(gca,'XTick',1:4,'XTickLabel',{'M0','M1','M2','M3'}); box off;
    title(whatS{iS});
    
    %
    subplot(2,2,iS+2);
    h = plot(sq(ALL_err(:,iS,:))');
    
    set(gca,'XTick',1:4,'XTickLabel',{'M0','M1','M2','M3'}); box off;
    title(whatS{iS});
    if iS == 1
        legend(h,'Location','Southwest');
    end
    
end

%% plot recovered parameters
mean_w = 0.5;

subplot(2,2,1);
m = nanmean(sq(ALL_param(:,1,1:2)));
plot([1-mean_w/2 1+mean_w/2],[m(1) m(1)],'k','LineWidth',2);
hold on;
plot([2-mean_w/2 2+mean_w/2],[m(2) m(2)],'k','LineWidth',2);
plot([1 2],sq(ALL_param(:,1,1:2))','o','MarkerSize',5);
m = nanmean(sq(ALL_param(:,2,1:2)));
plot([3-mean_w/2 3+mean_w/2],[m(1) m(1)],'k','LineWidth',2);
hold on;
plot([4-mean_w/2 4+mean_w/2],[m(2) m(2)],'k','LineWidth',2);
plot([3 4],sq(ALL_param(:,2,1:2))','o','MarkerSize',5);
plot([0.5 4.5],[0 0],'k--');
set(gca,'FontSize',24,'XTick',1:4,'XTickLabel',{'CW','CCW','CW','CCW'},'XLim',[0.5 4.5],'YLim',[-0.2 0.2],'TickDir','out'); box off
ylabel('gain \gamma');

subplot(2,2,2);
m = nanmean(sq(ALL_param(:,1,3)));
plot([3-mean_w/2 3+mean_w/2],[m m],'k','LineWidth',2);
hold on;
plot(3,sq(ALL_param(:,1,3))','o','MarkerSize',5);
m = nanmean(sq(ALL_param(:,2,3)));
plot([4-mean_w/2 4+mean_w/2],[m m],'k','LineWidth',2);
hold on;
plot(4,sq(ALL_param(:,2,3))','o','MarkerSize',5);
plot([0.5 4.5],[0 0],'k--');
set(gca,'FontSize',24,'XTick',1:4,'XTickLabel',{'','','std','laser'},'XLim',[0.5 4.5],'YLim',[-1.5 1.5],'TickDir','out','YAxisLocation','right'); box off
ylabel('drift \delta');

%% pass criterion
iS = 1; % std
iP = 3; % drift
drift_z = sq(ALL_param(:,iS,iP))./sq(ALL_paramSD(:,iS,iP));
keep_idx = find(abs(drift_z) < 3) % use this to input into analysis again