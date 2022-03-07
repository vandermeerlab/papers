function AMPX_Naris_pow_dist_corr(all_data_pre)
%% AMPX_Naris_pow_dist_corr: computes the correlation between the distance
% from the ventrolateral pole and the power in the two gamma bands based on
% the averages for each rat
%
%          Inputs:
%           - all_data_pre
%           -
%           -
%          Outputs:
%           -
%           -
%           -
%
% EC - 2016-07-02
global PARAMS;
%% set up the distance matrix
[x,y] = meshgrid(1.4:-.2:0);

dist = sqrt(((x).^2) +((y).^2));
dist_1D = reshape(dist,1,64);

%% get the correlation for each rat
band = {'lg', 'hg', 'lg_ran', 'hg_ran'};
for iBand = 1:4
    R061_pow.(band{iBand}) = nanmean(cat(3,all_data_pre.R061_2014_09_26.(band{iBand}).power.power_distrib_avg, all_data_pre.R061_2014_09_27.(band{iBand}).power.power_distrib_avg,all_data_pre.R061_2014_09_28.(band{iBand}).power.power_distrib_avg),3);
    R054_pow.(band{iBand}) = nanmean(cat(3,all_data_pre.R054_2014_10_10.(band{iBand}).power.power_distrib_avg,all_data_pre.R054_2014_10_13.(band{iBand}).power.power_distrib_avg),3);
    R045_pow.(band{iBand}) = nanmean(cat(3,all_data_pre.R045_2014_04_16.(band{iBand}).power.power_distrib_avg, all_data_pre.R045_2014_04_17.(band{iBand}).power.power_distrib_avg,all_data_pre.R045_2014_04_15.(band{iBand}).power.power_distrib_avg),3);
    R049_pow.(band{iBand}) = nanmean(cat(3,all_data_pre.R049_2014_02_07.(band{iBand}).power.power_distrib_avg, all_data_pre.R049_2014_02_08.(band{iBand}).power.power_distrib_avg,all_data_pre.R049_2014_02_10.(band{iBand}).power.power_distrib_avg),3);
end

%% normalize within each rat to allow for cross probe comparisons



for iBand = 1:4
    R061_pow_norm.(band{iBand}) =   R061_pow.(band{iBand})/max(max(R061_pow.(band{iBand})));
    R054_pow_norm.(band{iBand}) =   R054_pow.(band{iBand})/max(max(R054_pow.(band{iBand})));
    R045_pow_norm.(band{iBand}) =   R045_pow.(band{iBand})/max(max(R045_pow.(band{iBand})));
    R049_pow_norm.(band{iBand}) =   R049_pow.(band{iBand})/max(max(R049_pow.(band{iBand})));
    
    all_pow.(band{iBand}) = mean(cat(3,R061_pow_norm.(band{iBand}), R054_pow_norm.(band{iBand}), R045_pow_norm.(band{iBand}), R049_pow_norm.(band{iBand})), 3);
end

%%
% [rho, p] =corrcoef(dist, R61_pow.lg, 'rows', 'complete')
for iBand = 1:4
    R061_pow_norm.(band{iBand})(R061_pow_norm.(band{iBand})<=0.000001) = NaN;
    R054_pow_norm.(band{iBand})(R054_pow_norm.(band{iBand})<=0.000001) = NaN;
    R045_pow_norm.(band{iBand})(R045_pow_norm.(band{iBand})<=0.000001) = NaN;
    R049_pow_norm.(band{iBand})(R049_pow_norm.(band{iBand})<=0.000001) = NaN;
end

% R061_pow_norm.hg(R061_pow_norm.hg<=0.000001) = NaN;
% R054_pow_norm.hg(R054_pow_norm.hg<=0.000001) = NaN;
% R045_pow_norm.hg(R045_pow_norm.hg<=0.000001) = NaN;
% R049_pow_norm.hg(R049_pow_norm.hg<=0.000001) = NaN;
%% plot an example
h1 = figure(889)
subplot(121)
c_ord = linspecer(4);
hold on
plot(dist_1D, reshape(R061_pow_norm.lg,1,64), '*', 'color', c_ord(1,:))
plot(dist_1D, reshape(R054_pow_norm.lg,1,64), '*', 'color', c_ord(2,:))
plot(dist_1D, reshape(R045_pow_norm.lg,1,64), '*', 'color', c_ord(3,:))
plot(dist_1D, reshape(R049_pow_norm.lg,1,64), '*', 'color', c_ord(4,:))

plot(dist_1D, reshape(R061_pow_norm.lg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R054_pow_norm.lg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R045_pow_norm.lg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R049_pow_norm.lg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
legend('R1', 'R2', 'R3', 'R4','Ran')
xlabel('Distance from ventrolateral most electrode (mm)','FontSize', 16)
ylabel('Normalized power','FontSize', 16)
set(gca, 'fontsize', 16)
subplot(122)
c_ord = linspecer(4);
hold on
plot(dist_1D, reshape(R061_pow_norm.hg,1,64), '*', 'color', c_ord(1,:))
plot(dist_1D, reshape(R054_pow_norm.hg,1,64), '*', 'color', c_ord(2,:))
plot(dist_1D, reshape(R045_pow_norm.hg,1,64), '*', 'color', c_ord(3,:))
plot(dist_1D, reshape(R049_pow_norm.hg,1,64), '*', 'color', c_ord(4,:))
plot(dist_1D, reshape(R061_pow_norm.hg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R054_pow_norm.hg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R045_pow_norm.hg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
plot(dist_1D, reshape(R049_pow_norm.hg_ran,1,64), '*', 'color', [0.7 0.7 0.7])
legend('R1', 'R2', 'R3', 'R4', 'Ran')
xlabel('Distance from ventrolateral most electrode (mm)', 'FontSize', 16)
ylabel('Normalized power','FontSize', 16)
set(gca, 'fontsize', 16)
SetFigure([],gcf)
maximize
%%
% figure(888)
% hold on
% plot(dist_1D, reshape(all_pow.lg,1,64), '*', 'color', c_ord(1,:))
% plot(dist_1D, reshape(all_pow.hg,1,64), '*', 'color', c_ord(2,:))

% [rho, p] = corrcoef(dist, all_pow.lg, 'rows', 'complete')


%% cat the power for each rat with distance

all_dist = [dist; dist; dist; dist];
all_p.lg = [R061_pow_norm.lg; R054_pow_norm.lg; R045_pow_norm.lg; R049_pow_norm.lg]

[rho_all, p_all] = corrcoef(all_dist, all_p.lg, 'rows', 'complete')


%% try it with the session averages cat together
sess_list = fieldnames(all_data_pre)
all_sess.lg = [];
all_sess.hg = [];
all_sess.lg_ran = [];
all_sess.hg_ran = [];
all_sess_avg.lg = [];
all_sess_avg.hg = [];
all_sess_avg.lg_ran = [];
all_sess_avg.hg_ran = [];
all_dist = [];
all_dist_avg = [];
for iSess = 1:length(sess_list)
    for iBand = 1:4
        for ievt  = 1:length(all_data_pre.(sess_list{iSess}).(band{iBand}).power.power_distrib)
            temp_pow = all_data_pre.(sess_list{iSess}).(band{iBand}).power.power_distrib{ievt} / max(max(all_data_pre.(sess_list{iSess}).(band{iBand}).power.power_distrib{ievt}));
            temp_pow_avg = all_data_pre.(sess_list{iSess}).(band{iBand}).power.power_distrib_avg / max(max(all_data_pre.(sess_list{iSess}).(band{iBand}).power.power_distrib_avg));
            all_sess.(band{iBand}) = [all_sess.(band{iBand}); temp_pow];
            all_sess_avg.(band{iBand}) = cat(3,all_sess_avg.(band{iBand}), temp_pow_avg);
        end
    end
    all_dist_avg = [all_dist_avg; dist];
end
all_dist.lg = repmat(dist, length(all_sess.lg)/8,1);
all_dist.hg = repmat(dist, length(all_sess.hg)/8,1);
all_dist.lg_ran = repmat(dist, length(all_sess.lg_ran)/8,1);
all_dist.hg_ran = repmat(dist, length(all_sess.hg_ran)/8,1);


[sess_rho_lg, sess_p_lg] = corrcoef(all_dist.lg, all_sess.lg, 'rows', 'complete')
[sess_rho_hg, sess_p_hg] = corrcoef(all_dist.hg, all_sess.hg, 'rows', 'complete')
[sess_rho_lg_ran, sess_p_lg_ran] = corrcoef(all_dist.lg_ran, all_sess.lg_ran, 'rows', 'complete')
[sess_rho_hg_ran, sess_p_hg_ran] = corrcoef(all_dist.hg_ran, all_sess.hg_ran, 'rows', 'complete')


%% plot all
h1 = figure(889)
subplot(121)
% all events (used for 
all_sess_1D.lg = reshape(all_sess.lg, 1, numel(all_sess.lg));
all_sess_1D.hg = reshape(all_sess.hg, 1, numel(all_sess.hg));
all_sess_1D.lg_ran = reshape(all_sess.lg_ran, 1, numel(all_sess.lg_ran));
all_sess_1D.hg_ran = reshape(all_sess.hg_ran, 1, numel(all_sess.hg_ran));
% averages for wach session (used for plotting so that things don't get
% crazy)
all_sess_1D_avg.lg = nanmean(all_sess_avg.lg,3);
all_sess_1D_avg.hg = nanmean(all_sess_avg.hg,3);
all_sess_1D_avg.lg_ran = nanmean(all_sess_avg.lg_ran,3);
all_sess_1D_avg.hg_ran = nanmean(all_sess_avg.hg_ran,3);
all_dist_1D.lg = reshape(all_dist.lg, 1, numel(all_dist.lg));
all_dist_1D.hg = reshape(all_dist.hg, 1, numel(all_dist.hg));
all_dist_1D.lg_ran = reshape(all_dist.lg_ran, 1, numel(all_dist.lg_ran));
all_dist_1D.hg_ran = reshape(all_dist.hg_ran, 1, numel(all_dist.hg_ran));
all_dist_1D_avg = reshape(all_dist_avg, 1, numel(all_dist_avg));

hold on
% plot(dist,all_sess_1D_avg.lg, '*', 'color', c_ord(1,:))
% plot(all_dist_1D_avg, all_sess_1D.lg_ran, '*', 'color', [0.7 0.7 0.7])


% Fit the range fro 10 to 20 with a line.
ind = 1:length(all_sess_1D.lg);
k = ~isnan(all_sess_1D.lg);
coeffs = polyfit(all_dist_1D.lg(k),all_sess_1D.lg(k),1)
coeffs_ran = polyfit(all_dist_1D.lg_ran(k),all_sess_1D.lg_ran(k),1)
% coeffs = polyfit(all_dist_1D, all_sess_1D.lg, 1)
% Define the range where we want the line
xFitting = 0:19; % Or wherever...
yFitted = polyval(coeffs, xFitting);
yFitted_ran = polyval(coeffs_ran, xFitting);
% Plot the fitted line over the specified range.
hold on; % Don't blow away prior plot
plot(xFitting, yFitted, 'k', 'LineWidth', 2);
plot(xFitting, yFitted_ran, 'color', [0.7 0.7 0.7], 'LineWidth', 2);
% legend('Low gamma (avg)', 'Gamma Fit','Ran fit' );
xlim([0 2])
xlabel('Distance from VL (mm)','FontSize', 16)
ylabel('Normalized power','FontSize', 16)
set(gca, 'fontsize', 16)

subplot(122)
% plot(dist, all_sess_1D_avg.hg, '*', 'color', c_ord(3,:))

% Fit the range fro 10 to 20 with a line.
ind = 1:length(all_sess_1D.hg);
k = ~isnan(all_sess_1D.hg);
coeffs = polyfit(all_dist_1D.hg(k),all_sess_1D.hg(k),1)
coeffs_ran = polyfit(all_dist_1D.hg_ran(k),all_sess_1D.hg_ran(k),1)
% coeffs = polyfit(all_dist_1D, all_sess_1D.lg, 1)
% Define the range where we want the line
xFitting = 0:19; % Or wherever...
yFitted = polyval(coeffs, xFitting);
yFitted_ran = polyval(coeffs_ran, xFitting);
% Plot the fitted line over the specified range.
hold on; % Don't blow away prior plot
plot(xFitting, yFitted, 'k', 'LineWidth', 2);
plot(xFitting, yFitted_ran, 'color', [0.7 0.7 0.7], 'LineWidth', 2);
% legend('High gamma', 'Line Fit');
xlim([0 2])
xlabel('Distance from VL (mm)','FontSize', 16)
ylabel('Normalized power','FontSize', 16)

set(gca, 'fontsize', 16)
SetFigure([],gcf)
Square_subplots()
maximize
% %% same for the random
% % figure(899)
% subplot(1,2,1)
% hold on
% plot(all_dist_1D, all_sess_1D.lg_ran, '*', 'color', [0.7 0.7 0.7])
% 
% 
% % Fit the range fro 10 to 20 with a line.
% ind = 1:length(all_sess_1D.lg_ran);
% k = ~isnan(all_sess_1D.lg_ran);
% coeffs = polyfit(all_dist_1D(k),all_sess_1D.lg_ran(k),1)
% % coeffs = polyfit(all_dist_1D, all_sess_1D.lg, 1)
% % Define the range where we want the line
% xFitting = 0:19; % Or wherever...
% yFitted = polyval(coeffs, xFitting);
% % Plot the fitted line over the specified range.
% hold on; % Don't blow away prior plot
% plot(xFitting, yFitted, 'color', [0.7 0.7 0.7], 'LineWidth', 2);
% legend('Low gamma', 'Line Fit');
% xlim([0 2])
% xlabel('Distance from VL (mm)','FontSize', 16)
% ylabel('Normalized power','FontSize', 16)
% set(gca, 'fontsize', 16)
% 
% subplot(1,2,2)
% hold on
% plot(all_dist_1D, all_sess_1D.hg_ran, '*', 'color', [0.7 0.7 0.7])
% 
% % Fit the range fro 10 to 20 with a line.
% ind = 1:length(all_sess_1D.hg_ran);
% k = ~isnan(all_sess_1D.hg_ran);
% coeffs = polyfit(all_dist_1D(k),all_sess_1D.hg_ran(k),1)
% % coeffs = polyfit(all_dist_1D, all_sess_1D.lg, 1)
% % Define the range where we want the line
% xFitting = 0:19; % Or wherever...
% yFitted = polyval(coeffs, xFitting);
% % Plot the fitted line over the specified range.
% hold on; % Don't blow away prior plot
% plot(xFitting, yFitted, 'color', [0.7 0.7 0.7], 'LineWidth', 2);
% legend('High gamma', 'Line Fit');
% xlim([0 2])
% xlabel('Distance from VL (mm)','FontSize', 16)
% ylabel('Normalized power','FontSize', 16)
% 
% set(gca, 'fontsize', 16)
% SetFigure([],gcf)
% maximize

%% save the image
saveas(h1, cat(2,PARAMS.figure_dir,'\Fig6'), 'epsc');
saveas(h1, cat(2,PARAMS.figure_dir,'\Fig6'), 'fig')