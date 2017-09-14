function AMPX_Naris_PCA(cfg, all_data)
% INPUTS:
%   - "power" is a structure from the output of AMPX_pow_distrib
%   - data is the original data structure.  Just used for naming and
%   finding the correct directories.

global PARAMS
%% collect all power distributions from all sessions per rat.
bands = {'lg', 'hg', 'lg_ran', 'hg_ran', 'spindles'};

% convert to
%% get data ready for pca
% Rats = {'R061', 'R054', 'R049', 'R045'};
% for iRat = 1:length(Rats)
%     for ibands = 1:length(bands)
%         for iT = length(PrePCA.(Rats{iRat}).(bands{ibands})):-1:1
% %             NaN_idx = ~isnan(PrePCA.(Rats{iRat}).(bands{ibands}){iT}); % remove NaNs from the power distributions
%             all_trials.(Rats{iRat}).(bands{ibands})(iT,:) = reshape(PrePCA.(Rats{iRat}).(bands{ibands}){iT}, 1, 64);
%         end
%         % z score the power within each event.
%         all_trials.(Rats{iRat}).(bands{ibands}) = all_trials.(Rats{iRat}).(bands{ibands}) ./ repmat(nansum(all_trials.(Rats{iRat}).(bands{ibands}),2),[1 size(all_trials.(Rats{iRat}).(bands{ibands}),2)]);
%         % run the PCA per rat per band.
%         [PCA.(Rats{iRat}).(bands{ibands}).coeff, PCA.(Rats{iRat}).(bands{ibands}).score, PCA.(Rats{iRat}).(bands{ibands}).latent, PCA.(Rats{iRat}).(bands{ibands}).tsq] = princomp(all_trials.(Rats{iRat}).(bands{ibands}));
%     end
% end

%% test something
all_trials_z_out = [];
for iband = 1:5
    temp = all_data.R049_2014_02_07.(bands{iband}).power.power_distrib;
    
    data.(bands{iband}) = temp; % separate the data based on the bands
    co_length.(bands{iband}) = length(temp);
    all_trials_z = [];
    for iT = length(temp):-1:1
        
        non_nan_idx = ~isnan(temp{iT});
        
        all_trials_z(iT,:) = temp{iT}(non_nan_idx);
        
    end
    
    all_trials_z = all_trials_z ./ repmat(nansum(all_trials_z,2),[1 size(all_trials_z,2)]);
    all_trials.(bands{iband}) = all_trials_z;
    all_trials_z_out = [all_trials_z_out; all_trials_z];
end
[coeff, score, latent, tsq] = pca(all_trials_z_out);

co_eff.lg = coeff(1:co_length.lg);
co_eff.hg = coeff(co_length.lg+1:co_length.hg);
co_eff.lg_ran = coeff(co_length.hg+1:co_length.lg_ran);
co_eff.hg.ran = coeff(co_length.lg_ran+1:co_length.hg_ran);



%% sanity check - reconstruct input data from transformed data
for sg = 1:length(all_trials_z_out(1,:)) % # of dimensions
    both.mean(sg) = mean(all_trials_z_out(:,sg)); %avg obs for each dimension
end

% both.Xtransformed and both.Xmean should be same
Xtransformed = score*coeff'; % get mean centered input data from transformed data
for ii = 1:length(all_trials_z_out(1,:)) % for each dim (cell)
    for jj = 1:length(all_trials_z_out(:,1)) % for each obs (ts)
        both.Xmean(jj,ii) = all_trials_z_out(jj,ii) - both.mean(ii); % mean center the input data
    end
end
isequal(Xtransformed, both.Xmean)
Xtransformed(1:5,1:5)
both.Xmean(1:5,1:5)
%PCA.both.score and both.Xscore should be same
both.Xscore = both.Xmean / coeff'; % recreate transformed data using eigenvectors from PCA
both.Xscore(1:5,1:5)
score(1:5,1:5)

% end sanity check

%% transform the data

for iband = 1:5
    
    %get the mean for each input dimension (corresponds to the sites on the
    %probe which form the columns of the PCA input array "all_trials_z_out"
    for iDim = 1:length(all_trials.(bands{iband})(1,:))
        PCA_dim_mean.(bands{iband})(1,iDim) = mean(all_trials.(bands{iband})(1,:));
    end
    
    %mean center the data points using the mean of each input dimension
    for iDim = 1:length(all_trials.(bands{iband})(1,:))
        for iEvt = 1:length(all_trials.(bands{iband})(:,1))
            data_Xmean.(bands{iband})(iEvt,iDim) = all_trials.(bands{iband})(iEvt,iDim) - PCA_dim_mean.(bands{iband})(iDim); % sub mean of each dim for each obs in that dim
        end
    end
    
    Xscore.(bands{iband}) = data_Xmean.(bands{iband})/coeff';
    
end

%% plot PCs
figure(15)
c_ord = linspecer(length(bands));
for iband = 1:5
    hold on
    if iband ==5
        plot3(Xscore.(bands{iband})(:,1),Xscore.(bands{iband})(:,2),Xscore.(bands{iband})(:,3),'o', 'MarkerFaceColor', c_ord(iband,:), 'MarkerEdgeColor', c_ord(iband,:), 'MarkerSize', 2)
    elseif iband == 3 || iband ==4
        plot3(Xscore.(bands{iband})(:,1),Xscore.(bands{iband})(:,2),Xscore.(bands{iband})(:,3),'+', 'MarkerFaceColor', c_ord(iband,:), 'MarkerEdgeColor', [.2 .2 .2], 'MarkerSize', 1)
    else
        plot3(Xscore.(bands{iband})(:,1),Xscore.(bands{iband})(:,2),Xscore.(bands{iband})(:,3),'o', 'MarkerFaceColor', c_ord(iband,:), 'MarkerEdgeColor', c_ord(iband,:), 'MarkerSize', 2)
    end
    view([-15 15])
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3', 'fontsize', 18)
SetFigure([],gcf)
grid on
l = legend({'Low', 'High', 'Rand l', 'Rand h', 'HVS'}, 'location', 'NorthEast');
a=get(l,'children');
set(a([1:3:end]),'markersize',8); % This line changes the legend marker size
saveas(gcf, [PARAMS.figure_dir '\PCA'], 'png')
saveas(gcf, [PARAMS.figure_dir '\PCA'], 'fig')
saveas(gcf, [PARAMS.figure_dir '\PCA'], 'epsc')
%% quick check for gradients
for iT = 1:9^2
    subplot(9,9,iT)
    nan_imagesc_ec(all_data.R049_2014_02_07.hg_ran.power.power_distrib{iT});
end
