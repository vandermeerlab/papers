function err = HDerrfun_xval_2D(xtrain,xtest,particle_vals,obs_ahv,obs_sdf,tc,cfg_param)
% function err = HDerrfun_xval_2D(xtrain,xtest,particle_vals,obs_ahv,obs_sdf,tc,cfg_param)
%
% xtrain: timestamps to fit model with
% xtest: timestamps to test model fit on (obtain yfit)

global param;
global param_count;
global param_hist;

lb = [cfg_param.gain_l(1) cfg_param.gain_r(1) cfg_param.drift(1) cfg_param.hd0(1)];
ub = [cfg_param.gain_l(end) cfg_param.gain_r(end) cfg_param.drift(end) cfg_param.hd0(end)];

[~,~,train_idx] = intersect(xtrain,obs_ahv.tvec);
[~,~,test_idx] = intersect(xtest,obs_ahv.tvec);

% do first pass of getting error for all particles
nP = size(particle_vals,1);
for iP = nP:-1:1
   err(iP) = HDerrfun_mask(particle_vals(iP,:),obs_ahv,obs_sdf,tc,train_idx);
end

% keep track of marginals
clear this_errByBin;
for iDim = 1:4
    this_hist.xcenters{iDim} = cfg_param.histBins{iDim};
    d = median(diff(this_hist.xcenters{iDim}));
    this_hist.xedges{iDim} = this_hist.xcenters{iDim}(1)-d/2:d:this_hist.xcenters{iDim}(end)+d/2;
    
    [~,particle_idx] = histc(particle_vals(:,iDim),this_hist.xedges{iDim});
    
    % now compute average error in each bin
    for iB = 1:length(this_hist.xcenters{iDim})
        
        this_idx = find(particle_idx == iB);
        if isempty(this_idx) % no particles for this bin
            this_errByBin{iDim}(iB) = NaN;
        else
            this_errByBin{iDim}(iB) = nanmean(err(this_idx));
        end
            
    end % of bins
    
end

% do optimization pass with winningest particles as starting points
if cfg_param.nKeep < nP
    [~,sort_idx] = sort(err,'ascend');
    if cfg_param.addRandom
        rand_idx = randperm(nP); rand_idx = rand_idx(1:cfg_param.nKeep);
        particle_vals = particle_vals(unique([sort_idx(1:cfg_param.nKeep) rand_idx]),:);
    else
        particle_vals = particle_vals(sort_idx(1:cfg_param.nKeep),:);
    end
end

HDerrfunA = @(x)HDerrfun_mask(cat(2,x(1:4),0),obs_ahv,obs_sdf,tc,train_idx);
particle_valsO = particle_vals;

clear errO;
opts = optimoptions('fmincon');
opts.Display = 'off';
nP = size(particle_vals,1);
for iP = nP:-1:1
    
    [x,fval] = fmincon(HDerrfunA,particle_vals(iP,1:4),[],[],[],[],lb,ub,[],opts);
    errO(iP) = fval;
    particle_valsO(iP,1:4) = x;
    
end
[~,win_idx] = min(errO);
win_particle = particle_valsO(win_idx,:);

err = HDerrfun_mask(win_particle,obs_ahv,obs_sdf,tc,test_idx);
fprintf('Win gl %.4f gr %.4f d %.3f s %.2f (%.2f)\n',win_particle(1),win_particle(2),win_particle(3),win_particle(4),err);

param(param_count,:) = win_particle(1:4);
param_hist{param_count} = this_errByBin;
param_count = param_count + 1;

function err = HDerrfun_mask(params,obs_ahv,obs_sdf,tc,idx)
%

internal_hd = AHVtoHDfast(obs_ahv,params);
internal_hd = tsd(obs_ahv.tvec,internal_hd);

cfg = []; cfg.mode = 'interp';
predicted_sdf = GenerateSDFfromTC(cfg,internal_hd,tc);

err = (predicted_sdf.data(:,idx)-obs_sdf.data(:,idx)).^2;
%err = abs(predicted_sdf.data(:,idx)-obs_sdf.data(:,idx));
err = nanmean(err(:));