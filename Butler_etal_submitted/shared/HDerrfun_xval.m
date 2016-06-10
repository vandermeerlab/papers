function yfit = HDerrfun_xval(xtrain,ytrain,xtest,particle_vals,obs_ahv,obs_sdf,tc,cfg_param)
% function yfit = HDerrfun_xval(xtrain,ytrain,xtest,particle_vals,obs_ahv,obs_sdf,tc,cfg_param)
%
% xtrain: timestamps to fit model with
% ytrain: nCells x 1 observed SDFs, to fit model with
% xtest: timestamps to test model fit on (obtain yfit)

nKeep = 10; % best particles to keep

[~,~,train_idx] = intersect(xtrain,obs_ahv.tvec);
[~,~,test_idx] = intersect(xtest,obs_ahv.tvec);

% do first pass of getting error for all particles
nP = size(particle_vals,1);
for iP = nP:-1:1
   err(iP) = HDerrfun_mask(particle_vals(iP,:),obs_ahv,obs_sdf,tc,train_idx);
end

[~,sort_idx] = sort(err,'ascend');
particle_vals = particle_vals(sort_idx(1:nKeep),:);

% do optimization pass
HDerrfunA = @(x)HDerrfun_mask(cat(2,x(1:4),0),obs_ahv,obs_sdf,tc,train_idx);
lb = [cfg_param.gain_l(1) cfg_param.gain_r(1) cfg_param.drift(1) cfg_param.hd0(1)];
ub = [cfg_param.gain_l(end) cfg_param.gain_r(end) cfg_param.drift(end) cfg_param.hd0(end)];
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

% generate yfit
internal_hd = AHVtoHDfast(obs_ahv,win_particle);
internal_hd = tsd(obs_ahv.tvec,internal_hd);

cfg = []; cfg.mode = 'interp';
predicted_sdf = GenerateSDFfromTC(cfg,internal_hd,tc);
yfit = predicted_sdf.data(1,test_idx)';

if any(isnan(yfit))
    error;
end


function err = HDerrfun_mask(params,obs_ahv,obs_sdf,tc,train_idx)
%

internal_hd = AHVtoHDfast(obs_ahv,params);
internal_hd = tsd(obs_ahv.tvec,internal_hd);

cfg = []; cfg.mode = 'interp';
predicted_sdf = GenerateSDFfromTC(cfg,internal_hd,tc);

err = (predicted_sdf.data(train_idx)-obs_sdf.data(train_idx)).^2;
err = nanmean(err(:));