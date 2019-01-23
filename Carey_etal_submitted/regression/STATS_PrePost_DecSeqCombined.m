% compare motivational bias scores for pre and post
%
%% assemble data (sequenceless z)
% assumes PLOT_DecSeqCombinedShuf.m is run
true_pre_z = data.pre.all.median_z;
true_post_z = data.post.all.median_z;
true_type = data.pre.all.this_type;

motiv_bias = @(z,type) nanmean(z(type == 0)) - nanmean(z(type == 1)); % check if right

true_diff = motiv_bias(true_pre_z, true_type) - motiv_bias(true_post_z, true_type);

nShuf = 1000;
for iS = nShuf:-1:1
   
    swap = logical(randi(2, size(true_pre_z)) - 1);
    
    this_pre_z = true_pre_z;
    this_pre_z(swap) = true_post_z(swap);
    
    this_post_z = true_post_z;
    this_post_z(swap) = true_pre_z(swap);
    
    shuf_diff(iS) = motiv_bias(this_pre_z, true_type) - motiv_bias(this_post_z, true_type);
    
end

this_z = (true_diff - nanmean(shuf_diff)) ./ nanstd(shuf_diff);
fprintf('True pre/post diff %.2f, shuf %.2f +/- %.2f, p %.2e\n', true_diff, nanmean(shuf_diff), nanstd(shuf_diff), 2*normcdf(abs(this_z), 0, 1, 'upper'));

%% assemble data (frac)
% assumes PLOT_DecSeqCombinedShuf.m is run
true_pre_z = data.pre.all.fracL_evt;
true_post_z = data.post.all.fracL_evt;
true_type = data.pre.all.this_type;

motiv_bias = @(z,type) nanmean(z(type == 0)) - nanmean(z(type == 1)); % check if right

true_diff = motiv_bias(true_pre_z, true_type) - motiv_bias(true_post_z, true_type);

nShuf = 1000;
for iS = nShuf:-1:1
   
    swap = logical(randi(2, size(true_pre_z)) - 1);
    
    this_pre_z = true_pre_z;
    this_pre_z(swap) = true_post_z(swap);
    
    this_post_z = true_post_z;
    this_post_z(swap) = true_pre_z(swap);
    
    shuf_diff(iS) = motiv_bias(this_pre_z, true_type) - motiv_bias(this_post_z, true_type);
    
end

this_z = (true_diff - nanmean(shuf_diff)) ./ nanstd(shuf_diff);
fprintf('True pre/post diff %.2f, shuf %.2f +/- %.2f, p %.2e\n', true_diff, nanmean(shuf_diff), nanstd(shuf_diff), 2*normcdf(abs(this_z), 0, 1, 'upper'));

%% assemble data (sequences)
% needs data output from PLOT_MotivationalBias_SequencesShufEpochs(), make sure to specify cfg.epochs = {'pre','task','post')
true_pre_z = cat(1, data.pre.all.food_leftN, data.pre.all.water_leftN);
true_post_z = cat(1, data.post.all.food_leftN, data.post.all.water_leftN);
true_type = cat(1, zeros(size(data.post.all.food_leftN)), ones(size(data.post.all.water_leftN)));

motiv_bias = @(z,type) nanmean(z(type == 0)) - nanmean(z(type == 1)); % check if right

true_diff = motiv_bias(true_pre_z, true_type) - motiv_bias(true_post_z, true_type);

nShuf = 1000;
for iS = nShuf:-1:1
   
    swap = logical(randi(2, size(true_pre_z)) - 1);
    
    this_pre_z = true_pre_z;
    this_pre_z(swap) = true_post_z(swap);
    
    this_post_z = true_post_z;
    this_post_z(swap) = true_pre_z(swap);
    
    shuf_diff(iS) = motiv_bias(this_pre_z, true_type) - motiv_bias(this_post_z, true_type);
    
end

this_z = (true_diff - nanmean(shuf_diff)) ./ nanstd(shuf_diff);
fprintf('True pre/post diff %.2f, shuf %.2f +/- %.2f, p %.2e\n', true_diff, nanmean(shuf_diff), nanstd(shuf_diff), 2*normcdf(abs(this_z), 0, 1, 'upper'));