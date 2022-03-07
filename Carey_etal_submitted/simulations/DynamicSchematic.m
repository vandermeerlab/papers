%% set up task sequence
clear

% master parameters -- edit these if needed
epoch_len = [3 9 3]; % how long are the [pre task post] epochs? ALL NEED TO BE INTEGERS
epoch_break = 1; % how long is the break between epochs?
night_len = 30; % how long between sessions?
nCycles = 2; % how many food + water cycles to generate
%first_sess = 0; % start with water (0) or food (1) NOT CURRENTLY IMPLEMENTED
t0 = 5; % when first session starts
% end master parameters

% implied parameters
task_len = sum(epoch_len) + 2 * epoch_break;
cycle_len = 2 * (task_len + night_len); % how long until the pattern repeats?

% construct templates for single session
epoch_template = cat(1, ones(epoch_len(1), 1), zeros(epoch_break, 1), ones(epoch_len(2), 1), zeros(epoch_break, 1), ones(epoch_len(3), 1));
pre_template = cat(1, ones(epoch_len(1), 1), zeros(epoch_break, 1), zeros(epoch_len(2), 1), zeros(epoch_break, 1), zeros(epoch_len(3), 1));
task_template = cat(1, zeros(epoch_len(1), 1), zeros(epoch_break, 1), ones(epoch_len(2), 1), zeros(epoch_break, 1), zeros(epoch_len(3), 1));
post_template = cat(1, zeros(epoch_len(1), 1), zeros(epoch_break, 1), zeros(epoch_len(2), 1), zeros(epoch_break, 1), ones(epoch_len(3), 1));
night_template = cat(1, zeros(task_len, 1), ones(night_len, 1));

break_len = 0:length(epoch_len)-1; break_len = break_len * epoch_break;
epoch_centers = cumsum(cat(2, 0, epoch_len(1:2))) + break_len + ceil(epoch_len / 2); % center index of each epoch within template

% assemble templates into full task sequence
w_starts = t0:cycle_len:t0 + (nCycles - 1) * cycle_len;
f_starts = t0 + cycle_len/2:cycle_len:t0 + cycle_len/2 + (nCycles - 1) * cycle_len;
t1 = f_starts(end) + task_len + night_len;

vars = {'w_pre', 'w_task', 'w_post', 'w_night', 'f_pre', 'f_task', 'f_post', 'f_night', 'session'};
for iVar = 1:length(vars)
   task.(vars{iVar}) = nan(t1, 1); 
end

empties = {'pre_ctr', 'task_ctr', 'post_ctr'};
for iE = 1:length(empties)
   task.(empties{iE}) = []; 
end

sessionCount = 1;
for iS = 1:nCycles
    
   % label epochs 
   this_start = w_starts(iS);
   task.w_pre(this_start + find(pre_template)) = 1; task.w_task(this_start + find(task_template)) = 1; 
   task.w_post(this_start + find(post_template)) = 1; task.w_night(this_start + find(night_template)) = 1;
   task.pre_ctr = cat(1, task.pre_ctr, this_start + epoch_centers(1)); task.task_ctr = cat(1, task.task_ctr, this_start + epoch_centers(2)); task.post_ctr = cat(1, task.post_ctr, this_start + epoch_centers(3));
    
   % need to track what session number we are in: use for variable gain later
   task.session(this_start:this_start + task_len + night_len) = sessionCount;
   sessionCount = sessionCount + 1;
   
   this_start = f_starts(iS);
   task.f_pre(this_start + find(pre_template)) = 1; task.f_task(this_start + find(task_template)) = 1; 
   task.f_post(this_start + find(post_template)) = 1; task.f_night(this_start + find(night_template)) = 1;
   task.pre_ctr = cat(1, task.pre_ctr, this_start + epoch_centers(1)); task.task_ctr = cat(1, task.task_ctr, this_start + epoch_centers(2)); task.post_ctr = cat(1, task.post_ctr, this_start + epoch_centers(3));
    
   task.session(this_start:this_start + task_len + night_len) = sessionCount;
   sessionCount = sessionCount + 1;
   
end

% plot something
task.w_all = nansum(cat(2, task.w_pre, task.w_task, task.w_post), 2);
task.f_all = nansum(cat(2, task.f_pre, task.f_task, task.f_post), 2);

figure; 

for iS = 1:4
subplot(7, 1, (2*iS)-1);
plot(-0.7 * task.w_all, 'LineWidth', 2, 'Color', 'b'); hold on;
plot(0.7 * task.f_all, 'LineWidth', 2, 'Color', 'r');
set(gca, 'YLim', [-0.71 0.71]);
axis off
end

%% construct actual data
%night_bias = [0.3 0.6 0.8 0.5]; % how much to update replay content during the night
night_bias = [0.5 0.5 0.5 0.5]; % how much to update replay content during the night
master_gain = 0.05;

session_gain = zeros(size(task.session));
for iS = 1:length(night_bias) % set gain for each session to behavioral bias
    session_gain(find(task.session == iS)) = night_bias(iS);
end

f = master_gain .* 1 .* session_gain .* task.f_night; % actual gain of experience
w = master_gain .* -1 .* session_gain .* task.w_night;

data = nansum(cat(2, f, w), 2);
data(isnan(data)) = 0;
data(1) = 0.3; % could set different starting point
data(2:4) = 0.02;

data = cumsum(data);

%plot(data, 'k--'); hold on;
%plot(task.pre_ctr, data(task.pre_ctr), '.g', 'MarkerSize', 10);
%plot(task.task_ctr, data(task.task_ctr), '.k', 'MarkerSize', 10);
%plot(task.post_ctr, data(task.post_ctr), '.g', 'MarkerSize', 10);

%% construct some hypotheses
% experience
%behav_bias = [0.5 0.3 0.6 0.8]; % behavioral bias for each session -- must match total sessions used above (nCycles * 2)
behav_bias = [0.5 0.5 0.5 0.5]; % behavioral bias for each session -- must match total sessions used above (nCycles * 2)
master_gain = 0.1;

session_gain = zeros(size(task.session));
for iS = 1:length(behav_bias) % set gain for each session to behavioral bias
    session_gain(find(task.session == iS)) = behav_bias(iS);
end

f = master_gain .* 1 .* session_gain .* task.f_task; % actual gain of experience
w = master_gain .* -1 .* session_gain .* task.w_task;

subplot(711);
plot([1 length(w)], [0 0], ':', 'Color', [0.7 0.7 0.7]); hold on;

exp_f = nansum(cat(2, f, w), 2);
exp_f(isnan(exp_f)) = 0;
exp_f(1) = 0.25; % could set different starting point

exp_f = cumsum(exp_f);

plot(exp_f, 'k'); hold on;
plot(task.pre_ctr, exp_f(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, exp_f(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, exp_f(task.post_ctr), '.g', 'MarkerSize', 10);

plot(data, 'k--'); hold on;
plot(task.pre_ctr, data(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, data(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, data(task.post_ctr), '.g', 'MarkerSize', 10);

%% planning for next trial
subplot(713);
plot([1 length(w)], [0 0], ':', 'Color', [0.7 0.7 0.7]); hold on;

dist_to_w = DistToNextOne(task.w_task, 1:length(task.w_task), Inf);
dist_to_f = DistToNextOne(task.f_task, 1:length(task.f_task), Inf);

plan = nan(size(task.w_task));
plan(dist_to_w < dist_to_f) = -0.5;
plan(dist_to_w > dist_to_f) = 0.5;
plan(dist_to_w == dist_to_f) = -0.5;

plot(plan, 'k');
plot(task.pre_ctr, plan(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, plan(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, plan(task.post_ctr), '.g', 'MarkerSize', 10);

plot(data, 'k--'); hold on;
plot(task.pre_ctr, data(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, data(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, data(task.post_ctr), '.g', 'MarkerSize', 10);
%% delayed experience
subplot(715);
plot([1 length(w)], [0 0], ':', 'Color', [0.7 0.7 0.7]); hold on;

behav_bias = [0.5 0.3 0.6 0.8]; % behavioral bias for each session -- must match total sessions used above (nCycles * 2)
%behav_bias = [0.5 0.5 0.5 0.5]; % behavioral bias for each session -- must match total sessions used above (nCycles * 2)
master_gain = 0.1;

session_gain = zeros(size(task.session));
for iS = 1:length(behav_bias) % set gain for each session to behavioral bias
    session_gain(find(task.session == iS)) = behav_bias(iS);
end

f = master_gain .* 1 .* session_gain .* task.f_task; % actual gain of experience
w = master_gain .* -1 .* session_gain .* task.w_task;

exp_f = nansum(cat(2, f, w), 2);
exp_f(isnan(exp_f)) = 0;

delay = sum(epoch_len);
exp_f = cat(1, zeros(delay, 1), exp_f); exp_f = exp_f(1:end-delay);
exp_f(1) = 0.35; % could set different starting point

exp_f = cumsum(exp_f);

plot(exp_f, 'k');
plot(task.pre_ctr, exp_f(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, exp_f(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, exp_f(task.post_ctr), '.g', 'MarkerSize', 10);


%% (opposite) motivational state
subplot(717);
plot([1 length(w)], [0 0], ':', 'Color', [0.7 0.7 0.7]); hold on;

night_bias = [0.3 0.6 0.8 0.5]; % how much to update replay content during the night
master_gain = 0.04;

session_gain = zeros(size(task.session));
for iS = 1:length(night_bias) % set gain for each session to behavioral bias
    session_gain(find(task.session == iS)) = night_bias(iS);
end

f = master_gain .* 1 .* session_gain .* task.f_night; % actual gain of experience
w = master_gain .* -1 .* session_gain .* task.w_night;

motiv = nansum(cat(2, f, w), 2);
motiv(isnan(motiv)) = 0;
motiv(1) = 0.15; % could set different starting point
motiv(2:4) = 0.02;

motiv = cumsum(motiv);

plot(motiv, 'k'); hold on;
plot(task.pre_ctr, motiv(task.pre_ctr), '.g', 'MarkerSize', 10);
plot(task.task_ctr, motiv(task.task_ctr), '.k', 'MarkerSize', 10);
plot(task.post_ctr, motiv(task.post_ctr), '.g', 'MarkerSize', 10);

%%
