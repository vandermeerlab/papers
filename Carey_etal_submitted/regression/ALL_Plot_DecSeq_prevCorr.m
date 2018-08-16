% ALL_Plot_DeqSeq_prevCorr
%
% compute correlations between variables (behavior, sequence content)
% across sessions

%% load the behavior
Behavior_GenData; % loads 'all' variable
all_behav = all;

%% set up and load sequence decoding data
cfg.colormode = 'inventory3';
FontSize = 8;
cfg.input_fd = 'C:\temp'; %'D:\projects\AlyssaTmaze\resultsFiles';
cfg.output_fd = 'C:\temp\viz'; %'D:\projects\AlyssaTmaze\resultsFiles\viz';
cfg.showAllRatsText = 1; % do you want to show the text "all rats" on the combined data figures?
cfg.writeOutput = 0;
cfg.input_prefix = 'S1_'; % which files to load? assumes filenames are *DecSeq_$task-phase_all_out.mat
cfg.outbasefn = 'S1_'; % base filename for figure output
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;
biasfun = @(d) (d(1)-d(2))-(d(3)-d(4)); % computes bias measure as (food_left-food_right)-(water_left-water_right) sequence content proportions

% load the data
cd(cfg.input_fd)

all = load(cat(2,cfg.input_prefix,'DecSeq_all_all_eligible_out'));
pre = load(cat(2,cfg.input_prefix,'DecSeq_prerecord_all_eligible_out'));
task = load(cat(2,cfg.input_prefix,'DecSeq_taskrest_all_eligible_out'));
post = load(cat(2,cfg.input_prefix,'DecSeq_postrecord_all_eligible_out'));

%% collect some variables of interest (behavioral)

% behavior of current session (1-24)

behav = all_behav.pLeftBehav;
choice = all_behav.pLeftChoice;
correct_t1 = all_behav.sessionChoice(:,1)'; % trial 1 choice (1 = correct)
behav_sess = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6];
behav_rat = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));

behav_prev = nan(size(behav_sess)); behav_next = nan(size(behav_sess));
choice_prev = nan(size(behav_sess)); choice_next = nan(size(behav_sess));

% behavior of previous & next sessions
for iI = 1:length(behav_sess)
    
    this_sess = behav_sess(iI);
    this_rat = behav_rat(iI);
    
    fprintf('Rat %d, session %d behavior...\n',this_rat,this_sess);
    
    temp_idx_prev = find(behav_rat == this_rat & behav_sess == this_sess - 1);
    if ~isempty(temp_idx_prev)
        fprintf('Found idx %d for prev day behavior\n',temp_idx_prev);
        behav_prev(iI) = behav(temp_idx_prev); choice_prev(iI) = choice(temp_idx_prev);
    end
    
    temp_idx_next = find(behav_rat == this_rat & behav_sess == this_sess + 1);
    if ~isempty(temp_idx_next)
        fprintf('Found idx %d for next day behavior\n',temp_idx_next);
        behav_next(iI) = behav(temp_idx_next); choice_next(iI) = choice(temp_idx_next);
    end
    
end
%% collect some variables of interest (neural)
% want to know:
% prev_post: previous session's post
% prev_all: previous sessions's all
% curr_pre: current session's pre
% curr_all: current session's all

% collect the post data
post_data = post.data.all.ALL_sig_seq;
pre_data = pre.data.all.ALL_sig_seq;
all_data = all.data.all.ALL_sig_seq;

sess_label = repmat(1:6,[1 4]);
rat_label = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));

if max(pre_data.sess) == 19 % collector was run using 19 included sessions only (instead of all 24)
    fprintf('*** manual session IDs used.\n');
    temp_sess = [2 2 3 3 4 4 5 5 6 6 10 10 11 11 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24];
    pre_data.sess = temp_sess; post_data.sess = temp_sess; all_data.sess = temp_sess;  
end

sess_idx = sess_label(temp_sess);
rat_idx = rat_label(temp_sess);

% normalized (N) input data is redundant for left and right, so remove
sess_idx = sess_idx(1:2:end); rat_idx = rat_idx(1:2:end);
curr_pre = nan(24,1); curr_pre(temp_sess(1:2:end)) = pre_data.countN(1:2:end);
curr_post = nan(24,1); curr_post(temp_sess(1:2:end)) = post_data.countN(1:2:end);
curr_all = nan(24,1); curr_all(temp_sess(1:2:end)) = all_data.countN(1:2:end);

% get the prev-post data
prev_post = nan(size(curr_pre)); prev_all = nan(size(curr_pre));

for iI = 1:length(curr_pre)

    this_sess = sess_label(iI);
    this_rat = rat_label(iI);
    
    temp_idx = find(rat_label == this_rat & sess_label == this_sess - 1);
    
    if ~isempty(temp_idx)
        prev_post(iI) = curr_post(temp_idx);
        prev_all(iI) = curr_all(temp_idx);
    end
    
end

%% some correlations
% previous behavior predicts current pre content?
behav_prevB = max(cat(1,behav_prev,1-behav_prev)); % invert normalized behavior below 0.5 to above 0.5
curr_preB = max(cat(1,curr_pre',1-curr_pre')); % invert normalized behavior below 0.5 to above 0.5

keep = ~isnan(curr_preB) & ~isnan(behav_prevB);
[r,p] = corrcoef(behav_prevB(keep),curr_preB(keep))
plot(behav_prevB(keep),curr_preB(keep),'.');

% previous content predicts behavior?
keep = ~isnan(prev_post) & ~isnan(choice_t1)';
[r,p] = corrcoef(prev_post(keep),choice_t1(keep))

prev_postB = max(cat(1,prev_post',1-prev_post'));
[r,p] = corrcoef(prev_postB(keep),choice_t1(keep))
%% multiple regression/model comparison: is behavior or motiv state the best predictor?
[r,p] = corrcoef(behav_prev(keep),curr_preB(keep))
[r,p] = corrcoef(all_behav.sessionType(keep),curr_preB(keep))

% construct table
tbl = table(behav_prev(keep)',categorical(all_behav.sessionType(keep)),curr_pre(keep));
tbl.Properties.VariableNames = {'behav','restr','content'};

modelspec1 = 'content ~ 1 + behav';
modelspec2 = 'content ~ 1 + restr';

glm1 = fitglme(tbl,modelspec1);
glm2 = fitglme(tbl,modelspec2);

compare(glm1,glm2)

%% models of session-by-session behavior
tbl = table(behav_prev',behav',choice_prev',choice',categorical(rat_label)', ...
    categorical(all_behav.sessionType),prev_post,prev_all,curr_pre,curr_all);
    
tbl.Properties.VariableNames = {'behav_prev','behav','choice_prev','choice','ratID',...
    'restr','prev_post','prev_all','curr_pre','curr_all'};

model1 = 'choice ~ 1 + restr';
model2 = 'choice ~ 1 + restr + (1 | ratID)';

glm1 = fitglme(tbl,model1);
glm2 = fitglme(tbl,model2); % subject specific intercept doesn't help
compare(glm1,glm2)

model3 = 'choice ~ 1 + curr_all';
glm3 = fitglme(tbl,model3); % adding replay content doesn't help
compare(glm3,glm1)