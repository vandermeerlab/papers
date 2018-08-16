%%
Behavior_GenData;

%%
load('C:\users\mvdm\Dropbox\projects\Alyssa\2018-04-25\noSeq_data'); % load output from ALL_Generate_DecSeqCombined.m

%%
cd('C:\users\mvdm\Documents\GitHub\papers\Carey_etal_submitted\decoding_noSeq');

if 1
   for iFD = 1:length(fd)
       fd{iFD}(1) = 'C';
   end
end

PLOT_DecSeqCombined;

%% correlation between session n post and session n + 1 pre
post_median_z = data.post.all.median_z;
post_sess = [2 3 4 5 6 4 5 1 2 3 4 5 6 1 2 3 4 5 6];
post_rat = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));
post_rat = post_rat(data.post.all.this_sess);

pre_median_z = nan(size(post_median_z)); % this will hold next session pre data

% find idx into n + 1 session
for iI = 1:length(post_sess)

    this_sess = post_sess(iI);
    this_rat = post_rat(iI);
    
    temp_idx_out = find(post_rat == this_rat & post_sess == this_sess + 1);
    
    if ~isempty(temp_idx_out)
        pre_median_z(iI) = data.pre.all.median_z(temp_idx_out);
    end
    
end

keep = ~isnan(pre_median_z);
[r,p] = corrcoef(post_median_z(keep),pre_median_z(keep))
plot(post_median_z(keep),pre_median_z(keep),'.');

%% correlation between session n behavior and session n + 1 pre
pre_median_z = nan(size(all.pLeftBehav)); % this will hold pre decoded content
pre_sess = [2 3 4 5 6 4 5 1 2 3 4 5 6 1 2 3 4 5 6];
pre_rat = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));
pre_rat = pre_rat(data.post.all.this_sess);

behav = all.pLeftBehav;
choice = all.pLeftChoice;
behav_sess = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6];
behav_rat = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));

% find idx into next session
for iI = 1:length(behav_sess)
    
    this_sess = behav_sess(iI);
    this_rat = behav_rat(iI);
    
    fprintf('Rat %d, session %d behavior...\n',this_rat,this_sess);
    
    temp_idx_out = find(pre_rat == this_rat & pre_sess == this_sess + 1);
        
    if ~isempty(temp_idx_out)
        fprintf('Found idx %d for next day pre\n',temp_idx_out);
        pre_median_z(iI) = data.pre.all.median_z(temp_idx_out);
    end
    
end

keep = ~isnan(pre_median_z);
[r,p] = corrcoef(behav(keep),pre_median_z(keep))
plot(behav(keep),pre_median_z(keep),'.');

%% "bias score" that gets rid of food vs. water induced (anti)correlations
behav_bias = max(cat(1,behav,1-behav)); % invert normalized behavior below 0.5 to above 0.5
content_bias = abs(pre_median_z);

keep = ~isnan(content_bias);
[r,p] = corrcoef(behav_bias(keep),content_bias(keep))
plot(behav_bias(keep),content_bias(keep),'.');
%% could do a kind of multiple regression/model comparison here: is behavior or motiv state the best predictor?
keep = ~isnan(pre_median_z);
[r,p] = corrcoef(all.sessionType(keep),pre_median_z(keep))

% construct table
tbl = table(behav(keep)',choice(keep)',categorical(all.sessionType(keep)),pre_median_z(keep)',categorical(behav_rat(keep))');
tbl.Properties.VariableNames = {'behav','choice','restr','content','rat'};

modelspec1 = 'content ~ 1 + behav';
modelspec2 = 'content ~ 1 + restr';

glm1 = fitglme(tbl,modelspec1);
glm2 = fitglme(tbl,modelspec2);

compare(glm1,glm2)