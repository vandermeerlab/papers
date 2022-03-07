%% replay content scenarios
%
% approach: find all detected SWR candidates, and assign them hypothetical content based on various scenarios,
% such as "match last trial" or "match total experience"

%% preamble
clear all; pack

%MASTER_root = 'D:\My_Documents\Dropbox\projects\Alyssa'; % replace this with home folder of your project
MASTER_root = 'C:\Users\mvdm\Dropbox\projects\Alyssa'; % replace this with home folder of your project
cd(MASTER_root);
MASTER_path; % reset and then set up path for this project

%% master config
cfg_master = [];
cfg_master.dt = 0.025; % discretization for cumulative sum of trials over time -- useful to find trial for given event
cfg_master.sessions = {'food','water'};
cfg_master.fd_skip = [1 7 8 9 12]; % sessions to skip because not enough cells; cf. van der Meer et al. (2017) Hippocampus
cfg_master.useAllExperience = 0; % if 1, include pre-recording behavior (defined below); if 0, use recording behavior only
cfg_master.preRecordExperience(:,1) = cat(1,repmat(250,[6 1]),repmat(71,[6 1]),repmat(79,[6 1]), repmat(77,[6 1])); % left trials pre-recording
cfg_master.preRecordExperience(:,2) = cat(1,repmat(96,[6 1]),repmat(55,[6 1]),repmat(71,[6 1]), repmat(66,[6 1])); % right trials pre-recording
cfg_master.cand_frac = 0.2; % **TBD** what fraction of candidates to assign content (because not all candidates are L or R sequences in the data)

%% find data folders
fd_cfg = []; fd_cfg.requireCandidates = 1; % don't need previously saved sharp wave-ripple (SWR) candidates here
fd_cfg.requireEvents = 1;
fd = getTmazeDataPath(fd_cfg);

%% initialize variables to be populated
out.left_total = 0; % running totals for all trials for each rat
out.right_total = 0;

what = {'this_frac','all_frac','this_prev','this_fracC','this_prevC'};
epoch = {'pre','task','post'};
for iFD = 1:length(fd)
    for iW = 1:length(what)
        for iE = 1:length(epoch)
            out.(what{iW}).(epoch{iE}).L(iFD) = NaN; 
            out.(what{iW}).(epoch{iE}).R(iFD) = NaN; 
        end
    end
end

%% set up RNG to reproducible state
rng(0);

%% generate data
for iFD = 1:length(fd)
    
    cd(fd{iFD});
    fprintf('Entering session %d/%d...\n',iFD,length(fd));
    
    %% load data
    LoadExpKeys;
    LoadMetadata;
    LoadCandidates; % this loads evt variable with detected SWRs
    
    %% reset running totals if day 1
    if ExpKeys.day == 1
        if cfg_master.useAllExperience % start from known behavioral trial totals (pre-recording)
            out.left_total = cfg_master.preRecordExperience(iFD,1);
            out.right_total = cfg_master.preRecordExperience(iFD,2);
        else % start from zero
            out.left_total = 0; out.right_total = 0;
        end
        
        % also, day 1 means no "post" data is available from previous day,
        % so set that to be empty
        prev_post_frac.L = NaN; prev_post_frac.R = NaN;
        prev_post_prev.L = NaN; prev_post_prev.R = NaN;
        
    end
    
    %% collect this session's data
    left_session = zeros(size(evt.tvec)); % cumsum "experience" left trials for this session
    right_session = zeros(size(evt.tvec)); % cumsum "experience" right trials for this session
    
    left_idx = nearest_idx3(metadata.taskvars.trial_iv_L.tend,evt.tvec);
    left_session(left_idx) = 1;
    left_session = cumsum(left_session);
    
    right_idx = nearest_idx3(metadata.taskvars.trial_iv_R.tend,evt.tvec);
    right_session(right_idx) = 1;
    right_session = cumsum(right_session);
    
    left_total = left_session + out.left_total; % cumsum "experience" left across all sessions so far
    right_total = right_session + out.right_total;
    
    this_session.left_t = metadata.taskvars.trial_iv_L.tend; % trial times for this session
    this_session.right_t = metadata.taskvars.trial_iv_R.tend;
    
    out.session_type(iFD) = find(strcmp(ExpKeys.RestrictionType,cfg_master.sessions));
    
    %% update trial totals so far
    out.left_total = out.left_total + length(metadata.taskvars.trial_iv_L.tstart);
    out.right_total = out.right_total + length(metadata.taskvars.trial_iv_R.tstart);
    
    %% check if sequence content analysis needs to be skipped because not enough cells
    if any(cfg_master.fd_skip == iFD)
        % if we're skipping this session, no post data will be available to
        % carry over
        prev_post_frac.L = NaN; prev_post_frac.R = NaN;
        prev_post_prev.L = NaN; prev_post_prev.R = NaN;
        
        fprintf('Session skipped as per cfg.\n');
        continue;
    end
    
    %% get idxs mapping candidate times to pre, task, post
    pre_idx = evt.tstart < ExpKeys.TimeOnTrack;
    task_idx = evt.tstart >= ExpKeys.TimeOnTrack & evt.tstart < ExpKeys.TimeOffTrack;
    post_idx = evt.tstart >= ExpKeys.TimeOffTrack;
    
    %% find fraction L / (L + R) for all candidate times
    cand_idx = nearest_idx3(evt.tstart,evt.tvec);
    this_frac = left_session(cand_idx) ./ (left_session(cand_idx) + right_session(cand_idx));
    all_frac = left_total(cand_idx) ./ (left_total(cand_idx) + right_total(cand_idx));
    
    %% find previous trial for all candidate times (currently inefficient); 1 is left, 0 is right
    cfg_ff = []; cfg_ff.mode = 'prev';
    for iEvt = length(evt.tstart):-1:1
        [~,fieldname] = FindFieldTime(cfg_ff,this_session,evt.tstart(iEvt));
        if isempty(fieldname) % no preceding event
            this_prev(iEvt) = NaN;
        else
            switch fieldname{1}
                case 'left_t'
                    this_prev(iEvt) = 1;
                case 'right_t'
                    this_prev(iEvt) = 0;
            end
        end
    end
    
    %% convert fractions into sequence counts for each scenario and epoch
    what = {'this_frac','all_frac','this_prev'};
    epoch = {'pre','task','post'};
    for iW = 1:length(what)
          
        for iE = 1:length(epoch)
           
            temp_frac = eval(what{iW});
      
            this_idx = eval(cat(2,epoch{iE},'_idx')); % idxs into events for this epoch
            temp_frac = temp_frac(this_idx);
            
            % now sample randomly from all candidates
            sampleSize = round(length(temp_frac).*cfg_master.cand_frac);
            temp_frac = datasample(temp_frac,sampleSize,'replace',false);
            
            r = rand(size(temp_frac));
            this_R = sum(temp_frac < r); this_L = sum(temp_frac >= r); % assign L/R to event probabilistically
            
            out.(what{iW}).(epoch{iE}).L(iFD) = this_L;
            out.(what{iW}).(epoch{iE}).R(iFD) = this_R;
                    
        end
              
    end
    
    %% create fracC and prevC scenarios that copy over post epoch from previous session into pre epoch
    
    % see if for this session, there is a previous post epoch to copy from
    % needs to set previous_postdata.L and .R to zero when entering day 1, or when a day is skipped
    out.this_fracC.pre.L(iFD) = prev_post_frac.L; out.this_fracC.pre.R(iFD) = prev_post_frac.R;
    out.this_prevC.pre.L(iFD) = prev_post_prev.L; out.this_prevC.pre.R(iFD) = prev_post_prev.R;
    
    % also copy the other fields (task and post)
    out.this_fracC.task.L(iFD) = out.this_frac.task.L(iFD); out.this_fracC.task.R(iFD) = out.this_frac.task.R(iFD);
    out.this_fracC.post.L(iFD) = out.this_frac.post.L(iFD); out.this_fracC.post.R(iFD) = out.this_frac.post.R(iFD);
    
    out.this_prevC.task.L(iFD) = out.this_prev.task.L(iFD); out.this_prevC.task.R(iFD) = out.this_prev.task.R(iFD);
    out.this_prevC.post.L(iFD) = out.this_prev.post.L(iFD); out.this_prevC.post.R(iFD) = out.this_prev.post.R(iFD);
    
    %% update prev_post.L and .R to carry over to next session
    prev_post_frac.L = out.this_frac.post.L(iFD); prev_post_frac.R = out.this_frac.post.R(iFD);
    prev_post_prev.L = out.this_prev.post.L(iFD); prev_post_prev.R = out.this_prev.post.R(iFD);
    
end % of sessions

%% condition data into plottable format
what = {'this_frac','all_frac','this_prev','this_fracC','this_prevC'};
epoch = {'pre','task','post'};
for iW = 1:length(what)
    
    for iE = 1:length(epoch)
        
        food_left = out.(what{iW}).(epoch{iE}).L(out.session_type == 1);
        food_right = out.(what{iW}).(epoch{iE}).R(out.session_type == 1);
        
        water_left = out.(what{iW}).(epoch{iE}).L(out.session_type == 2);
        water_right = out.(what{iW}).(epoch{iE}).R(out.session_type == 2);
        
        food_leftN = food_left ./ (food_left + food_right);
        food_rightN = food_right ./ (food_left + food_right);
        
        water_leftN = water_left ./ (water_left + water_right);
        water_rightN = water_right ./ (water_left + water_right);
        
        out.(what{iW}).(epoch{iE}).food_left = nansum(food_left); out.(what{iW}).(epoch{iE}).food_leftN = nanmean(food_leftN);
        out.(what{iW}).(epoch{iE}).food_right = nansum(food_right); out.(what{iW}).(epoch{iE}).food_rightN = nanmean(food_rightN);
        out.(what{iW}).(epoch{iE}).water_left = nansum(water_left); out.(what{iW}).(epoch{iE}).water_leftN = nanmean(water_leftN);
        out.(what{iW}).(epoch{iE}).water_right = nansum(water_right); out.(what{iW}).(epoch{iE}).water_rightN = nanmean(water_rightN);
    end
end

%% plot data
cfg.colormode = 'inventory3';
FontSize = 8;
%cfg.input_fd = 'D:\projects\AlyssaTmaze\resultsFiles\';
%cfg.output_fd = 'D:\projects\AlyssaTmaze\resultsFiles\viz';
cfg.input_fd = 'C:\temp\';
cfg.output_fd = 'C:\temp\';
cfg.showAllRatsText = 1; % do you want to show the text "all rats" on the combined data figures?
cfg.writeOutput = 1;
cfg.outbasefn = 'Scenario'; % base filename for figure output
colors = TmazeColors(cfg.colormode);
originalFolder = pwd;

%% F1: raw sequence counts
for iW = 1:length(what)
    
    cfg.ylim = [600 300]; cfg.ylimtick = [150 75]; % ALL - overall and single lims
    %cfg.ylim = [3000 1500]; cfg.ylimtick = [750 375]; % ALL - overall and single lims
    
    cd(cfg.input_fd)
    cfg.output_fn = cat(2,cfg.outbasefn,'_',what{iW},'_counts');
    
    pre = out.(what{iW}).(epoch{1});
    task = out.(what{iW}).(epoch{2});
    post = out.(what{iW}).(epoch{3});
    
    ylab = {'Number of significant'; 'sequences'};
    ylimsall = [0 cfg.ylim(1)];
    yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
    ylimssing = [0 cfg.ylim(2)];
    ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);
    
    location = [1 2 4 5]; % where to place the bar
    xlims = [0 location(4)+1];
    
    %% get the accumulated data subplotted first
    rats = {'all'};
    iRat = 1;
    
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    figure; hold on
    
    subplot(4,6,[1 2 7 8]) % for all rats PRE
    d = [pre.food_left pre.food_right pre.water_left pre.water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall)
    xlabel('  food                  water','FontSize',FontSize)
    ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
    title(cat(2,what{iW},' PRERECORD'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
    subplot(4,6,[3 4 9 10]) % for all rats TASKREST
    d = [task.food_left task.food_right task.water_left task.water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',[])
    xlabel('  food                  water','FontSize',FontSize)
    title(cat(2,what{iW},' TASK'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
    subplot(4,6,[5 6 11 12]) % for all rats POST
    d = [post.food_left post.food_right post.water_left post.water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',[])
    xlabel('  food                  water','FontSize',FontSize)
    title(cat(2,what{iW},' POSTRECORD'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
end

%% F2: proportional sequence counts
for iW = 1:length(what)
    
    cfg.ylim = [1 1]; cfg.ylimtick = [0.25 0.25]; % overall and single lims

    cd(cfg.input_fd)
    cfg.output_fn = cat(2,cfg.outbasefn,'_',what{iW},'_props');
    
    pre = out.(what{iW}).(epoch{1});
    task = out.(what{iW}).(epoch{2});
    post = out.(what{iW}).(epoch{3});
    
    ylab = {'Number of significant'; 'sequences'};
    ylimsall = [0 cfg.ylim(1)];
    yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
    ylimssing = [0 cfg.ylim(2)];
    ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);
    
    location = [1 2 4 5]; % where to place the bar
    xlims = [0 location(4)+1];
    
    %% get the accumulated data subplotted first
    rats = {'all'};
    iRat = 1;
    
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    figure; hold on
    
    subplot(4,6,[1 2 7 8]) % for all rats PRE
    d = [pre.food_leftN pre.food_rightN pre.water_leftN pre.water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall)
    xlabel('  food                  water','FontSize',FontSize)
    ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
    title(cat(2,what{iW},' PRERECORD'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
    subplot(4,6,[3 4 9 10]) % for all rats TASKREST
    d = [task.food_leftN task.food_rightN task.water_leftN task.water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',[])
    xlabel('  food                  water','FontSize',FontSize)
    title(cat(2,what{iW},' TASK'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
    subplot(4,6,[5 6 11 12]) % for all rats POST
    d = [post.food_leftN post.food_rightN post.water_leftN post.water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',[])
    xlabel('  food                  water','FontSize',FontSize)
    title(cat(2,what{iW},' POSTRECORD'))
    box off
    set(gca,'Layer','top')
    if cfg.showAllRatsText
        txt = 'all rats';
        text(0.58,0.9,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',FontSize)
    end
    
end

%% export
cfg.output_fn = 'scenario_thisFracC';
print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-painters','-dpdf',[cfg.output_fn,'.pdf']);
print(gcf,'-painters','-depsc',[cfg.output_fn,'.eps']);