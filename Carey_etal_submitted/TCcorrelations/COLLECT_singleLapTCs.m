% COLLECT_singleLapTCs
%
% analyze relationship between L and R trial tuning curves (Figure S2c-d)
%
% uses output from GENERATE_singleLapTCs.m (run that first)

%% init
fd = getTmazeDataPath([]);
cfg.prefix = 'R0_'; % specifies which DecSeq data file (from Generate_DecSeq.m) to load -- only used for finding choice point!

nFD = 24; nMaxCells = 150; nMaxLaps = 20;
left_sameC = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps); % C matrices for tracking correlations between tuning curves
right_sameC = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps);
left_differentC = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps);
left_sameD = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps); % D matrices for euclidean distances between tuning curves
right_sameD = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps);
left_differentD = nan(nFD,nMaxCells,nMaxLaps,nMaxLaps);
%%
for iFD = 1:length(fd)
    
    this_tc = ALL_TC{iFD}; % get this from GENERATE_singleLapTCs.m
    
    if isempty(this_tc) % session was skipped in generate process
        continue;
    end
    
    % need to find linpos cp
    cd([fd{iFD},'\files'])
    [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
    this_file = FindFiles([cfg.prefix,sessionID,'-DecSeq_data.mat']);
    
    if isempty(this_file)
        fprintf('Session %s: no DecSeq file found, skipping...\n',fd{iFD});
        continue;
    end
    
    load(this_file{1}); cd .. % loads 'out' variable saved by ALL_Generate_DecSeq.m
    left_cp_bin = out.expCond(1).cp_bin; right_cp_bin = out.expCond(2).cp_bin;
    
    % first TC -- left
    nLapsL = length(this_tc.left);
    for iL = 1:nLapsL 
        
        this_leftL = this_tc.left{iL};

        nCells = size(this_leftL.tc,1);
        
        % second TC -- same
        nLapsK = length(this_tc.left);  
        for iK = iL+1:nLapsK 
            
            this_leftK = this_tc.left{iK};
            
            for iC = 1:nCells
                
                % restrict to pre-CP
                this_left_tcL = this_leftL.tc(iC,1:left_cp_bin);
                this_left_tcK = this_leftK.tc(iC,1:left_cp_bin);
                
                % correlate TC for each cell in this lap with TCs of all other laps
                [left_sameC(iFD,iC,iL,iK),left_sameD(iFD,iC,iL,iK)] = cc_skipnan(this_left_tcL,this_left_tcK);
                
            end % of cells
            
        end % of K-laps
        
        % second TC -- different
        nLapsK = length(this_tc.right);
        for iK = 1:nLapsK
            
            this_rightK = this_tc.right{iK};
            
            for iC = 1:nCells
                
                % find matching cell
                k_idx = strmatch(this_leftL.label{iC},this_rightK.label);
                
                if ~isempty(k_idx)
                    
                    % restrict to pre-CP
                    this_right_tcK = this_rightK.tc(k_idx,1:right_cp_bin);
                    
                    % correlate TC for each cell in this lap with TCs of all other laps
                    [left_differentC(iFD,iC,iL,iK),left_differentD(iFD,iC,iL,iK)] = cc_skipnan(this_left_tcL,this_right_tcK);
                    
                end
                
            end % of cells
            
        end % of K-laps
        
    end % of L-laps
    
    % first TC -- right
    nLapsL = length(this_tc.right);
    for iL = 1:nLapsL 
        
        this_rightL = this_tc.right{iL};

        nCells = size(this_rightL.tc,1);
        
        % second TC -- same
        nLapsK = length(this_tc.right);  
        for iK = iL+1:nLapsK 
            
            this_rightK = this_tc.right{iK};
            
            for iC = 1:nCells
                
                % restrict to pre-CP
                this_right_tcL = this_rightL.tc(iC,1:right_cp_bin);
                this_right_tcK = this_rightK.tc(iC,1:right_cp_bin);
                
                % correlate TC for each cell in this lap with TCs of all other laps
                [right_sameC(iFD,iC,iL,iK),right_sameD(iFD,iC,iL,iK)] = cc_skipnan(this_right_tcL,this_right_tcK);
                
            end % of cells
            
        end % of K-laps
                
    end % of L-laps
    
    
    % some quick plots within session
    Sleft_sameC = sq(nanmean(sq(left_sameC(iFD,:,:,:)))); Sleft_sameD = sq(nanmean(sq(left_sameD(iFD,:,:,:))));
    Sright_sameC = sq(nanmean(sq(right_sameC(iFD,:,:,:)))); Sright_sameD = sq(nanmean(sq(right_sameD(iFD,:,:,:))));
    Sleft_differentC = sq(nanmean(sq(left_differentC(iFD,:,:,:)))); Sleft_differentD = sq(nanmean(sq(left_differentD(iFD,:,:,:))));
    fprintf('Session %d: left-same %.2f, right-same %.2f, different %.2f (corr)\n',iFD,nanmean(Sleft_sameC(:)),nanmean(Sright_sameC(:)),nanmean(Sleft_differentC(:)));
  
    % may need to normalize within session here?
    
end % of sessions

%% need to reshape so that all cells can be averaged
lsc = reshape(left_sameC,[nFD nMaxCells nMaxLaps.^2]);
rsc = reshape(right_sameC,[nFD nMaxCells nMaxLaps.^2]);
ldc = reshape(left_differentC,[nFD nMaxCells nMaxLaps.^2]);

% average over laps 
lsc = sq(nanmean(lsc,3));
rsc = sq(nanmean(rsc,3));
ldc = sq(nanmean(ldc,3));

% by-session average
lscS = sq(nanmean(lsc,2));
rscS = sq(nanmean(rsc,2));
ldcS = sq(nanmean(ldc,2));
sS = nanmean(cat(2,lscS,rscS),2);
nS = sum(~isnan(sS));

[h,p] = ttest2(sS,ldcS);

% make session scatterplot 
subplot(221);
plot(sS,ldcS,'.k','MarkerSize',5); grid on;
set(gca,'LineWidth',1,'TickDir','out','FontSize',18,'XLim',[-0.2 1],'YLim',[-0.2 1]);
hold on;
plot([-0.2 1],[-0.2 1],'k--');
xlabel('same'); ylabel('different');
title(sprintf('n = %d sessions, same %.2f +/- %.2f, different %.2f +/- %.2f, p = %1.2e',nS,nanmean(sS),nanstd(sS),nanmean(ldcS),nanstd(ldcS),p));

% reshape to average over cells
lscC = reshape(lsc,[nFD*nMaxCells 1]);
rscC = reshape(rsc,[nFD*nMaxCells 1]);
ldcC = reshape(ldc,[nFD*nMaxCells 1]);
sC = nanmean(cat(2,lscC,rscC),2);
nC = sum(~isnan(sC));

[h,p] = ttest2(sS,ldcS);

% make cell scatterplot 
subplot(222);
plot(sC,ldcC,'.k','MarkerSize',5); grid on;
set(gca,'LineWidth',1,'TickDir','out','FontSize',18,'XLim',[-0.2 1],'YLim',[-0.2 1]);
hold on;
plot([-0.2 1],[-0.2 1],'k--');
xlabel('same'); ylabel('different');
title(sprintf('n = %d cells, same %.2f +/- %.2f, different %.2f +/- %.2f, p = %1.2e',nC,nanmean(sC),nanstd(sC),nanmean(ldcC),nanstd(ldcC),p));


%%
function [cc,pd] = cc_skipnan(in1,in2)
% return correlation coefficient and euclidean distance bewteen two input
% vectors after removing nans

in1_nan = find(~isnan(in1));
in2_nan = find(~isnan(in2));
keep_idx = intersect(in1_nan,in2_nan);

in1 = in1(keep_idx); in2 = in2(keep_idx);

cc = corrcoef(in1,in2);
cc = cc(1,2);

pd = pdist(cat(1,in1,in2),'euclidean');

end