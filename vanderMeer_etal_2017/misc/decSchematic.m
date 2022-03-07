%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% produce schematic to illustrate different data splits %%%
%%%                                                       %%%
%%%     Figure 3 in van der Meer et al. (submitted)       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

%%
nTrials = 6;
fz = 20; % font size

% construct matrices
m1 = eye(nTrials);
m2 = ones(nTrials,nTrials);

m3 = circshift(m1,[1 0]);
m4 = -m1+1;
% m5 = toeplitz(mod(0:nTrials-1,2)); % odd/even

% black/white colormap
cm = colormap(gray); 
cm = cm(end:-1:1,:);

%%
what = {'m1','m2','m3','m4'};
ttll = {'same-trial tautology','all-all tautology', ...
    'next-lap generalization','leave-one-out generalization'};

for iPlot = 1:length(what)
    
    subplot(2,2,iPlot);
    
    this_mat = eval(what{iPlot});
    
    imagesc(this_mat); colormap(cm)
    set(gca,'LineWidth',1,'XTick',0.5:1:nTrials,'YTick',0.5:1:nTrials,'TickDir','out','XTickLabel',[],'YTickLabel',[],'FontSize',fz,'YDir','normal');
    grid on;
    ylabel('decoding trial'); xlabel('encoding trial');
    title(ttll{iPlot});
    caxis([0 1]);
    
end