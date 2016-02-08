
close all, clear all,
dataFolder = 'D:\Julien_VdmLab\Data\0choice'
FigFolder =  'D:\Julien_VdmLab\Data\Figures'

cd(dataFolder) % cd('D:\Julien_VdmLab\Data\incoming');
Folders = dir;

EpochList= {'Task', 'PostTask'}
% EpochList= {'Task'}

Normalization =1; 

for iEpoch = 1 : max(size(EpochList))
    
    allP1 = [];  allP2 = [];  allP3 = []; allP4 = [];
    allC1 = [];  allC2 = [];  allC3 = [];
    Ns = 0;
    
    for nF =1:max(size(Folders))
        if Folders(nF).isdir==1 & Folders(nF).name(1)=='R' & max(size(Folders(nF).name)) == 15; 
            Ns = Ns+1; 
            disp(Folders(nF).name)
            
            cd ([dataFolder '\' Folders(nF).name])
            mkdir('Figures_png');
            
            load Epoch_list.mat
            run(FindFile('*keys.m'));
            SessionID = ExpKeys.SessionID;
            
            vStr1_csc = loadCSC_JC(ExpKeys.goodGamma_vStr{1});
            vStr2_csc = loadCSC_JC(ExpKeys.goodGamma_vStr{end});
            PFC_csc = loadCSC_JC(ExpKeys.goodGamma_PFC{1});
            HC_csc = loadCSC_JC(ExpKeys.goodGamma_HC{1});
            
            
            
            %% restrict to Task
            vStr1 = Restrict (vStr1_csc, Epoch.(EpochList{iEpoch}).Start,  Epoch.(EpochList{iEpoch}).End);
            vStr2 = Restrict (vStr2_csc, Epoch.(EpochList{iEpoch}).Start,  Epoch.(EpochList{iEpoch}).End);
            PFC   = Restrict (PFC_csc,   Epoch.(EpochList{iEpoch}).Start,  Epoch.(EpochList{iEpoch}).End);
            HC    = Restrict (HC_csc,    Epoch.(EpochList{iEpoch}).Start,  Epoch.(EpochList{iEpoch}).End);
            %
            vStr1 = Data(vStr1); vStr2 = Data(vStr2); PFC = Data(PFC); HC= Data(HC);
            % Next we can compute the PSDs for each signal in the familiar manner, as well as the coherence between signal pairs of interest:
            
            Fs = 2000; wsize = 2048;
            
            % compute PSDs and coherences for each signal and each pair respectively
            [P1,F] = pwelch(vStr1,hanning(wsize),wsize/2,2*wsize,Fs);
            [P2,F] = pwelch(vStr2,hanning(wsize),wsize/2,2*wsize,Fs);
            [P3,F] = pwelch(PFC,hanning(wsize),wsize/2,2*wsize,Fs);
            [P4,F] = pwelch(HC,hanning(wsize),wsize/2,2*wsize,Fs);
            
            if Normalization == 1
                for i=1:max(size(P1))
                    P1(i) = P1(i) .* F(i)^2;
                    P2(i) = P2(i) .* F(i)^2;
                    P3(i) = P3(i) .* F(i)^2;
                    P4(i) = P4(i) .* F(i)^2;
                end
            end
                  
            allP1 = [allP1 P1];
            allP2 = [allP2 P2];
            allP3 = [allP3 P3];
            allP4 = [allP4 P4];
            
            [C1,F] = mscohere(vStr1,vStr2,hanning(wsize),wsize/2,2*wsize,Fs);
            [C2,F] = mscohere(vStr1,PFC,hanning(wsize),wsize/2,2*wsize,Fs);
            [C3,F] = mscohere(vStr1,HC,hanning(wsize),wsize/2,2*wsize,Fs);
            
            allC1 = [allC1 C1];
            allC2 = [allC2 C2];
            allC3 = [allC3 C3];
            
            % plot
            f = figure,
            subplot(121)
            h(1) = plot(F,10*log10(P1),'b','LineWidth',2); hold on;
            h(2) = plot(F,10*log10(P2),'c','LineWidth',2); hold on;
            h(3) = plot(F,10*log10(P3),'r','LineWidth',2); hold on;
            h(4) = plot(F,10*log10(P4),'y','LineWidth',2); hold on;
            h(1) = plot(F,10*log10(P1),'b','LineWidth',2); hold on;
            h(2) = plot(F,10*log10(P2),'c','LineWidth',2); hold on;
            h(3) = plot(F,10*log10(P3),'r','LineWidth',2); hold on;
            
            if Normalization == 1; 
                set(gca,'YLim',[40 60],'YTick',40:2:60, 'XLim',[0 140],'XTick',0:20:140,'FontSize',12); grid on;
            else
                set(gca,'YLim',[0 50],'YTick',0:10:40, 'XLim',[0 160],'XTick',0:20:160,'FontSize',12); grid on;
            end
            legend(h,{'vStr1','vStr2','PFC','HC'},'Location','Northeast'); legend boxoff;
            xlabel('Frequency (Hz)'); ylabel('Power (dB)');
            title([Folders(nF).name ' ' EpochList{iEpoch}])
            
            subplot(122)
            h(1) = plot(F,C1, 'b', 'LineWidth',2); hold on;
            h(2) = plot(F,C2, 'r', 'LineWidth',2); hold on;
            h(3) = plot(F,C3, 'y', 'LineWidth',2); hold on;
            h(1) = plot(F,C1, 'b', 'LineWidth',2); hold on;
            h(2) = plot(F,C2, 'r', 'LineWidth',2); hold on;
            set(gca,'YLim',[0 1.5],'YTick',0:0.5:1, 'XLim',[10 110],'XTick',0:20:120,'FontSize',12); grid on;
            legend(h,{'vStr1-vStr2','vStr1-PFC', 'vStr1-HC'},'Location','Northeast'); legend boxoff;
            xlabel('Frequency (Hz)'); ylabel('Coherence');
            
            if Normalization == 1; 
                saveas(f, [dataFolder  '/' Folders(nF).name '/Figures_png/' SessionID '_PSDnorm_COH_vStr_CPF_HC_' EpochList{iEpoch}],'png')
                saveas(f, [dataFolder  '/' Folders(nF).name '/Figures_png/' SessionID '_PSDnorm_COH_vStr_CPF_HC_' EpochList{iEpoch}],'ai')
            else
                saveas(f, [dataFolder  '/' Folders(nF).name '/Figures_png/' SessionID '_PSD_COH_vStr_CPF_HC_' EpochList{iEpoch}],'png')
                saveas(f, [dataFolder  '/' Folders(nF).name '/Figures_png/' SessionID '_PSD_COH_vStr_CPF_HC_' EpochList{iEpoch}],'ai')
            end
        end
    end
 
    %%
    mP1 = mean(allP1'); %stdP1 = std(allP1',2);
    mP2 = mean(allP2'); %stdP2 = std(allP2',2);
    mP3 = mean(allP3'); %stdP3 = std(allP3',2);
    mP4 = mean(allP4'); %stdP4 = std(allP4',2);
        
    mC1 = mean(allC1'); %stdC1 = std(allC1',2);
    mC2 = mean(allC2'); %stdC2 = std(allC2',2);
    mC3 = mean(allC3'); %stdC3 = std(allC3',2);
    
    
    f2 = figure,
    subplot(121)
    h(1) = plot(F,10*log10(mP1),'b','LineWidth',2); hold on;
    h(2) = plot(F,10*log10(mP2),'c','LineWidth',2); hold on;
    h(3) = plot(F,10*log10(mP3),'r','LineWidth',2); hold on;
    h(4) = plot(F,10*log10(mP4),'y','LineWidth',2); hold on;
    h(1) = plot(F,10*log10(mP1),'b','LineWidth',2); hold on;
    h(2) = plot(F,10*log10(mP2),'c','LineWidth',2); hold on;
    h(3) = plot(F,10*log10(mP3),'r','LineWidth',2); hold on;
    
    
    if Normalization == 1;
        set(gca,'YLim',[40 60],'YTick',40:2:60, 'XLim',[3 130],'XTick',0:30:120,'FontSize',12); grid on;
    else
        set(gca,'XLim',[0 160],'XTick',0:20:160,'FontSize',12); grid on;
    end
    
    legend(h,{'vStr1','vStr2','PFC','HC'},'Location','Northeast'); legend boxoff;
    xlabel('Frequency (Hz)'); ylabel('mean Power (dB)');
    title(['mean over ' num2str(Ns) ' Sessions '  EpochList{iEpoch} '  ' ExpKeys.Task])
    
    subplot(122)
    h(1) = plot(F,mC1, 'b', 'LineWidth',2); hold on;
    h(2) = plot(F,mC2, 'r', 'LineWidth',2); hold on;
    h(3) = plot(F,mC3, 'y', 'LineWidth',2); hold on;
    h(1) = plot(F,mC1, 'b', 'LineWidth',2); hold on;
    h(2) = plot(F,mC2, 'r', 'LineWidth',2); hold on;
    
    set(gca,'YLim',[0 1.5],'YTick',0:0.5:1,'XLim',[0 120],'XTick',0:20:120,'FontSize',12); grid on;
    legend(h,{'vStr1-vStr2','vStr1-PFC', 'vStr1-HC'},'Location','Northeast'); legend boxoff;
    xlabel('Frequency (Hz)'); ylabel('mean Coherence');
    
    if Normalization == 1; 
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSDnorm_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.ai'],'ai')
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSDnorm_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.png'],'png')
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSDnorm_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.fig'],'fig')
    else    
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSD_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.ai'],'ai')
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSD_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.png'],'png')
        saveas(f2, [ 'D:\Julien_VdmLab\Data\Figures\Average_PSD_COH_' EpochList{iEpoch} '_' ExpKeys.Task '.fig'],'fig')
    end
end


