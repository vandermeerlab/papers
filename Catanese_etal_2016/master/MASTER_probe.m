%%
fd = 'D:\data\Minirator-2016-02-04-gamma';

%% define sitemap -- everything is in Neuronexus world
% http://www.neuronexus.com/images/Electrode%20Site%20Map/H32_Maps_20150121.pdf
% neuronexus probe map, looking into connector, Neuronexus side up (as
% shown)
probe_map = [18 27 28 29 17 30 31 32  1  2  3 16  4  5  6 15 ...
             20 21 22 23 19 24 25 26  7  8  9 14 10 11 12 13];

% http://neuralynx.com/manuals/HS-36_Manual.pdf
% neuralynx map, lining up Omnetics connector print with Neuronexus side         
nlx_map =   [17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ...
             1   2  3  4  5  6  7  8  9 10 11 12 13 14 15 16];

% idea is to convert NanoZ impedances and Neuralynx channels into the
% common language of NeuroNexus site IDs. so based on the above, Neuralynx
% channel 1 is probe site 20.
         
% list of defective sites on probe
defective_sites = [5 10 16 19 21 28];

%% load data
cd(fd);
csc = LoadCSC([]); % load all CSCs

for iL = 1:length(csc.label)
    nlx_ch = str2num(csc.label{iL}(4:5))-32; % convert labels to Neuralynx channel labels
    csc.label{iL} = probe_map(nlx_map(nlx_ch)); % now refers to pins on NeuroNexus connector
end
all_ch = cell2mat(csc.label);

%% load event markers
cfg_evt = [];
cfg_evt.fc = 'Events.nev';
evt = LoadEvents(cfg_evt);
evt.label

%% load impedances -- idx of this matches all_ch, so to find impedance of csc.data(iCh), do imp.mag(all_ch(iCh))
fn = 'Probe_4x4_JC_test_feb_2016_test_2.txt';
cfg_z = []; %cfg_z.sitemap = nlx_map;
imp = ReadNanoZ(cfg_z,fn);

%% frequencies and amplitudes for loaded events (manually coded for now based on evt.label)
fvec = [100 40  40  63  63  80  80];
avec = [100 200 100 200 100 100 200];

%%
f = [40 63 80 100];
% f = [40 63 80];
a = 100;
Fs = 2000;

phi = nan(length(f),length(all_ch),length(all_ch)); % phase lag
Zmag1 = nan(length(f),length(all_ch),length(all_ch));
Zmag2 = nan(length(f),length(all_ch),length(all_ch));
Zang1 = nan(length(f),length(all_ch),length(all_ch));
Zang2 = nan(length(f),length(all_ch),length(all_ch));

for iF = 1:length(f)
    
    % restrict to data of interest (recording "trials" are 3 min)
    this_f = f(iF);
    this_evt_idx = find(fvec == this_f & avec == a);
    
    this_data = restrict(csc,evt.t{this_evt_idx},evt.t{this_evt_idx}+180);
    
    % compute phase lags
    for iRef = 1:length(all_ch)
        
        %ref = this_data.data(find(all_ch == reference_site),:);
        ref = this_data.data(iRef,:);
        
        if ~isempty(intersect(all_ch(iRef),defective_sites))
            continue;
        end
        
       % for iCh = length(all_ch):-1:1
       for iCh = iRef+1:length(all_ch)
            
           if ~isempty(intersect(all_ch(iCh),defective_sites))
               continue;
           end
           
            this_ch = this_data.data(iCh,:);
            [Pxy,W] = cpsd(ref,this_ch,hanning(2048),[],[],Fs);
            
            Pxy = angle(Pxy)*(360/(2*pi)); % convert to degrees
            
            fprintf('f %d, pair %d-%d...\n',this_f,iRef,iCh);
            
            phi(iF,iRef,iCh) = Pxy(nearest_idx(this_f,W));
            Zmag1(iF,iRef,iCh) = imp.mag(all_ch(iRef));
            Zmag2(iF,iRef,iCh) = imp.mag(all_ch(iCh));
            Zang1(iF,iRef,iCh) = imp.phase(all_ch(iRef));
            Zang2(iF,iRef,iCh) = imp.phase(all_ch(iCh));
            
        end
        
    end
    
end

%% exmaple fig with one site as ref
figure;
iRef = 1;
for iF = 1:length(f)
   
    this_f = nan(32,1);
    this_f(all_ch) = sq(phi(iF,iRef,:));
    
    this_f(defective_sites) = NaN;
    
    this_f = reshape(this_f,[8 4]);

    subplot(2,2,iF);
    imagesc(this_f);
    
    set(gca,'FontSize',12);
    caxis([-1 1]); axis off; 
    h = colorbar; set(h,'FontSize',18);
    
    title(f(iF));
end

%%
phiR = reshape(phi,[4 30*30]);

figure; subplot(221);

plot(f,phiR,'Color',[0.7 0.7 0.7]);
hold on;
phi_mean = nanmean(phiR,2);
phi_std = nanstd(phiR,[],2);

plot(f,phi_mean,'.k','MarkerSize',20);
plot(f,phi_mean,'k','LineWidth',2);
errorbar(f,phi_mean,phi_std,'k','LineWidth',1);

set(gca,'FontSize',18,'XLim',[35 105],'XTick',40:20:100);
box off;
xlabel('Frequency (Hz)'); ylabel('phase lag (deg)');

%% add mean phase lag and phase slope here: diff(y)/diff(x) plotted in intermediate x values (plotyy)
subplot(222);
fd = f(1:3)+0.5*diff(f);
ps = diff(phi_mean)./diff(f');

plot(fd,ps,'k.','MarkerSize',20); hold on;
plot(fd,ps,'k','LineWidth',2);
set(gca,'FontSize',18,'XLim',[35 105],'XTick',40:20:100,'YLim',[-0.1 0.1]);
box off;
xlabel('Frequency (Hz)'); ylabel('phase slope (deg/Hz)');

%% stats - raw
fR = repmat(f',[1 30*30]);
ZR = (reshape(Zmag1,[4 30*30])+reshape(Zmag2,[4 30*30]))./2;

[P,T,STATS,TERMS] = anovan(phiR(:),{fR(:),ZR(:)},'model','full','varnames',{'frequency','impedance'});

subplot(222);
plot(ZR(:),phiR(:),'.');

keep = ~isnan(phiR(:));
[r,p] = corrcoef(phiR(keep),ZR(keep))

%% stats - abs
fR = repmat(f',[1 30*30]);
ZR = (reshape(Zmag1,[4 30*30])+reshape(Zmag2,[4 30*30]))./2;

[P,T,STATS,TERMS] = anovan(abs(phiR(:)),{fR(:),ZR(:)},'model','full','varnames',{'frequency','impedance'});

subplot(222);
plot(ZR(:),abs(phiR(:)),'.');

keep = ~isnan(phiR(:));
[r,p] = corrcoef(abs(phiR(keep)),ZR(keep))
