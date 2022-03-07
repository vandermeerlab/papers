%% load
cd('D:\data\R016\R016-2012-10-01_promoted'); % has a cell or 2
cd('D:\data\R016\R016-2012-10-03_promoted'); % nice
S = LoadSpikes(cat(1,FindFiles('*.t'),FindFiles('*._t')));

evt = getEvents;

%% get events
all_np = cat(2,evt.n0,evt.n1);
all_cue = cat(2,evt.c1,evt.c5);

cfg = []; cfg.mode = 'prev';

np1 = []; np5 = [];
for iNP = 1:length(all_np)
    [prev,fn] = FindFieldTime(cfg,evt,all_np(iNP));
    switch fn{1}
        case 'c1'
            np1 = cat(2,np1,all_np(iNP));
        case 'c5'
            np5 = cat(2,np5,all_np(iNP));
    end
end
    

%% plot
for iC = 1:length(S)
   
    % np
    
    subplot(211)
    
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},all_np,'window',w,'dt',dt);
    
    tbin = linspace(w(1), w(2), diff(w)/dt+1);
    
    m = histc(a, tbin);
	bar(tbin,m/dt/length(all_np));
    
    hold on;
    
    % 1p
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},np1,'window',w,'dt',dt);
     
    m = histc(a, tbin);
    plot(tbin,m/dt/length(np1),'r');
    
    % 5p
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},np5,'window',w,'dt',dt);
      
    m = histc(a, tbin);
    plot(tbin,m/dt/length(np5),'g');

    % cue
    
    subplot(212)
    
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},all_cue,'window',w,'dt',dt);
    
    tbin = linspace(w(1), w(2), diff(w)/dt+1);
    
    m = histc(a, tbin);
	bar(tbin,m/dt/length(all_cue));
    
    hold on;
    
    % 1p
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},evt.c1,'window',w,'dt',dt);
     
    m = histc(a, tbin);
    plot(tbin,m/dt/length(evt.c1),'r');
    
    % 5p
    w = [-2 5]; dt = 0.1;
    [a,b] = spikePETH(S{iC},evt.c5,'window',w,'dt',dt);
      
    m = histc(a, tbin);
    plot(tbin,m/dt/length(evt.c5),'g');
    
    pause;
    clf;
    
    
end