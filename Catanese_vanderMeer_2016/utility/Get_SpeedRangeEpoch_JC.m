function [Epoch] = Get_SpeedRangeEpoch_JC(DataFolder, epoch, MinSpeed, MaxSpeed, FixSpeed, SecxBin, varargin)
% function [Tps Idx Epoch] = Get_SpeedRangeTimes_JC(DataFolder, epoch, MinSpeed, MaxSpeed, SecxBin, varargin)
% INPUT : define - Epoch = One of this 3 strings = 'Task' , 'PostTask', 'PreTask';
%                - SpeedRange = min and max speed values in cm/sec : (eg 0 < speed < 1 cm/s)
%                - FixSpeed = 'MinSpeed' or 'MaxSpeed'; if Min it will reduce Max until less than 30% of task selected.  
%                - SecxBin = 0.2 (second per bin) will determine bin size.  
%                - varargin = plot or not;  default = 0; type 1 for display plots;  
% OUTPUT: return - the Times (Timestamps.* 10^-6) corresponding to the defined speed range.
%                - the Index from Data(x)(position tsd).
% by J. Catanese, Dec 2013, mvdmlab.

AddStd = 1

if nargin == 5
    fplot = 0
else
    fplot = varargin{1}
end

cd(DataFolder); % eg :  DataFolder = 'D:\Julien_VdmLab\Data\incoming\R039-2013-08-15'
load('Epoch_list.mat');

%% Get Position Data
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT(FindFile('*.nvt'), [1 1 1 1 1 1], 1, 1, [] );

Timestamps = Timestamps .* 10^-6;
x = tsd(Timestamps,X');
y = tsd(Timestamps,Y');

Xd = Data(x);
Xt = Range(x);
Yd = Data(y);
Yt = Range(y);

%% Remove and interpolate GAP due to loose of tracking (x=0 y=0)
GAP = find(Xd == 0 & Yd == 0);
XdGAP= Xd;
XtGAP=Xt;
XdGAP(GAP) = [];
XtGAP(GAP) = [];
disp(['global percentage of GAP = ' num2str(round(max(size(GAP))/max(size(Xd))*100)) '%'])

XdNEW = interp1(XtGAP, XdGAP, Xt);

YdGAP= Yd;
YtGAP=Yt;
YdGAP(GAP) = [];
YtGAP(GAP) = [];

YdNEW = interp1(YtGAP, YdGAP, Yt); % linear interp
% figure, subplot(122), plot(XdNEW, YdNEW);  subplot(121), plot(Xd, Yd);

x_new = tsd(Xt, XdNEW');
y_new = tsd(Yt, YdNEW');

%% Restrict to the selected Epoch
x = Restrict(x_new, eval(['Epoch.' epoch '.Start']), eval(['Epoch.' epoch '.End']));
y = Restrict(y_new, eval(['Epoch.' epoch '.Start']), eval(['Epoch.' epoch '.End']));

%% Get the Speed.tsd
s = GetLinSpd(x,y);
Sd = Data(s) ;
St= Range(s) ;

%% Calculate mean speed over Xsec bin
Fs = 1/mean(diff(St))
binsize = round(SecxBin*Fs) ;
binvec = [1:binsize:max(size(St))];
% centerbin = [];
Nbin = max(size(binvec));
Sbin = NaN(1,Nbin); SdevS = NaN(1,Nbin);
for ib = 1:Nbin-1
    Sbin(ib) = mean(Sd(binvec(ib):binvec(ib+1)));
    SdevS(ib) = std(Sd(binvec(ib):binvec(ib+1))) ;
end

if fplot
    figure, hold on,
    plot(St,Sd, 'g');
    plot(St(binvec), Sbin, 'r');
    plot(St(binvec), Sbin+SdevS, '--k');
    legend('inst speed',['mean speed (' num2str(SecxBin) 's/bin)'],'std dev')
end

if AddStd == 1
    Sbin = Sbin+SdevS ;
end

%% Get Idx and Tps in the speed range (OUTPUT)
Idx = find(Sbin >= MinSpeed & Sbin <= MaxSpeed);
Dat = Sbin(Idx);
Tps = St(binvec(Idx))';

%% Adjust the treshold (MaxSpeed) in order to select around 20% of the time with lower speed.
MaxSpeed2=MaxSpeed;
MinSpeed2=MinSpeed;

W =0;  
while length(Idx) >= 0.3*length(Sbin);
    W=W+1
    if FixSpeed == 'MinSpeed'
        MaxSpeed2 = MaxSpeed2-0.2  
        Idx = find(Sbin >= MinSpeed2 & Sbin <= MaxSpeed2);
    elseif FixSpeed == 'MaxSpeed'
        MinSpeed2 = MinSpeed2+0.2  
        Idx = find(Sbin >= MinSpeed2 & Sbin <= MaxSpeed2);
    end
end
MinSpeed=MinSpeed2;
MaxSpeed = MaxSpeed2;
Dat = Sbin(Idx);
Tps = St(binvec(Idx))';

%% Find the Edges (Epoch.start and Epoch.end) to Restrict the Data :
% Just to understand you can draw a simple example :
% T =  0 0 0 x x x x 0 0 0 x x x 0 0 0 0       (x = missing data)
% Q =  1 2 3         4 5 6       7 8 9 10
% D =  1 1 5 1 1 4 1 1 1
% A =  [3 6]
% this is the epoch we want :
% Q( 1 : A(1) )
% Q( A(1)+1 : A(2))
% Q( A(2)+1 : end )

%% define Variables
CRIT = 1
T = Tps;
Q = find(T);
D = diff(Idx); % (could have use T himself but then have to convert my criterion in Time)
A = find(D > CRIT );


%% initialize
start = NaN(length(A)+1,1) ;
fin = NaN(length(A)+1,1);

start(1) = T(Q(1));
fin(end) = T(Q(end)) ;

%% Get Epoch.start & Epoch.end
for ai = 1:length(A)
    start(ai+1) = T(Q(A(ai)+1)) ;
    fin(ai) = T(Q(A(ai)));
end

clear Epoch
Epoch.Start = start;
Epoch.End = fin;
Epoch.EpochID = [epoch '-Speed-' num2str(MinSpeed) '-to-' num2str(MaxSpeed)]
Epoch.MinSpeed = MinSpeed
Epoch.MaxSpeed = MaxSpeed
Epoch.SecxBin = SecxBin
Epoch.speedunit = 'pixel/sec'

%% Restrict and Plot speed

A=[];
Nepoch = max(size(Epoch.Start));
figure, hold on,
plot(s,'Color', [0.9 0.9 0.9]),
plot(St(binvec), Sbin, 'r'),
plot(St(binvec), ones(1,max(size(binvec)))*MaxSpeed, '--k'),
plot(St(binvec), ones(1,max(size(binvec)))*MinSpeed, '-.k'),

for ie= 1:Nepoch
    A = [A Restrict(s, Epoch.Start(ie), Epoch.End(ie))];
    plot(A(ie), 'Color', 'g'), hold on;
end

if AddStd == 1
    legend('inst speed',['mean speed + Std (' num2str(SecxBin) 's/bin)'], 'treshold Max', 'treshold Min' , 'selected speed' )
else
    legend('inst speed',['mean speed (' num2str(SecxBin) 's/bin)'], 'treshold' , 'selected speed' )
end
% plot(St(binvec), Sbin, 'r'),
% plot(St(binvec), ones(1,max(size(binvec)))*MaxSpeed, '--k'),
% plot(St(binvec), ones(1,max(size(binvec)))*MinSpeed, '-.k'),

title([ num2str(round(length(Idx)/length(Sbin)*100)) '% of ' epoch ' with speed from '  num2str(MinSpeed)  ' to ' num2str(MaxSpeed) 'px/s (= ' num2str(round(length(Idx)*SecxBin/60)) 'min)' ])

if fplot == 1
    saveas(gcf,[DataFolder '\Figures_jpg\Select_' Epoch.EpochID 'PxS_' num2str(Epoch.SecxBin) 'SxB.png'])
end

%% Plot (Optional)
Xd = Data(x_new);
Yd = Data(y_new);
Xd2 = Data(Restrict(x, Epoch.Start,  Epoch.End )) ;
Yd2 = Data(Restrict(y, Epoch.Start,  Epoch.End)) ;

if fplot == 1;
    figure, hold on,
    plot(Yd,Xd,'.');hold on,
    plot(Yd2,Xd2,'.r');
    legend('Full Session',[num2str(MinSpeed) '< px/s <' num2str(MaxSpeed) '==>'  num2str(round(max(size(Xd2))/max(size(Xd))*100)) '%']);
end

