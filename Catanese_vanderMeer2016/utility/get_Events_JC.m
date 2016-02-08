function [events] = get_Events_JC (dataFolder, varargin)
% [events] = get_Events_JC (dataFolder, varargin)
% output a structure events containing the time of each events
% varagin = 1 to plot; default = no plot;
% JC 01 Novembre 2013

if nargin == 1
    Plot_evts = 0 ;
elseif nargin == 2
    Plot_evts = varargin{1} ;
end

dataFolder
cd(dataFolder);

fn = FindFile('*.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);
% EVTimeStamps_converted = 10^-6*(EVTimeStamps - EVTimeStamps(1));

events = [];

%% Create List of interestings events :
eventList = {...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000)','p0',...  % pokeOUT
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0001)','p1',...  % npoke1
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0002)','p2',...  % npoke2
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004)','p3',...  % npoke3
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0008)','p4',...  % npoke4
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0010)','p5',...  % Centralpoke5
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0000)','f0',... % FeederOUT
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0001)','f1',... % Feeder1
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0002)','f2',... % Feeder2
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004)','f3',... % Feeder3
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0008)','f4'}; % Feeeder4

%% Create events structures


for iE = 1:2:length(eventList),
    
    ev_string = eventList{iE};
    ev_target = eventList{iE+1};
    
    ev_id = strncmp(ev_string,EventStrings,length(ev_string));
    ev_t = EVTimeStamps(ev_id);
    
    events.(ev_target) = ev_t;
end

%%
if Plot_evts == 1
    
    [Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT(FindFile('*.nvt'), [1 1 1 1 1 1], 1, 1, [] );
%     Timestamps = 10^-6*(Timestamps -  EVTimeStamps(1));
    
    %%
    figure(10), close(10), figure(10),
    toPlot = {'p1', 'p2', 'p5'};
    cols = {[1 0 0],[0 1 0] ,[0 0 1]};
    
    plot(X,Y,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
    hold on;
    for iP = 1:length(toPlot)
        
        temp = getfield(events,toPlot{iP});
        temp_x = interp1(Timestamps,X,temp,'linear');
        temp_y = interp1(Timestamps,Y,temp,'linear');
        
        h = plot(temp_x,temp_y,'o','MarkerSize',20);
        set(h,'MarkerFaceColor',cols{iP},'MarkerEdgeColor',cols{iP});
    end
end


% if Plot_evts == 1
%     pos = FindFiles('*vt.mat');
%     load(pos{1});
%     figure(1), close(1), figure(1),
%     toPlot = {'p1', 'p2', 'p5'};
%     cols = {[1 0 0],[0 1 0] ,[0 0 1]};
%
%     plot(Data(x),Data(y),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
%     hold on;
%     for iP = 1:length(toPlot)
%
%         temp = getfield(events,toPlot{iP});
%         temp_x = interp1(Range(x),Data(x),temp,'linear');
%         temp_y = interp1(Range(y),Data(y),temp,'linear');
%
%         h = plot(temp_x,temp_y,'o','MarkerSize',20);
%         set(h,'MarkerFaceColor',cols{iP},'MarkerEdgeColor',cols{iP});
%     end
%
%
% end



