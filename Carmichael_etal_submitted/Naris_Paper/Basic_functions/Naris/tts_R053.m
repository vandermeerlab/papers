function tetrodes = tts(tetrode);
% functional purpose: to call up the channels for requested tetrode

list_of_channels{1} = [22 24 26 28];
list_of_channels{2} = [21 23 25 27];
list_of_channels{3} = [14 16 18 20];
list_of_channels{4} = [13 15 17 19];
list_of_channels{5} = [6 8 10 12];
list_of_channels{6} = [5 7 9 11];
list_of_channels{7} = [1 2 3 4];
list_of_channels{8} = [67 68 69 70]; %fill space
list_of_channels{9} = [54 56 58 60];
list_of_channels{10} = [53 55 57 59];
list_of_channels{11} = [46 48 50 52];
list_of_channels{12} = [45 47 49 51];
list_of_channels{13} = [38 40 42 44];
list_of_channels{14} = [37 39 41 43];
list_of_channels{15} = [30 32 34 36];
list_of_channels{16} = [29 31 33 35];

tetrodes.channels_array = [];

for iT = 1:1:length(tetrode)
    tetrodes.identity{iT} = tetrode(iT);
    tetrodes.channels{iT} = list_of_channels{tetrode(iT)};
    tetrodes.channels_array = cat(2, tetrodes.channels_array, tetrodes.channels{iT});
end

