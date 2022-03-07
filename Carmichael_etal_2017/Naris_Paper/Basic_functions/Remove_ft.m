function Remove_ft()
% check if field trip is being used, if so then remove it from the path.

out_path = which('filtfilt');

if strfind(out_path,'fieldtrip')
    idx = strfind(out_path,'fieldtrip');
    disp('FT found in path, removing...')
    g = genpath(out_path(1:idx+9));
    rmpath(g)
end


% out_path = which('hann');
% 
% if strfind(out_path,'FieldTrip')
%     idx = strfind(out_path,'FieldTrip');
%     disp('FT found in path, removing...')
%     g = genpath(out_path(1:idx+9));
%     rmpath(g)
%     disp('FT found for the Hann function')
% end
% check again

out_path = which('filtfilt');

if strfind(out_path,'fieldtrip')
    idx = strfind(out_path,'fieldtrip');
    disp('FT found in path, removing...')
    g = genpath(out_path(1:idx+9));
    rmpath(g)
else
    disp('No FT in path')
end
% out_path = which('hann');
% 
% if strfind(out_path,'FieldTrip')
%     idx = strfind(out_path,'FieldTrip');
%     disp('FT found in path, removing...')
%     g = genpath(out_path(1:idx+9));
%     rmpath(g)
%     disp('FT found for the Hann function')
% end