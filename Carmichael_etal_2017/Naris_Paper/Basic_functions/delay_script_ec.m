function delay_script_ec(time, script_to_run)
%% runs a script at a later time.
% inputs:
%    time {HH:MM} = time you want it to run example time = '17:30'
%    script_to_run {str}


%%
x= datestr(now);
time_hr = x(end-8:end-6);
time_min = x(end-4:end-3);
fprintf(['Current time: ' x(end-8:end-3) '\nTarget time: ' time(1:2) ':' time(3:4) '\n'])
while str2double(time_hr) <= (str2double(time(1:2))) && str2double(time_min) <=str2double(time(4:5)) ==1
    disp(datestr(now))
    pause(60)
    x= datestr(now);
    time_hr = x(end-8:end-6);
    time_min = x(end-4:end-3);
end
%%
disp('finished')

disp(script_to_run)
end