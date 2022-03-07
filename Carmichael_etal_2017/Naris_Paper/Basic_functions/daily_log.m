function daily_log()
%% Wrap up script used to log daily progress and put it into the google project log.  
ii_l = 10;
formatOut = 'dddd dd mmmm yyyy';
day = datestr(datenum(date),formatOut);

today = cell(10,1); tomorrow = cell(10,1);
for ii = 1:ii_l
today{ii} = input('Today: ', 's');
if isempty(today{ii})
    break
end
end
for ii = 1:ii_l
tomorrow{ii} = input('Tomorrow: ', 's');
if isempty(tomorrow{ii})
    break
end
end
fprintf([day '\n' '\n'])
fprintf('Today: \n')
for ii= 1:length(today)
    if isempty(today{ii})
        continue
    end
    fprintf(['- ' today{ii} '\n'])
end
fprintf('\n')
fprintf('Tomorrow: \n')
for ii= 1:length(tomorrow)
    if isempty(tomorrow{ii})
        continue
    end
    fprintf(['- ' tomorrow{ii} '\n'])
end
fprintf('__________________________________________________________________________\n')
end