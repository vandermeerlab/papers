function [ahv,hd_unwrapped] = HDtoAHVnoflip(hd)
% function [ahv,hd_unwrapped] = HDtoAHVnoflip(hd)
%
flip_threshold = 180;

ahv = zeros(1,length(hd.data)-1);
d = zeros(1,length(hd.data)-1);

for iT = 2:length(hd.data)
    
    d(iT-1) = diffang(hd.data(iT),hd.data(iT-1));
    
    if d(iT-1) > flip_threshold
        hd.data(iT) = add360(hd.data(iT),180);
        d(iT-1) = add360(d(iT-1),180);
    end
    
	ahv(iT-1) = d(iT-1)/(hd.tvec(iT)-hd.tvec(iT-1));


end

if nargout == 2
    d = [hd.data(1) d];
    hd_unwrapped = cumsum(d);
    hd_unwrapped = tsd(hd.tvec,hd_unwrapped);
end

ahv = [0 ahv];
ahv = tsd(hd.tvec,ahv);