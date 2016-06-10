function [hd,hd_unwrapped] = AHVtoHDfast(ahv,state)
% function [hd,hd_unwrapped] = AHVtoHDfast(ahv,state)
%
% ahv: input AHV tsd in deg/s
% state is state vector, [gain_l gain_r drift hd0 hd0_idx]

%cfg_def = [];
%cfg_def.gain = [1 1]; % L, R gain
%cfg_def.drift = 0; % in deg/s
%cfg_def.noisetype = 'none'; % 'none','gauss'
%cfg_def.noisesd = 1;
%cfg_def.hd0 = 0; % starting HD (deg)
%cfg_def.hd0_idx = []; % if specified, then force HD at this idx to be cfg_def.hd0 




dt = diff(cat(2,0,ahv.tvec));
dahv = zeros(size(ahv.data));

left_idx = ahv.data < 0;
right_idx = ahv.data >= 0;

dahv(left_idx) = dt(left_idx).*ahv.data(left_idx).*state(1);
dahv(right_idx) = dt(right_idx).*ahv.data(right_idx).*state(2);

dahv = dahv + dt*state(3);

% cumulative HD (without circular wrap-around)
hd_unwrapped = state(4) + cumsum(dahv);

% do the circwrap
hd_wrapped = nan(size(hd_unwrapped));
plus_idx = hd_unwrapped >= 0;
hd_wrapped(plus_idx) = rem(hd_unwrapped(plus_idx),360);
minus_idx = hd_unwrapped < 0;
hd_wrapped(minus_idx) = rem(hd_unwrapped(minus_idx),-360)+360;

% shift if requested
if state(5) ~= 0
   shft = hd_wrapped(state(5));
   hd_temp = hd_wrapped - shft;
   
   plus_idx = hd_temp >= 0;
   hd_wrapped(plus_idx) = rem(hd_temp(plus_idx),360);
   minus_idx = hd_temp < 0;
   hd_wrapped(minus_idx) = rem(hd_temp(minus_idx),-360)+360;
end

hd = hd_wrapped;

