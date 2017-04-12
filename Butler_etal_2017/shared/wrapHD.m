function hd_wrapped = wrapHD(hd_unwrapped)
% function hd_wrapped = wrapHD(hd_unwrapped)
%

hd_wrapped = nan(size(hd_unwrapped));
plus_idx = hd_unwrapped >= 0;
hd_wrapped(plus_idx) = rem(hd_unwrapped(plus_idx),360);
minus_idx = hd_unwrapped < 0;
hd_wrapped(minus_idx) = rem(hd_unwrapped(minus_idx),-360)+360;