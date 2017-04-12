function d = diffang_rad(a1,a2)
%

d = a1 - a2;

idx = find(d < -pi);
d(idx) = d(idx) + 2*pi;

idx = find(d >= pi);
d(idx) = d(idx) - 2*pi;
