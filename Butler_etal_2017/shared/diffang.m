function d = diffang(a1,a2)
%

d = a1 - a2;

idx = find(d < -180);
d(idx) = d(idx) + 360;

idx = find(d >= 180);
d(idx) = d(idx) - 360;
