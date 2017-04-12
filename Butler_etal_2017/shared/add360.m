function out = add360(in1,in2)
% function out = add360(in1,in2)
%
% adds two angles on a {0,360} domain together
% in1 can be a vector

out = zeros(length(in1),1);
for i = 1:length(in1)

if (in1(i)+in2) > 360
	over = (in1(i)+in2)-360;
	out(i) = 0+over;
elseif (in1(i)+in2) < 0
	under = in1(i)+in2;
	out(i) = 360+under;
else
	out(i) = in1(i)+in2;
end

end