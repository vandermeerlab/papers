function out = add2pi(in1,in2)
% function out = add2pi(in1,in2)
%
% adds two angles on a {-pi,pi} domain together
% in1 can be a vector

out = zeros(length(in1),1);
for i = 1:length(in1)

if (in1(i)+in2) > pi
	over = (in1(i)+in2)-pi;
	out(i) = -pi+over;
elseif (in1(i)+in2) < -pi
	under = in1(i)+in2+pi;
	out(i) = pi+under;
else
	out(i) = in1(i)+in2;
end

end