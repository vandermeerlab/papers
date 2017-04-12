function [csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm, xc_cond, FORCE)
% [csd, CSDelecinds] = csd5pt(xcpots, elec_sep_mm, xc_cond*, FORCE*)
% (* = optional argument)
%
% 5-point estimation of the CSD using linear electrode recordings, using a
% Hamming window smoothing (Reference: Rappelsberger et al., Pfluegers
% Archiv 1981). FORCE flag prevents function from treating larger dimension
% as time dimension. If xc_cond is set to 1, csd is returned in units of
% mV/mm^2.
%
% Created by E.W.Schomburg, March 2014.

if (nargin < 3) || isempty(xc_cond) || (xc_cond <= 0)
	xc_cond = 1/3333; % S/mm = 3e-3 S/cm
end
if (nargin < 4)
    FORCE = 0;
end

if (size(xcpots,2) > size(xcpots,1)) && ~FORCE
    xcpots = xcpots';
    DIMFLIP = true;
else
    DIMFLIP = false;
end
N = size(xcpots,2);

% There also must be at least five electrodes
if (N < 5)
    error('Voltages from at least five electrodes needed')
end

% Create 2nd order derivative matrix; for N electrodes, this is size (N-4) x N.
% The square of the electrode separation is in the denominator.
D = zeros(N-4,N);
D(1,1:5) = [0.23 0.08 -0.62 0.08 0.23]/elec_sep_mm^2; % mm^-2 (Rappelsberger et al, 1981)
for i=2:N-4
    D(i,:) = circshift(D(i-1,:),[0,1]);
end

% CSD = -(conductivity)*(2nd spatial derivative of voltage along electrode line)
CSDelecinds = 3:(size(xcpots,2)-2);
csd = (-xc_cond*D*xcpots')'; % (ohm-mm)^-1*(mV/mm^2) = mA/mm^3

if DIMFLIP
    csd = csd';
end
