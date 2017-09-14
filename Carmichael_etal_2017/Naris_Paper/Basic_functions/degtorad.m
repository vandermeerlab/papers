function angleInRadians = degtorad(angleInDegrees)
% DEGTORAD Convert angles from degrees to radians
%
%   angleInRadians = DEGTORAD(angleInDegrees) converts angle units from
%   degrees to radians.
%
%   Example
%   -------
%   % Compute the tangent of a 45-degree angle
%   tan(degtorad(45))
%
%   See also: fromDegrees, fromRadians, toDegrees, toRadians, radtodeg.

% Copyright 2009 The MathWorks, Inc.

angleInRadians = (pi/180) * angleInDegrees;
