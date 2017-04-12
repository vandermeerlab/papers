function angleInDegrees = rad2deg(angleInRadians)

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2012/01/16 16:51:28 $

%#codegen

angleInDegrees = (180/pi) * angleInRadians;
