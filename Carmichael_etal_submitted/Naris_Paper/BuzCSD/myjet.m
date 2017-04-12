% MYJET             modified jet with extreme values in pure R,B.
%
% call              MAP = MYJET( M )
%
% gets              M       colormap length
%
% returns           MAP     the colormap ( M x 3 )
%
% does              creates a colormap where maximal values are red (not
%                   brown) and minimal values are blue (and not black-blue)
%
% calls             nothing
%
% called by         NANPAD

% 17-apr-04 ES

function map = myjet( m )

if nargin < 1 | isempty( m ), m = 64; end

n = ceil( ( m - 2 ) / 3 );
u = [ ( 1 : 1 : n ) / n  ones( 1, n - 2 ) ( n : -1 : 1 ) / n ]';
g = ( 1 : length( u ) )';
r = g + n + 1;
b = g - floor( 3 * n / 4 ) + 1;
g = g + 2;
g( g > m ) = [];
r( r > m ) = [];
b( b < 1 ) = [];
map = zeros(m,3);
map( r, 1 ) = u( 1 : length( r ) );
map( g, 2 ) = u( 1 : length( g ) );
map( b, 3 ) = u( end - length( b ) + 1 : end );

return