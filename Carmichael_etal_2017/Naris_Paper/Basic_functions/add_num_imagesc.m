function [h, hStrings] = add_num_imagesc(h, x, n_dec, font_size)
%% add_num_imagesc:  adds the numerical values of a matrix used to generate an imagesc.
%
%          Inputs:
%           - h: axis handle
%           - x: matrix
%           - font_size: default is what ever the figure already had
%          Outputs:
%           - h: axis handle
%
% EC - 2016-10-31
if nargin <3
    n_dec = 3;
end
if nargin <4
    font_size = get(gca, 'fontsize');
end
if n_dec == 2
    textStrings = num2str(x(:),'%0.2f');  %# Create strings from the matrix values
elseif n_dec == 3
    textStrings = num2str(x(:),'%0.3f');  %# Create strings from the matrix values
elseif n_dec == 4
    textStrings = num2str(x(:),'%0.4f');  %# Create strings from the matrix values
elseif n_dec == 5
    textStrings = num2str(x(:),'%0.5f');  %# Create strings from the matrix values
end
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:length(x));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
set(hStrings,'fontsize',font_size)