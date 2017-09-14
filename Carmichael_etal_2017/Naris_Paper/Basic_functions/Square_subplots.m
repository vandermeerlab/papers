function Square_subplots()    
% makes subplots square.


axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
end