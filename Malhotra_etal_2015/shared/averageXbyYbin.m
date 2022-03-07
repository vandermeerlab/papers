function x_avg = averageXbyYbin(x,y,y_edges)
%function x_avg = averageXbyYbin(x,y,y_edges)

[~,idx] = histc(y,y_edges);

x_avg = zeros(size(y_edges));
for iBin = length(y_edges):-1:1
    
    if sum(idx == iBin) ~= 0 % at least one sample in this bin
        x_avg(iBin) = nanmean(x(idx == iBin));
    end
    
end
x_avg = x_avg(1:end-1);

% fast but memory-intensive
%bins = min(max(idx,1),length(y_edges)-1); % eliminate last edge
%ypos = ones(size(bins,1),1);
%ns = sparse(bins,ypos,1);
%xsum = sparse(bins,ypos,x);
%x_avg = full(xsum)./(full(ns));


