function [tc_corr,tc_shift] = TCcorr(tc1,tc2)
% function [tc_corr,tc_shift] = TCcorr(tc1,tc2)
%
% output is how much to shift tc1 to be most correlated with tc2


nCells = size(tc1,1);
nBins = size(tc1,2);

all_shift = 0:nBins-1;
for iShift = length(all_shift):-1:1
    
    this_shift = all_shift(iShift);
    this_tc1 = circshift(tc1,[0 this_shift]);
    
    tc1_data = this_tc1(:); tc2_data = tc2(:);
    keep_idx = find(~isnan(tc1_data) & ~isnan(tc2_data));
    
    this_corr = corrcoef(tc1_data(keep_idx),tc2_data(keep_idx));
    all_corr(iShift) = this_corr(1,2);
    
end

[tc_corr,idx] = max(all_corr);
tc_shift = all_shift(idx);
