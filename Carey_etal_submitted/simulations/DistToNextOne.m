function dist = DistToNextOne(data, idx, def)
% compute how many idxs away the next 1 is in data from idx (can be a vector)
% if no 1 found, set distance to def

for iI = length(idx):-1:1
    
    I = idx(iI);
    this_data = data(I:end);
    is_one = find(this_data == 1);
    
    if isempty(is_one) % no next one found, set to default
        dist(iI) = def;
    else
        dist(iI) = min(is_one);
    end
end
