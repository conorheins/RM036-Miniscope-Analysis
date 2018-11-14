function [cleaned_cell_to_index_map,idx2del] = calculate_global_excludeIDX( cell_to_index_map,badneur_array )
%calculate_global_excludeIDX Cell array of day-specific neurons to delete
% from global mapping
%   Detailed explanation goes here

num_days = size(cell_to_index_map,2);
idx2del = [];
for day = 1:num_days
    [~,temp] = ismember(badneur_array{day},cell_to_index_map(:,day));
    idx2del = [idx2del,temp];
end

idx2del = unique(idx2del);
cleaned_cell_to_index_map = cell_to_index_map;
cleaned_cell_to_index_map(idx2del,:) = [];

end

