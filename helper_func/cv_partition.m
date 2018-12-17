function CV_cell_array = cv_partition(response_variable,P,num_folds,class_names)
%cv_partition Wrapper function for cvpartition function in MATLAB
%   Returns cell array of training and test indices, relative to the original group 'response_variable' 


CV_cell_array = cell(1,num_folds);
C1_idx_relative = find(strcmp(response_variable,class_names{1}));
C2_idx_relative = find(strcmp(response_variable,class_names{2}));

% find the class that occurs the least, and ensure that other class is sampled with same probability in training/test sets
[~,rarest_class] = min([length(C1_idx_relative),length(C2_idx_relative)]); 

for f_i = 1:num_folds
    
    if rarest_class == 1
        fold_idx = randperm(length(C2_idx_relative),length(C1_idx_relative));
        group = [response_variable(C1_idx_relative);response_variable(C2_idx_relative(fold_idx))];
        full_relative_idx = [C1_idx_relative;C2_idx_relative(fold_idx)];
    elseif rarest_class == 2
        fold_idx = randperm(length(C1_idx_relative),length(C2_idx_relative));
        group = [response_variable(C1_idx_relative(fold_idx));response_variable(C2_idx_relative)];
        full_relative_idx = [C1_idx_relative(fold_idx),C2_idx_relative];
    end
    temp = cvpartition(group,'HoldOut',P);
    CV_cell_array{1,f_i} = {full_relative_idx(temp.training), full_relative_idx(temp.test)};
    
end


end

