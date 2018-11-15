function n = numel_cell_array(A)
    n = 0;
    for i=1:numel(A)
        if iscell(A{i})
            n = n + numel_cell_array(A{i});
        else
            n = n + numel(A{i});
        end
    end
end

