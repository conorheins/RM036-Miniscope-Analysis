function [ lower_tri_vect ] = get_lower_tri( matrix )
%get_lower_tri Pulls out lower-triangle of a matrix. Super useful

    A = 1:length(matrix);
    B = rot90(bsxfun(@plus,A,A(:)))-1;
    lower_tri_vect = matrix(B<length(matrix));

end

