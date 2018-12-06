function [ lower_tri_vect ] = get_lower_tri( matrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    A = 1:length(matrix);
    B = rot90(bsxfun(@plus,A,A(:)))-1;
    lower_tri_vect = matrix(B<length(matrix));

end

