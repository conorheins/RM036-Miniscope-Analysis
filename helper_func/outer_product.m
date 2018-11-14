function [ rank_one_tensor ] = outer_product( a,b,c )
%outer_product computes outer product of three vectors such that entry
%i,j,k of output tensor ('rank_one_tensor') is a-i * b-j * c-k
%   Detailed explanation goes here

[xx, yy, zz] = ndgrid(1:length(a), 1:length(b), 1:length(c));

rank_one_tensor = a(xx) .* b(yy) .* c(zz);

end

