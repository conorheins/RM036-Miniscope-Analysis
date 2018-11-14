function [ boots_interval_array ] = bootci_array(nboot,func,nd_array,DIM)
% bootci_array Wrapper function to apply bootci (native MATLAB function) to
% ND arrays
%  Reshapes nd_array to put the dimension over which bootstrap statistic is computed
%  at the end 

dim_idx = 1:ndims(nd_array);
dim_idx(DIM) = []; dim_idx = [dim_idx,DIM]; % re-sort dim_idx (for use in permutation below)
nd_array = permute(nd_array,dim_idx);  % permute to put the dimension over which the statistics are computed at the end

array_dims = size(nd_array); % get sizes of each dimension, for initializing result matrix
array_dims(DIM) = []; % dimension that statistic is computed over will be replaced with 2-D matrix for boostrap intervals
array_dims = [array_dims,2];

boots_interval_array = zeros(array_dims);

% only handle 2-D and 3-D cases at the moment, need to expand into a
% recursive function at some point to handle ND_arrays of arbitrary
% dimension

if length(array_dims(1:end-1)) == 1
    for ii = 1:array_dims(1)
        boots_interval_array(ii,:) = bootci(nboot,{func,squeeze(nd_array(ii,:))},'type','per');
    end
elseif length(array_dims(1:end-1)) == 2
    for ii = 1:array_dims(1)
        for jj = 1:array_dims(2)
            boots_interval_array(ii,jj,:) = bootci(nboot,{func,squeeze(nd_array(ii,jj,:))},'type','per');
        end
    end
end

