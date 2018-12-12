function [ binned ] = bin_patterns_overTime( pattern_matrix,bin_size )
%bin_patterns_overTime This function takes a pattern matrix (N x T) and
%takes a moving window average of adjacent patterns
%   Inputs: pattern_matrix: N x T matrix of patterns, where N is the number
%           of variables in one pattern, and T is the number of patterns (or
%           temporal 'slices')
%           bin_size: number of adjacent patterns to average into a single pattern  

[n,T] = size(pattern_matrix);

bin_segs = 1:bin_size:T;

binned = zeros(n,length(bin_segs)-1);
for seg_i = 1:length(bin_segs)-1
    binned(:,seg_i) = mean(pattern_matrix(:,bin_segs(seg_i):bin_segs(seg_i+1)),2);
end
    

end

