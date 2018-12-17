function [ binned_Rtrials,binned_NRtrials,bin_centers ] = bin_distances(distances_R_trials,R_times,distances_NR_trials,NR_times,bin_width)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

T = max([max(R_times), max(NR_times)]);

bin_edges = 0:bin_width:T;
bin_edges(end) = T;

binned_Rtrials = zeros(1,length(bin_edges)-1);
binned_NRtrials = zeros(1,length(bin_edges)-1);

for bin_i = 1:length(bin_edges)-1
    
    curr_bin = bin_edges(bin_i):bin_edges(bin_i+1);
    
    [~,Rbin_idx] = ismember(curr_bin,R_times);
    Rbin_idx = Rbin_idx(Rbin_idx~=0);
    
    binned_Rtrials(bin_i) = mean(distances_R_trials(Rbin_idx));
    
    [~,NRbin_idx] = ismember(curr_bin,NR_times);
    NRbin_idx = NRbin_idx(NRbin_idx~=0);
    
    binned_NRtrials(bin_i) = mean(distances_NR_trials(NRbin_idx));
    
end

bin_centers = diff(bin_edges)/2 + bin_edges(1:end-1);
    

end

