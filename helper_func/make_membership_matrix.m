function [ mem_matrix ] = make_membership_matrix( assembly_info )
%make_membership_matrix: Turn assembly info struct into membership matrix
%   Detailed explanation goes here

assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,~] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

mem_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    mem_matrix(neuron,assembly_info(neuron).membership) = 1;
end


end

