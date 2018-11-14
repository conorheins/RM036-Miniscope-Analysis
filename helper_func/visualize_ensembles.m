function visualize_ensembles(data2use,assembly_info)
%visualize_ensembles Visualizes time traces of each ensemble in addition to
%their within-ensemble correlation matrices
%   INPUTS: data2use -- any N x T matrix of timeseries
%           assembly_info -- a structure array with two fields:
%           'membership' and 'r_sq'

assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,~] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    membership_matrix(neuron,assembly_info(neuron).membership) = 1;
end

shift_width = 1;

for i = 1:length(unique_ensembles)
    
    ens = unique_ensembles(i);
    
    ensemble_neurons = data2use(find(membership_matrix(:,i)),:);
    
    num_ens = size(ensemble_neurons,1);
    
    figure;
    subplot(121)
    max_norm = bsxfun(@rdivide,bsxfun(@minus,ensemble_neurons,min(ensemble_neurons,[],2)),(max(ensemble_neurons,[],2) -  min(ensemble_neurons,[],2)));
    shifts = repmat( (1:shift_width: (shift_width*num_ens))',1,size(max_norm,2));
    
    plot( (max_norm + shifts)','b-');
    title(sprintf('Single cell activities from ensemble no. %d',ens))
    
    subplot(122)
    imagesc(corrcoef(max_norm')); caxis([0 1]); colorbar;
    title(sprintf('Correlation matrix for neurons in ensemble no. %d',ens))
    
    pause; close gcf;
    
end

end

