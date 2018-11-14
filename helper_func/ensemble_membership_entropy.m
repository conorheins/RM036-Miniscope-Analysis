function ensemble_membership_entropy(trial_reshaped,binWidth,assembly_info)
%ensemble_membership_entropy Relate neuronal ensemble membership to
%spike-timing variability/entropy


assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,~] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    membership_matrix(neuron,assembly_info(neuron).membership) = 1;
end

[entropies,~] = compute_var(trial_reshaped,binWidth);

figure;

colors = flag(length(unique_ensembles));

for i = 1:length(unique_ensembles)
    
    ensemble_members = find(membership_matrix(:,unique_ensembles(i)) == 1);
    
    scatter(entropies(ensemble_members),[assembly_info(ensemble_members).r_sq],35,repmat(colors(i,:),length(ensemble_members),1),...
        'DisplayName',sprintf('Ensemble No. %d',unique_ensembles(i)));
    xlabel('Entropy of within-trial firing patterns')
    ylabel('R^2 of assembly membership');
    title('Entropy vs. R^2 for Different ensembles');

    pause; hold on;
end

legend('show')

    
    
    
    
    



end

