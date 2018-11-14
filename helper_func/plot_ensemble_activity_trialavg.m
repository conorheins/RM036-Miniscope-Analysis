function  plot_ensemble_activity_trialavg(est_factors,data2use,trial_length,num_trials,trial_indices,assembly_info,event_names,events_reshaped,events2show)
%plot_ensemble_activity Plot joint activity of different neuronal ensembles
% identified with TCA
%   Plot trial-averaged activity of neuronal ensembles identified with TCA


assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,~] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    membership_matrix(neuron,assembly_info(neuron).membership) = 1;
end

event_tmsp = zeros(1,length(events2show));
for i = 1:length(events2show)
    event_idx = find(cellfun(@(x) ~isempty(x), strfind(event_names,events2show{i})));
    [~,event_tmsp(i)] = max(sum(events_reshaped(event_idx,:,:),3));
end


figure('Position',[200 400 1000 800])

for i = 1:length(unique_ensembles)
    subplot(length(unique_ensembles),1,i)
    temp_ensemble_neurons = data2use(membership_matrix(:,i) == 1,:);
    trial_reshaped = reshape(temp_ensemble_neurons,size(temp_ensemble_neurons,1),trial_length,num_trials);
    trial_reshaped = trial_reshaped(:,trial_indices,:);
    ensemble_recon = re_estimate_latent(trial_reshaped,est_factors.U{1}(membership_matrix(:,i) == 1,i),est_factors.U{3}(:,i));
    plot(ensemble_recon./max(ensemble_recon),'LineWidth',1.5);
    hold on;
    plot(est_factors.U{2}(:,unique_ensembles(i))./max(est_factors.U{2}(:,unique_ensembles(i))),'LineWidth',1.5);
    limz = ylim;
    for jj = 1:length(events2show)
        plot([event_tmsp(jj) event_tmsp(jj)],[0 1],'-','LineWidth',2,'DisplayName',events2show{jj});
        hold on;
    end
    axis tight
    title(sprintf('Estimated and reconstructed activity for ensemble %d',unique_ensembles(i)),'FontSize',15)
%     if i == 1
%         legend('show')
%         legend_names = [{'Average estimated ensemble activity','Estimated trial trajectory'},events2show];
%         h=legend(legend_names);
%         %     rect=[.85 (1 - i/length(unique_ensembles)) .05 .05];
%         %     set(h,'Position',rect)
%         set(h,'Location','northeastoutside');
%         
%     end

    legend('show')
    legend_names = [{'Average estimated ensemble activity','Estimated trial trajectory'},events2show];
    h=legend(legend_names);
    set(h,'Location','northeastoutside');
    
end

