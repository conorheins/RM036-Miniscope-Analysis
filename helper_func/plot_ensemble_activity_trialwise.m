function plot_ensemble_activity_trialwise(est_factors,trial_reshaped,events_reshaped,event_names,which_events,which_trials,assembly_info)
%plot_ensemble_activity Plot joint activity of different neuronal ensembles
% identified with TCA
%   Plot trial-by-trial plots of activity of neuronal ensembles identified with TCA

assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,srt] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    membership_matrix(neuron,assembly_info(neuron).membership) = 1;
end

[neural_mem,ensemble_ids] = find(membership_matrix);

event_names = event_names(which_events);

for ens = 1:length(unique_ensembles)
    ensembles = trial_reshaped(neural_mem(ensemble_ids == unique_ensembles(ens)),:,which_trials);
    for trial = 1:length(which_trials)
        subplot(121)
        imagesc(ensembles(:,:,which_trials(trial)));
        subplot(122)
        plot(est_factors.U{2}(:,ens)*est_factors.U{3}(which_trials(trial),ens),'DisplayName',sprintf('Ensemble # %d',ens)); hold on;
        this_trial_events = squeeze(events_reshaped(which_events,:,which_trials(trial)));
        [~,max_inds] = max(this_trial_events,[],2);
        [~,srt] = sort(max_inds,'ascend');
        this_trial_events = this_trial_events(srt,:);
        event_names = event_names(srt);
        limz = ylim;
        for event_i = 1:length(event_names)
            plot(repmat(find(this_trial_events(event_i,:),1),2,1),[0;limz(2)],'-','DisplayName',event_names{event_i})
        end
        legend('show');
        pause; hold off;
    end
end
    

