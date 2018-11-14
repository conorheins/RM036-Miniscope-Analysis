function ensemble_recon = re_estimate_latent(trial_reshaped,neuron_weights,trial_weights)

% re_estimate_latent: re-estimate within-trial trajectories using data and
% other two latent factors.
% INPUTS: trial_reshaped: N x TrialLength x nTrials tensor
%         neuron_weights: N x 1 neuron weight vector (used as
%                           neuron-specific weights in the average)
%         trial_weights: nTrials x 1 trial weight vector (used as
%                           trial-specfic weights in the average)

ensemble_recon = zeros(1,size(trial_reshaped,2));

for trial = 1:length(trial_weights)
    ensemble_recon = ensemble_recon + trial_weights(trial)*sum(squeeze(trial_reshaped(:,:,trial)).*repmat(neuron_weights,1,size(trial_reshaped,2)),1);
end