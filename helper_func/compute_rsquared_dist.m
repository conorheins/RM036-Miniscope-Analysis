function [ neurons_rsquared ] = compute_rsquared_dist(est_factors,data)
%compute_rsquared_dist Given a tensor factorization model est_factors and
%some data, reconstruct the data and see how well it recapitulates single
%neuron activity.
%   Detailed explanation goes here


n_neurons = size(data,1); n_timesteps = size(data,2); nTrials = size(data,3);

data = reshape(double(data),n_neurons,n_timesteps*nTrials);

recon = reshape(double(full(est_factors)),n_neurons,n_timesteps*nTrials);

neurons_rsquared = zeros(n_neurons,1);

for neuron = 1:n_neurons
    neurons_rsquared(neuron) = corr(data(neuron,:)',recon(neuron,:)').^2;
end

