function [ C_decon,spikes,noise_vals ] = trial_wise_deconvolve(C,sortedIndices,num_trials)
% This function deconvolves calcium from raw fluorescence traces using
% OASIS with an AR1 model
%   Uses parameters (such as time-constant, baseline, etc.) estimated from each cell's entire trace to
%   deconvolve individual calcium traces to obtain calcium estimate C_decon
%   and sparse spiking activity spikes. Also returns vector of noise
%   standard deviation estimated for every neuron


C = C - min(C,[],2); % subtract minimum

num_neurons = size(C,1);
C_decon = zeros(size(C));
spikes = zeros(size(C));
noise_vals = zeros(num_neurons,1);

trial_idx = unique(sortedIndices(:,1));

for neuron = 1:num_neurons
    
    if ~isempty(find(isnan(C(neuron,:)), 1))
        fprintf('Skipping deconvolution of neuron %d due to NaN entries\n',neuron);
    else
       
        fprintf('Currently deconvolving %d of %d total neurons\n',neuron,num_neurons);
        [~,~,params] = deconvolveCa(C(neuron,:),'ar1','optimize_pars',1,'optimize_b',1,'method','thresholded');
        
        noise_vals(neuron) = params.sn;
        
        for trial = 1:num_trials
            trial_id = trial_idx(trial);
            trial_inds = find(sortedIndices(:,1) == trial_id);
            trial_data = C(neuron,trial_inds);
            [C_decon(neuron,trial_inds),spikes(neuron,trial_inds),~] = deconvolveCa(trial_data,'ar1','pars',params.pars,'b',params.b,'sn',noise_vals(neuron),'method','thresholded');
        end
    end
    
end

end

