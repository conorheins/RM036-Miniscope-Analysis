function [entropies, spike_time_var  ] = compute_var(trial_reshaped,binWidth)
%compute_var Compute spike-timing variability 
%  uses two metrics -- entropy of different spike-timing states (including the
%  zero state/no firing), plus spike-timing variability during those
%  patterns that it did indeed fire

max_matrix = bsxfun(@eq,trial_reshaped,max(trial_reshaped,[],2)); % find all the times that the neuron was firing maximally

% clean up max_matrix 

for ii = 1:size(max_matrix,3)
    for jj = 1:size(max_matrix,1)        
        where_elements = find(max_matrix(jj,:,ii));        
        if numel(where_elements) == size(max_matrix,2)
            max_matrix(jj,:,ii) = 0;
        end        
        where_elements = find(max_matrix(jj,:,ii));        
        if numel(where_elements) > 1
            max_matrix(jj,where_elements(2:end),ii) = 0;
        end        
    end
end

entropies = zeros(1,size(max_matrix,1));
spike_time_var = zeros(1,size(max_matrix,1));

numtimesteps = size(max_matrix,2);
T = floor(numtimesteps / binWidth);

for neuron = 1:size(max_matrix,1)
    
    firing_pattern_dist = mean(max_matrix(neuron,:,:),3);
    
    % bin probabilities
    
    binned = zeros(1,T);    
    for t = 1:T       
      iStart = binWidth * (t-1) + 1;
      iEnd   = binWidth * t;      
      binned(t) = sum(firing_pattern_dist(iStart:iEnd));
    end
    
    prob_dist = [binned, (1-sum(binned))]; % add probability of no firing probability (complement of total firing probability)
    P = prob_dist(prob_dist > 0); % only include positive probabilities (due to -Inf when you call log(0))
    entropies(neuron) = -dot(P,log(P)); % calculation of entropy of spike-pattern distribution
    
    firing_times = squeeze(max_matrix(neuron,:,:));
    firing_times = firing_times(:,sum(firing_times,1) > 0);
   
    [~,maxInds] = max(firing_times,[],1);
    
    spike_time_var(neuron) = std(maxInds); % standard deviation of spike times
    
end
    

end

