function [ seq ] = GPFA_preprocess(spikes,num_trials,trial_length,trial_indices,useSqrt)
%GPFA_preprocess Pre-processing of spike matrix for the GPFA analysis in the package
%'NeuralTraj' by Byron Yu & John Cunningham

% INPUTS:
%
% spikes         - array of size neurons x T with the neuron data (spike
%                   counts)
% num_trials    - number of trials to reshape spike array into, before splitting up into 
%                 trial-specific struct entries
%
% trial_length  - length of trial, in timebins/frames
%
% trial_indices - specific trial indices to look at (in case you don't want
%               to look at the beginning/ending of the trial, for instance
%
% useSqrt       - flag for whether to square-root transform the spike
%                 counts
%
% OUTPUTS:
%
% seq         - data structure, whose nth entry (corresponding to
%               the nth experimental trial) has fields
%                 trialId      -- unique trial identifier
%                 T (1 x 1)    -- number of timesteps
%                 y (yDim x T) -- neural data 

trial_reshaped = reshape(spikes,size(spikes,1),trial_length,num_trials);
trial_reshaped = trial_reshaped(:,trial_indices,:);

seq = [];

for n = 1:num_trials
    
    data = trial_reshaped(:,:,n);
    
    seq(n).trialId = n;
    seq(n).T       = size(data,2);
    
    if useSqrt
        seq(n).y   = sqrt(data);
    else
        seq(n).y = data;
    end

end

