function [ C_smooth ] = smooth_dFF(C,sortedIndices,num_trials,Fs,window_length,win_overlap)
%smooth_dFF trial-wise smoothing of calcium data, uses sortedIndices to
% separate trial timestamps

C_smooth = zeros(size(C));

trial_idx = unique(sortedIndices(:,1));
for trial = 1:num_trials
    trial_id = trial_idx(trial);
    temp = C(:,sortedIndices(:,1)==trial_id);
    for neuron = 1:size(temp,1)
        C_smooth(neuron,sortedIndices(:,1)==trial_id) = locsmooth(temp(neuron,:),Fs,window_length,win_overlap); % moving window of 2 frames with 1 frame overlap
    end
end

end

