function [frame_by_event_matrix,empiricalFs,sortedIndices] = create_filter_matrix(sortedIndices,tmsp_folder,num_trials,trial_starts,event_timestamps_array)
% [create_filter_matrix] Creates logical matrix of frames X event_types

num_events = length(event_timestamps_array);
total_timestamps = numel_cell_array(event_timestamps_array);
frame_by_event_matrix = spalloc(size(sortedIndices,1),num_events,total_timestamps);
trialwise_Fs = zeros(num_trials,1);

% check if the timestamp folder is comprised of sub-folders, and if so,
% move all the timestamps within the subfolders into the main folder, and
% re-name appropriately


if isempty(dir(fullfile(tmsp_folder,'*.txt')))
   trial_idx = rename_timestamps(tmsp_folder);
else
   trial_idx = obtain_trial_idx(tmsp_folder);
end
    
for i = 1:length(trial_idx)
    
    trial = trial_idx(i);
    frame_grabs = trial_starts(trial);
    neural_mapper_name = fullfile(tmsp_folder,['image_timestamp_exp_',num2str(trial),'.txt']);
    NM_timestamps = read_in_timestamps(neural_mapper_name);
    trialwise_Fs(trial) = 1000 / mean(diff(NM_timestamps)); % frame-rate calculation, accounting for conversion from milliseconds to seconds
    
%     temp = [95:95:(95*799)]'; in case timestamps aren't available, use this artificial substitute and just assume 95 ms
%     between frames
    
    frame_grabs = [frame_grabs;trial_starts(trial) + NM_timestamps./1000];
    frame_grabs = frame_grabs(2:end);
    start_t = min(frame_grabs); end_t = max(frame_grabs);

    for event = 1:num_events
        event_times = event_timestamps_array{event};
        event_times = event_times(event_times > start_t & event_times < end_t);
        for event_time_i = 1:length(event_times)
            [~,align_ind] = min(abs(frame_grabs - event_times(event_time_i)));
            tmsp_ind_i = sortedIndices(:,1) == trial & sortedIndices(:,2) == align_ind;
            frame_by_event_matrix(tmsp_ind_i,event) = 1;
        end
    end
end

% get rid of trials with no neural activity
[indices2keep,~] = ismember(sortedIndices(:,1),trial_idx);
[indices2keep_forFs,~] = ismember(unique(sortedIndices(:,1)),trial_idx);

sortedIndices(~indices2keep,:) = [];
frame_by_event_matrix(~indices2keep,:) = [];
trialwise_Fs(~indices2keep_forFs) = [];

empiricalFs = mean(trialwise_Fs); % estimate the 'true' frame rate of the camera via the average of the trial-wise average imaging rates

