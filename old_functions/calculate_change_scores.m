function [change_scores] = calculate_change_scores(sortedIndices,event_names,event_matrix,data2use,surround_time)

% calculate_change_scores calculates cell array (num_events,1) of all the
% change-score distributions for each event -- each cell_array is a 2 x 1
% cell array with the change_score_real (the actual difference in firing
% rate induced by the event) and change_score_rand (randomly sampled
% timestamp as a null distribution).
%   Calculate as: (FRevent - FRbaseline)/(FRevent + FRbaseline)
%   compare against: (FRrandomsample - FRbaseline)/(FRrandomsample +
%   FRbaseline)

num_events = length(event_names);
change_scores = cell(num_events,1);
num_trials = max(unique(sortedIndices(:,1)));

preframes = surround_time(1);
postframes = surround_time(2);

for event_i = 1:num_events
    
    event_times = find(event_matrix(:,event_i));
    
    between_events = diff(event_times);
    throw_outs = find(between_events <= min(surround_time));
    event_times(throw_outs + 1) = [];
    
    change_score_real = [];
    change_score_rand = [];
    
    for tmsp_i = 1:length(event_times)
        
        timestamp = event_times(tmsp_i);
        
        curr_trial = sortedIndices(timestamp,1);
        trial_firstframe_ind = find(sortedIndices(:,1) == curr_trial & sortedIndices(:,2) == 1);
        trial_lastframe_ind = find(sortedIndices(:,1) == curr_trial & sortedIndices(:,2) == 799);
        look_back = min(preframes,timestamp-trial_firstframe_ind);
        look_ahead = min(postframes,trial_lastframe_ind-timestamp);
        
        event_locked = mean(data2use(:,timestamp:(timestamp+look_ahead)),2);
        baseline = mean(data2use(:,(timestamp-look_back):(timestamp-1)),2);
        
        change_score_real = [change_score_real, (event_locked - baseline)./(event_locked + baseline)];
        
        rand_trial_id = randperm(num_trials,1);
        rand_trial = data2use(:,sortedIndices(:,1)==rand_trial_id);
        rand_align_ind = randi([(preframes + 1),(size(rand_trial,2) - postframes)],1);
        
        rand_locked = mean(rand_trial(:,rand_align_ind:(rand_align_ind + surround_time)),2);
        rand_baseline = mean(rand_trial(:,(rand_align_ind - surround_time):(rand_align_ind-1)),2);
        change_score_rand = [change_score_rand, (rand_locked - rand_baseline)./(rand_locked + rand_baseline)];
        
    end
    
    change_score_real(isnan(change_score_real)) = 0;
    change_score_rand(isnan(change_score_rand)) = 0;
    
    change_scores{event_i} = {change_score_real,change_score_rand};
end

