function [stats_results] = permutation_test(event_locked,surround_time,event_tmsp,num_shuffles,alpha)

% PERMUTATION_TEST -- Uses permutation_testing to see if
% mean of two periods surrounding event-period are significantly
% different from eachother. Then we look at which mean (pre- or post-event) is higher and then
% assigns directionality of modulation (+ or -) based on which is higher

num_events = length(event_locked);
num_neurons = size(event_locked{1},1);

stats_results.SigMatrix = zeros(num_neurons,num_events);
stats_results.pVals = zeros(num_neurons,num_events);
stats_results.Modulation_direction = zeros(num_neurons,num_events);
stats_results.StatTest = 'permutation';

for event_id = 1:num_events
    
    event_id_data = event_locked{event_id};
    
    if size(event_id_data,3) >= 5
    
        look_behind = surround_time(1);
        look_ahead = surround_time(2);
        
        start_idx = max(1,event_tmsp-look_behind);
        end_idx = min(size(event_id_data,2),event_tmsp+look_ahead);
        
        pre_event = squeeze(mean(event_id_data(:,start_idx:(event_tmsp-1),:),2));
        post_event = squeeze(mean(event_id_data(:,event_tmsp:end_idx,:),2));
        
        pre_event(isnan(pre_event)) = 0;
        post_event(isnan(post_event)) = 0;
        
        labels = [zeros(1,size(pre_event,2)),ones(1,size(post_event,2))];
        
        all_dat = [labels;[pre_event,post_event]];
        test_statistic = mean(all_dat(2:end,labels==1),2) - mean(all_dat(2:end,labels==0),2);
        %     test_statistic = all_dat(2:end,labels==1)./all_dat(2:end,labels==0);
        test_statistic(isnan(test_statistic)) = 0;
        
        null_dist = zeros(num_neurons,num_shuffles);
        
        for shuff = 1:num_shuffles
            shuffled = [labels;all_dat(2:end,randperm(length(labels)))]; % this is where the shuffling happens
            %         null_dist(:,shuff) = shuffled(2:end,labels==1)./shuffled(2:end,labels==0);  % other option for test statistic, ratio
            null_dist(:,shuff) = mean(shuffled(2:end,labels==1),2) - mean(shuffled(2:end,labels==0),2); % test statistic of choice, difference of means
        end
        
        for neuron = 1:num_neurons
            if test_statistic(neuron) < mean(null_dist(neuron,:))
                stats_results.pVals(neuron,event_id) = length(find(null_dist(neuron,:) < test_statistic(neuron)))/num_shuffles;
                if stats_results.pVals(neuron,event_id) < alpha
                    stats_results.SigMatrix(neuron,event_id) = 1;
                    stats_results.Modulation_direction(neuron,event_id) = -1;
                else
                    stats_results.pVals(neuron,event_id) = NaN;
                    stats_results.Modulation_direction(neuron,event_id) = NaN;
                end
            elseif test_statistic(neuron) > mean(null_dist(neuron,:))
                stats_results.pVals(neuron,event_id) = length(find(null_dist(neuron,:) > test_statistic(neuron)))/num_shuffles;
                if stats_results.pVals(neuron,event_id) < alpha
                    stats_results.SigMatrix(neuron,event_id) = 1;
                    stats_results.Modulation_direction(neuron,event_id) = 1;
                else
                    stats_results.pVals(neuron,event_id) = NaN;
                    stats_results.Modulation_direction(neuron,event_id) = NaN;
                end
            end
        end
    else
        stats_results.pVals(:,event_id) = nan(num_neurons,1);
        stats_results.SigMatrix(:,event_id) = nan(num_neurons,1);
        stats_results.Modulation_direction(:,event_id) = nan(num_neurons,1);
    end
            
end
