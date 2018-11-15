function [ stats_results ] = test_firingrates(test_type,event_locked,surround_time,event_tmsp)
% Non-parametric tests of pre- and post-event distributions, aligned to
% behavioral event indexed relative to segments in event_locked, index
% given in event_tmsp. Size of event-relative periods determined by
% event_tmsp - surround_time(1) and event_tmsp + surround_time(2)
%   test_type options: 'wilcoxon' and 'ks' (Kolgomorov-Smirnov)

num_events = length(event_locked);
num_neurons = size(event_locked{1},1);

stats_results.SigMatrix = zeros(num_neurons,num_events);
stats_results.pVals = zeros(num_neurons,num_events);
stats_results.Modulation_direction = zeros(num_neurons,num_events);
stats_results.StatTest = test_type;

for event_id = 1:num_events
    
    event_id_data = event_locked{event_id};
    
    look_behind = surround_time(1);
    look_ahead = surround_time(2);
    
    start_idx = max(1,event_tmsp-look_behind);
    end_idx = min(size(event_id_data,2),event_tmsp+look_ahead);
    
    pre_event_avg = squeeze(mean(event_id_data(:,start_idx:(event_tmsp-1),:),3));
    post_event_avg = squeeze(mean(event_id_data(:,event_tmsp:end_idx,:),3));
    
    pre_event_avg(isnan(pre_event_avg)) = 0;
    post_event_avg(isnan(post_event_avg)) = 0;
    
    for neuron = 1:size(event_id_data,1)
        switch lower(test_type)
            case 'ks'
                [sig_2sided,pVal_2sided] = kstest2(pre_event_avg(neuron,:),post_event_avg(neuron,:)); %tests firstly whether the two distributions are different, irrespective of directionality
                stats_results.SigMatrix(neuron,event_id) = sig_2sided;
                if sig_2sided
                    [sig_larger,pVal_larger] = kstest2(pre_event_avg(neuron,:),post_event_avg(neuron,:),'tail','larger'); % tests whether second argument/distribution is larger than the first
                    if ~sig_larger
                        [sig_smaller,pVal_smaller] = kstest2(pre_event_avg(neuron,:),post_event_avg(neuron,:),'tail','smaller'); % tests whether second argument/distribution is smaller than the first
                        if sig_smaller
                            stats_results.pVals(neuron,event_id) = pVal_smaller;
                            stats_results.Modulation_direction(neuron,event_id) = -1;
                        elseif ~sig_smaller
                            stats_results.pVals(neuron,event_id) = pVal_2sided;
                            stats_results.Modulation_direction(neuron,event_id) = NaN;
                        end
                    elseif sig_larger
                        stats_results.pVals(neuron,event_id) = pVal_larger;
                        stats_results.Modulation_direction(neuron,event_id) = 1;
                    end
                elseif ~sig_2sided
                    stats_results.pVals(neuron,event_id) = NaN;
                    stats_results.Modulation_direction(neuron,event_id) = NaN;
                end
            case 'wilcoxon'
                [pVal_2sided,sig_2sided] = ranksum(pre_event_avg(neuron,:),post_event_avg(neuron,:));
                stats_results.SigMatrix(neuron,event_id) = sig_2sided;
                if sig_2sided
                    [pVal_larger,sig_larger] = ranksum(pre_event_avg(neuron,:),post_event_avg(neuron,:),'tail','left'); % tests whether second argument/distribution is larger than the first
                    if ~sig_larger
                        [pVal_smaller,sig_smaller] = ranksum(pre_event_avg(neuron,:),post_event_avg(neuron,:),'tail','right'); % tests whether second argument/distribution is smaller than the first
                        if sig_smaller
                            stats_results.pVals(neuron,event_id) = pVal_smaller;
                            stats_results.Modulation_direction(neuron,event_id) = -1;
                        elseif ~sig_smaller
                            stats_results.pVals(neuron,event_id) = pVal_2sided;
                            stats_results.Modulation_direction(neuron,event_id) = NaN;
                        end
                    elseif sig_larger
                        stats_results.pVals(neuron,event_id) = pVal_larger;
                        stats_results.Modulation_direction(neuron,event_id) = 1;
                    end
                elseif ~sig_2sided
                    stats_results.pVals(neuron,event_id) = NaN;
                    stats_results.Modulation_direction(neuron,event_id) = NaN;
                end
        end
    end
    
end

