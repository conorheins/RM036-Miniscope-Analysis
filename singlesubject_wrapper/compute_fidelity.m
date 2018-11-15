function [ fidelity_matrix ] = compute_fidelity(event_locked,event_tmsp,surround_time,stats_results)
% Compute response fidelity of neurons that have been deemed 'significant responders' by statistical test 
%   Detailed explanation goes here


fidelity_matrix = nan(size(stats_results.SigMatrix));

for i = 1:length(event_locked)
    
    pos_sig_idx = stats_results.Modulation_direction(:,i) == 1;
    neg_sig_idx = stats_results.Modulation_direction(:,i) == -1;
    
    if size(event_locked{i},3) >= 5
        
        % first do positively-modulated neurons
        traces = event_locked{i}(pos_sig_idx,:,:);
        
        num_trials = size(traces,3);
        
        start_idx = max(1,event_tmsp - surround_time(1));
        end_idx = min(size(traces,2),event_tmsp + surround_time(2));
        
        pre_event = squeeze(mean(traces(:,start_idx:(event_tmsp-1),:),2));
        post_event = squeeze(mean(traces(:,event_tmsp:end_idx,:),2));
        
        if length(find(pos_sig_idx)) == 1
            response_fidelity = sum((post_event > 3 .* pre_event),1)./num_trials;
        else
            response_fidelity = sum((post_event > 3 .* pre_event),2)./num_trials;
        end
        
        fidelity_matrix(pos_sig_idx,i) = response_fidelity;

        
        % now do negatively-modulated neurons
        traces = event_locked{i}(neg_sig_idx,:,:);
        
        num_trials = size(traces,3);
        
        start_idx = max(1,event_tmsp - surround_time(1));
        end_idx = min(size(traces,2),event_tmsp + surround_time(2));
        
        pre_event = squeeze(mean(traces(:,start_idx:(event_tmsp-1),:),2));
        post_event = squeeze(mean(traces(:,event_tmsp:end_idx,:),2));
        
        if length(find(neg_sig_idx)) == 1
            response_fidelity = sum((pre_event > 3 .* post_event),1)./num_trials;
        else
            response_fidelity = sum((pre_event > 3 .* post_event),2)./num_trials;
        end
       
        fidelity_matrix(neg_sig_idx,i) = response_fidelity;
        
    end

end

