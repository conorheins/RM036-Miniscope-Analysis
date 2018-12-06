function [ stats_results] = kstest_change_scores(change_scores,event_names,num_neurons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

stats_results = struct('SigMatrix',[]);
stats_results = setfield(stats_results,'SigMatrix',zeros(num_neurons,length(event_names)));
stats_results = setfield(stats_results,'pVals',zeros(num_neurons,length(event_names)));
stats_results = setfield(stats_results,'Modulation_direction',zeros(num_neurons,length(event_names)));
stats_results = setfield(stats_results,'StatTest','Kolmogorov-Smirnov');


for event_id = 1:length(event_names)
    
    change_score_real = change_scores{event_id}{1};
    change_score_rand = change_scores{event_id}{2};
    
    for neuron = 1:num_neurons
        
        [stats_results.SigMatrix(neuron,event_id),stats_results.pVals(neuron,event_id)] = kstest2(change_score_real(neuron,:),change_score_rand(neuron,:));
        
        if stats_results.SigMatrix(neuron,event_id) == 1
            [stats_results.SigMatrix(neuron,event_id),stats_results.pVals(neuron,event_id)] = kstest2(change_score_real(neuron,:),change_score_rand(neuron,:),'tail','larger');
            if stats_results.SigMatrix(neuron,event_id) == 0
                [stats_results.SigMatrix(neuron,event_id),stats_results.pVals(neuron,event_id)] = kstest2(change_score_real(neuron,:),change_score_rand(neuron,:),'tail','smaller');
                if stats_results.SigMatrix(neuron,event_id) == 1
                   stats_results.Modulation_direction(neuron,event_id) = 1;
                else
                    stats_results.Modulation_direction(neuron,event_id) = 0;
                end
            else
                stats_results.Modulation_direction(neuron,event_id) = -1;
            end
        else
            stats_results.Modulation_direction(neuron,event_id) = 0;
        end       
    end
    
end




end

