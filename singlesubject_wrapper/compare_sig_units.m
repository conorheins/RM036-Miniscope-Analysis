function [ distance_matrix,distance_matrix_shuffled ] = compare_sig_units(session_objects,event_type)
%compare_sig_units Look at how response selectivity persists between days
%   Use the hamming distance between units' selectivity vectors

num_neurons = size(session_objects{1}.C,1);
num_sessions = length(session_objects);

if ~exist('event_type','var') || isempty(event_type)
    event_type = session_objects{end}.event_names;
end

if iscell(event_type)
    event_idx = zeros(length(event_type),1);
    for i = 1:length(event_type)
        event_idx(i) = find(cellfun(@(x) ~isempty(x), strfind(session_objects{end}.event_names,event_type{i})));
    end
else
    event_idx = find(cellfun(@(x) ~isempty(x), strfind(session_objects{end}.event_names,event_type)));
end


distance_matrix = zeros(num_sessions,num_sessions,num_neurons);
distance_matrix_shuffled = zeros(num_sessions,num_sessions,num_neurons);

for i = 1:length(session_objects)
    
    sigMatrix_i = session_objects{i}.stats_results.SigMatrix(:,event_idx);
    
    for j = 1:length(session_objects)
        
        sigMatrix_j = session_objects{j}.stats_results.SigMatrix(:,event_idx);
        
        shuffled_j = sigMatrix_j;
        
        for n = 1:num_neurons
            shuffled_j(n,:) = sigMatrix_j(n,randperm(length(event_idx)));
        end
        
        distances = pdist2(sigMatrix_i,sigMatrix_j,'hamming');
        distances_shuff = pdist2(sigMatrix_i,shuffled_j,'hamming');
        
        distance_matrix(i,j,:) = reshape(diag(distances),[1,1,num_neurons]);
        distance_matrix_shuffled(i,j,:) = reshape(diag(distances_shuff),[1,1,num_neurons]);
        
    end
    
end
        
        
        
        



        

