function [ population_vectors ] = populationVectors_accum(across_sess_array,sess_idx,cell_map,event_nam,surround_idx,type)
% populationCorr_overTime Look at correlation between population vectors on
% specific days, aligned to specific events, in specific trial types (i.e. failed or
% successful)

across_sess_array = across_sess_array(:,sess_idx);
cell_map = cell_map(:,sess_idx);

before_align = surround_idx(1); after_align = surround_idx(2);

common_neurons = find(sum(cell_map > 0,2) == length(sess_idx));
num_common = length(common_neurons);

population_vectors = [];

for i = 1:length(sess_idx)
    
    sess_nam = sess_idx(i);
    this_sess_data = across_sess_array(:,i);
    event_names = this_sess_data{2}{1};
    
    [~,T] = size(this_sess_data{1});
    num_trials = T/799;
    
    common_neuron_tempidx = cell_map(common_neurons,i);

    event_idx = find(cellfun(@(x) ~isempty(x),strfind(event_names,event_nam)));
    
    press_idx = find(cellfun(@(x) ~isempty(x),strfind(event_names,'press'))); % use to find out whether trial is a failed/successful trial
    
    reshaped_events = permute(reshape(full(this_sess_data{2}{2}),[],num_trials,length(event_names)),[1,3,2]);
    failed_trial_idx = sum(squeeze(reshaped_events(:,press_idx,:)),1) == 0;
    
    reshaped_data = this_sess_data{1};
    reshaped_data = reshape(reshaped_data(common_neuron_tempidx,:), num_common, 799, num_trials);
    
    for trial = 1:num_trials
        currTrialneur = squeeze(reshaped_data(:,:,trial));
        event_tmsp_trial = find(squeeze(reshaped_events(:,event_idx,trial)));
        if failed_trial_idx(trial) == 1
            if strcmp(type,'average')
                population_vectors = [population_vectors,...
                    [sess_nam;0;mean(currTrialneur(:,(event_tmsp_trial-before_align):(event_tmsp_trial+after_align)),2)]];
            else
                chunk2add = [sess_nam*ones(1,before_align+after_align+1);zeros(1,before_align+after_align+1);...
                    currTrialneur(:,(event_tmsp_trial-before_align):(event_tmsp_trial+after_align))];
                population_vectors = [population_vectors,chunk2add];
            end
        else
            if strcmp(type,'average')
                population_vectors = [population_vectors,...
                    [sess_nam;1;mean(currTrialneur(:,(event_tmsp_trial-before_align):(event_tmsp_trial+after_align)),2)]];
            else
                
                chunk2add = [sess_nam*ones(1,before_align+after_align+1);ones(1,before_align+after_align+1);...
                    currTrialneur(:,(event_tmsp_trial-before_align):(event_tmsp_trial+after_align))];
                population_vectors = [population_vectors,chunk2add];
            end
        end
    end
    
end
        
        
        
        
        
    
  
    
    
    
    
    

    
    




end

