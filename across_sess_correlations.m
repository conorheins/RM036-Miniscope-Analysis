function corrMat = across_sess_correlations(across_sess_array,sess_idx,cell_map,event_nam,surround_idx,type,display_flag,session_names)
% across_sess_correlations Generates correlation image of population
% activity timed to desired experimental events over sessions -- also
% splits into successful vs. failed trials
%   Detailed explanation goes here

across_sess_array = across_sess_array(:,sess_idx);
cell_map = cell_map(:,sess_idx);
session_names = session_names(sess_idx);

before_align = surround_idx(1); after_align = surround_idx(2);

common_neurons = find(sum(cell_map > 0,2) == length(sess_idx));
num_common = length(common_neurons);

if strcmp(type,'average')
    population_vectors = zeros(length(sess_idx),num_common,2);
elseif strcmp(type,'spike_time')
    population_vectors = zeros(length(sess_idx),num_common,before_align+after_align+1,2);
end

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
    
    
    if strcmp(type,'average')
        
        failed_avg = zeros(num_common,1);
        succ_avg = zeros(num_common,1);
        
        [when_in_trial,which_trial] = find(squeeze(reshaped_events(:,event_idx,:)));
        
        for j = 1:length(when_in_trial)
            if failed_trial_idx(which_trial(j)) == 1
                failed_avg = failed_avg + mean(reshaped_data(:,(when_in_trial(j) - before_align):(when_in_trial(j) + after_align),which_trial(j)),2);
            elseif failed_trial_idx(which_trial(j)) == 0
                succ_avg = succ_avg + mean(reshaped_data(:,(when_in_trial(j) - before_align):(when_in_trial(j) + after_align),which_trial(j)),2);
            end
        end
        
        population_vectors(i,:,1) = failed_avg./length(when_in_trial);
        population_vectors(i,:,2) = succ_avg./length(when_in_trial);
        
    elseif strcmp(type,'spike_time')
        
        failed_arrayz = zeros(num_common,before_align+after_align+1);
        succ_arrayz = zeros(num_common,before_align+after_align+1);
        
        [when_in_trial,which_trial] = find(squeeze(reshaped_events(:,event_idx,:)));
        
        for j = 1:length(when_in_trial)
            if failed_trial_idx(which_trial(j)) == 1
                failed_arrayz = failed_arrayz + reshaped_data(:,(when_in_trial(j) - before_align):(when_in_trial(j) + after_align),which_trial(j));
            elseif failed_trial_idx(which_trial(j)) == 0
                succ_arrayz = succ_arrayz + reshaped_data(:,(when_in_trial(j) - before_align):(when_in_trial(j) + after_align),which_trial(j));
            end
        
        
        end
        
        population_vectors(i,:,:,1) = failed_arrayz./length(when_in_trial);
        population_vectors(i,:,:,2) = succ_arrayz./length(when_in_trial);
        
    end
    
        
end

corrMat = zeros(length(sess_idx),length(sess_idx),2);

if strcmp(type,'average')
    
    corrMat(:,:,1) = corrcoef(population_vectors(:,:,1)');
    corrMat(:,:,2) = corrcoef(population_vectors(:,:,2)');
    
elseif strcmp(type,'spike_time')
    
    for i = 1:length(sess_idx)
        
        pattern_i_f = squeeze(population_vectors(i,:,:,1));
        pattern_i_s = squeeze(population_vectors(i,:,:,2));
        
        for j = 1:length(sess_idx)
            
            
            pattern_j_f = squeeze(population_vectors(j,:,:,1));
            pattern_j_s = squeeze(population_vectors(j,:,:,2));
            
            corrMat(i,j,1) = corr(pattern_i_f(:),pattern_j_f(:));
            corrMat(i,j,2) = corr(pattern_i_s(:),pattern_j_s(:));
            
        end
        
    end
    
end
            
            
            
if display_flag
    
    cmin = 0;
   
%     cmax = max([max(get_lower_tri(corrMat(:,:,1))),max(get_lower_tri(corrMat(:,:,2))]);
    cmax = 1;
    
%     for_display_succ = corrMat(:,:,2) - eye(length(sess_idx));
    for_display_succ = corrMat(:,:,2);
%     for_display_failed = corrMat(:,:,1) - eye(length(sess_idx));
    for_display_failed = corrMat(:,:,1);
    
    
    figure;
    
    subplot(121);
    imagesc(for_display_succ); caxis([cmin, cmax]); colormap gray; colorbar;
    yticks(1:length(sess_idx));
    yticklabels(session_names);
    
    xticks(1:length(sess_idx));
    xticklabels(session_names);

    title('Ensemble similarity between subsequent days (Successful Trials)')

    
    subplot(122);
    imagesc(for_display_failed); caxis([cmin, cmax]); colormap gray; colorbar;
    
    yticks(1:length(sess_idx));
    yticklabels(session_names);
    
    xticks(1:length(sess_idx));
    xticklabels(session_names);

    title('Ensemble activation similarity between subsequent days (Failed Trials)')
    
end

    

end
        

