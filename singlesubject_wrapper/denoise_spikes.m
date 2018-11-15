function spikes_denoised = denoise_spikes(sortedIndices,fluorescence_array,point_process_array,dFF_thr,look_back_f,look_ahead_f,counts_flag,excludeIDs)
%% do event detection 
% idea of event-detection: look at identified spiking events (usually estimates from an algorithm like OASIS),
% then look at relative fluorescence values before and after spike event,
% then detect an event if either or both of the following is true:
% a) mean fluorescence after spike is at least dF/F0_thr, with F0 is
%   defined as lower 5th percentile of fluorescence values across [look_back] frames
%   before the spike, and dF is defined as average fluorescence value over
%   [look_ahead] frames after the spike, minus F0;
% b) average fluorescence value over [look_ahead] frames after the spike is greater
%   [global_thr], defined as the median + 3.5 * median-absolute-devian of
%   the neuron's entire fluorescence trace.
%  in this way, event detection leverages 1) spiking estimates from OASIS,
%  2) local DF/F0 gradients, and 3) global information about the neuron's
%  overall calcium activity

C = fluorescence_array;
all_spikes = point_process_array;

spikes_denoised = zeros(size(all_spikes));
[neuron_idx,spk_idx] = find(all_spikes > 0);
active_neurs = unique(neuron_idx);

for neuron = 1:length(active_neurs)
    
    if ismember(active_neurs(neuron),excludeIDs)
        fprintf('Skipping spike-denoising for neuron #: %d\n',active_neurs(neuron))
        spikes_denoised(active_neurs(neuron),:) = 0;
    else
        fprintf('Currently denoising spikes from neuron #: %d\n',active_neurs(neuron));
        
        neur_data = C(active_neurs(neuron),:);
        neur_data = neur_data - min(neur_data,[],2);
        
        spk_times = spk_idx(neuron_idx == active_neurs(neuron));
        
        ISIs = diff([0;spk_times;length(neur_data)]);
        
        global_thr = median(neur_data) + 3.5*mad(neur_data);
        
        if counts_flag
            good_spks_counts = zeros(length(spk_times),1);
        end
        
        good_spks = false(length(spk_times),1);
        
        for spk_i = 1:length(spk_times)
            spk = spk_times(spk_i);
            
            if and(ISIs(spk_i) > prctile(ISIs,20), ISIs(spk_i+1) > prctile(ISIs,20))
                good_spks(spk_i) = false;
            else
                
                curr_trial = sortedIndices(spk,1);
                trial_firstframe_ind = find(sortedIndices(:,1) == curr_trial & sortedIndices(:,2) == 1);
                trial_lastframe_ind = find(sortedIndices(:,1) == curr_trial & sortedIndices(:,2) == 799);
                look_back = min(look_back_f,spk-trial_firstframe_ind);
                look_ahead = min(look_ahead_f,trial_lastframe_ind-spk);
                
                if and(look_ahead >= 5, look_back >= 5)
                    
                    dat_segment = neur_data( (spk - look_back) : (spk + look_ahead) );
                    dat_segment = dat_segment - min(dat_segment); % subtracting the minimum, super important
                    
                    baseline = prctile(dat_segment(1:look_back),5);
                    change = mean(dat_segment( (end - look_ahead) : end));
                    
                    if or((change-baseline)/baseline >= dFF_thr, change >= global_thr)
                        good_spks(spk_i) = true;
                        if counts_flag
                            counts_ahead = min(5,trial_lastframe_ind-spk);
                            good_spks_counts(spk_i) = sum(spk_times >= spk & spk_times <= spk+counts_ahead);
                        end
                    end
                end
            end
        end
        
        if counts_flag
            spikes_denoised(active_neurs(neuron),spk_times(good_spks)) = good_spks_counts(good_spks);
        else
            spikes_denoised(active_neurs(neuron),spk_times(good_spks)) = 1;
        end
    end
end

end

