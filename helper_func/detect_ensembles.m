function [ assembly_info ] = detect_ensembles(est_factors,trial_reshaped,num_components,num_factor_combos)
% visualize_ensembles 
% Function to detect cell ensemble activity identified by non-negative TCA
%   Searches combinations of factors to best explain single-neuron
%   activity. Those factors that contribute to single-neuron activity
%   are deemed the 'ensembles' that the cell is a participant of

[nNeurons,nT,nTrials] = size(trial_reshaped);

cell_assemblies = est_factors.U{1}(:,1:num_components);
trial_dyn = est_factors.U{2}(:,1:num_components);
session_dyn = est_factors.U{3}(:,1:num_components);

full_data = reshape(trial_reshaped,nNeurons,nT*nTrials);
          
factor_combos = cell(1,num_factor_combos);
factor_combo_r_sq = cell(1,num_factor_combos);

for num_fac = 1:num_factor_combos
    
    ensemble_combos = nchoosek(1:num_components,num_fac); 
    
    factor_combos{num_fac} = ensemble_combos;
    
    factor_combo_r_sq{num_fac} = zeros(nNeurons,size(ensemble_combos,1));
    
end
    
   
for num_fac = 1:num_factor_combos
    
    ensemble_combos = factor_combos{num_fac};
    
    for jj = 1:size(ensemble_combos,1)
        
        data_recon = zeros(nNeurons,nT,nTrials);

        factors2use = ensemble_combos(jj,:);
  
        % reconstruct the data with tensor product
        for t = 1:nTrials           
            for kk = 1:length(factors2use)
                data_recon(:,:,t) = data_recon(:,:,t) + cell_assemblies(:,factors2use(kk))*trial_dyn(:,factors2use(kk))' * session_dyn(t,factors2use(kk))';
            end              
        end
        
        data_recon = reshape(data_recon,nNeurons,nT*nTrials);
        
        units2check = find(sum(cell_assemblies(:,factors2use),2));
        
        for neur = 1:length(units2check)           
            unit_ID = units2check(neur);
            factor_combo_r_sq{num_fac}(unit_ID,jj) = corr(full_data(unit_ID,:)',data_recon(unit_ID,:)').^2;
        end
          
    end
end


max_r_sq_matrix = zeros(nNeurons,num_factor_combos);
best_combo_idx = zeros(nNeurons,num_factor_combos);
for neur = 1:nNeurons
    for num_fac = 1:num_factor_combos
        [max_r_sq_matrix(neur,num_fac),best_combo_idx(neur,num_fac)] = max(factor_combo_r_sq{num_fac}(neur,:));
    end
end

ensemble_neurons = find(max(max_r_sq_matrix,[],2) > 0.01);

assembly_info = [];

for neur = 1:nNeurons
    if ~ismember(neur,ensemble_neurons)
        assembly_info(neur).membership = [];
        assembly_info(neur).r_sq = max(max_r_sq_matrix(neur,:));
    else
        dR = diff(max_r_sq_matrix(neur,:),1);
        if ~isempty(find(max_r_sq_matrix(neur,:)./max_r_sq_matrix(neur,1) > 1.25, 1))
            rise_idx = find(dR == max(dR));
            after_rise_idx = find(dR == max(dR)) + 1;
            if dR(after_rise_idx) > 0.5 * dR(rise_idx)
                factor_idx = after_rise_idx + 1;
            else
                factor_idx = after_rise_idx;
            end
            assembly_info(neur).membership = factor_combos{factor_idx}(best_combo_idx(neur,factor_idx),:);
            assembly_info(neur).r_sq = max_r_sq_matrix(neur,factor_idx);
        else
            assembly_info(neur).membership = factor_combos{1}(best_combo_idx(neur,1));
            assembly_info(neur).r_sq = max_r_sq_matrix(neur,1);
        end
    end    
end

            
end

