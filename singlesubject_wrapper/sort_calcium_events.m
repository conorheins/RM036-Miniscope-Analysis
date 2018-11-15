function [sorted_spikes] = sort_calcium_events(spk_array,fluorescence_array,options)
    
    [all_neuron_idx, all_spk_times] = find(spk_array);
    
    unit_nams = unique(all_neuron_idx);
    
    ISI_thresh = options.ISI_thresh;
    before_win = options.before_win;
    after_win = options.after_win;
    K = options.K;
    
    sorted_spikes = zeros(size(spk_array));
    
    for ii = 1:length(unit_nams)
        
        all_spks = all_spk_times(all_neuron_idx==unit_nams(ii));
        
        spk_2_test = [all_spks(1); all_spks(find(diff(all_spks) > ISI_thresh) + 1)];
        
        features = zeros(length(spk_2_test),2);
        
        transients = zeros( length(spk_2_test), (before_win + after_win + 1) );
        edge_jj = [];
        
        for jj = 1:length(spk_2_test)
            
            if or( (spk_2_test(jj) - before_win) <= 0, (spk_2_test(jj) + after_win) > size(fluorescence_array,2))
                edge_jj = [edge_jj,jj];
            else
                transients(jj,:) = fluorescence_array(unit_nams(ii),(spk_2_test(jj)-before_win):(spk_2_test(jj)+after_win));
            end

        end
        
        transients(edge_jj,:) = [];
        spk_2_test(edge_jj) = [];
        features(edge_jj,:) = [];
        
        transients = transients - min(transients,[],2);
        transients = bsxfun(@rdivide,transients,max(transients,[],2));
        
        features(:,1) = prctile(transients(:,(before_win +1):end),80,2)./prctile(transients(:,(1:before_win)),5,2);
        features(:,2) = mean(transients(:,(before_win +1)),2)./mean(transients(:,(1:before_win)),2);
        
        [labels,centroidz] = kmeans(features,K);
               
        clust_scores = zeros(K,1);
        for clust = 1:K
            clust_scores(clust) = sum(sum(bsxfun(@gt,centroidz(clust,:),centroidz),1));
        end
            
        [~,true_clust] = max(clust_scores);
            
        good_spikes = spk_2_test(labels==true_clust);
        
        sorted_spikes(unit_nams(ii),good_spikes) = 1;
     
    end
        
end