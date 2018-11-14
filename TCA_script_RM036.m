%% RM036 -- Fitting of tensor components analysis (via tensor decomposition) to the Neurons x TrialLength x NumTrials data

clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);

%% Initialize instance of SubjectObj for user-chosen rat and load data

rat_id = 21;
behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';
session_name = 'SA1';

currRat = SubjectObj;

currRat.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);

%% remove bad trials, if any

% currRat.remove_trials(40); % for instance, Rat21, Ext3 had a big artifact at Trial 40

%% Process fluorescence traces and estimate spiking (deconvolution + spike-thresholding, etc.)

% smooth data with LLR

win_len = 2; win_overlap = ceil(win_len/2);
currRat.smooth_C(win_len,win_overlap);

% deconvolve data

currRat.deconvTrials()

% denoise spikes estimated via deconvolution with a local/global
% thresholding process

dFF_thr = 6;
look_back_f = 10;
look_ahead_f = 50;
counts_flag = true;
currRat.denoiseSpikes(currRat.C,threshold_timeseries(currRat.spikes,'mean',1.5),dFF_thr,look_back_f,look_ahead_f,counts_flag)

% visualize results of spike-denoising

currRat.overlay_spikes_calc(currRat.C,double(currRat.spikes_denoised > 0))

% convolve spikes

spike_data = sqrt(currRat.spikes_denoised); % square root transform before convolution
transpose_flag = 0; 
kernel = normpdf(-3:3,0,1.5);
currRat.convSpikes(spike_data,transpose_flag,kernel)


%% tensor components analysis (TCA) via canonical polyadic tensor decomposition

trial_indices = [10:780];

data2use = currRat.spikes_conv;
% normalized = bsxfun(@minus,data2use,mean(data2use,2));
% normalized = bsxfun(@rdivide,data2use - min(data2use,[],2),(max(data2use,[],2) - min(data2use,[],2)));
normalized = data2use;
trial_reshaped = reshape(normalized,size(normalized,1),799,currRat.num_trials);
trial_reshaped = trial_reshaped(:,trial_indices,:);
clear normalized data2use;

data = tensor(trial_reshaped);

n_init = 20;
r_vector = 1:10;
errors = zeros(n_init,length(r_vector));
similarity_scores = zeros(n_init,length(r_vector));

single_neuron_rsquared = zeros(n_init,size(data,1),length(r_vector));

for ii = 1:length(r_vector)
    tic
    
    r = r_vector(ii);

    fprintf('Running models with rank %d\n',r)
    
    all_models = cell(1,n_init);
    
    for jj = 1:n_init
%         est_factors = cp_als(data,r);
        all_models{jj} = ncp(data,r);
        errors(jj,ii) = norm(full(all_models{jj}) - data)/norm(data);
    end
    
    [~,best_model_idx] = min(errors(:,ii));
    
    best_model = all_models{best_model_idx};
    for jj = 1:n_init
        compare_model = all_models{jj};
        [~,similarity_scores(jj,ii)] = Hungarian(pdist2(best_model.U{1}',compare_model.U{1}','cosine'));
        single_neuron_rsquared(jj,:,ii) = compute_rsquared_dist(all_models{jj},data);
    end

%     for jj = 1:n_init
%         if r > 8
%             similarity_scores(jj,ii) = score(all_models{best_model_idx},all_models{jj},'greedy',true);
%         else
%             similarity_scores(jj,ii) = score(all_models{best_model_idx},all_models{jj},'greedy',false);
%         end
%     end
    
    clear all_models;
    
    fprintf('Time taken to run all rank %d models: %.2f seconds\n',r,toc)
    
end

figure(1);
plot(r_vector,mean(errors,1));
hold on; plot(r_vector,mean(errors,1) + iqr(errors),'r--'); 
plot(r_vector,mean(errors,1) - iqr(errors),'r--');


figure(2);
scatter(reshape(repmat(r_vector,n_init,1),n_init*length(r_vector),1),reshape(similarity_scores,n_init*length(r_vector),1))

figure(3);
r_squared_neuronavg = squeeze(mean(single_neuron_rsquared,2));
plot(r_vector,mean(r_squared_neuronavg,1));
hold on; plot(r_vector,mean(r_squared_neuronavg,1) + iqr(r_squared_neuronavg),'r--');
plot(r_vector,median(r_squared_neuronavg,1) - iqr(r_squared_neuronavg),'r--');

%% settle on some number of factors

est_factors = ncp(data,5);
assembly_info = detect_ensembles(est_factors,trial_reshaped,5,5);

%% relate trial-to-trial variability of neuronal responses to ensemble-membership

[entropies, spike_time_var ] = compute_var(trial_reshaped,2);

[membership_matrix] = make_membership_matrix(assembly_info);

scatter(entropies(sum(membership_matrix>0,2) == 0),[assembly_info(sum(membership_matrix>0,2) == 0).r_sq],'co')
hold on;
scatter(entropies(sum(membership_matrix>0,2) == 1),[assembly_info(sum(membership_matrix>0,2) == 1).r_sq],'ro')
hold on;
scatter(entropies(sum(membership_matrix>0,2) == 2),[assembly_info(sum(membership_matrix>0,2) == 2).r_sq],'bo')
hold on;
scatter(entropies(sum(membership_matrix>0,2) == 3),[assembly_info(sum(membership_matrix>0,2) == 3).r_sq],'ko')

ensemble_membership_entropy(trial_reshaped,2,assembly_info)


%% ways to visualize ensembles

trial_indices = 10:780;

events_reshaped = reshape(full(currRat.event_matrix)',size(currRat.event_matrix,2),799,100);
which_trials = [1:100];
events_reshaped = events_reshaped(:,trial_indices,which_trials);
events2show = {'HLON','press','cueOFF','bbk','pellets','HLOFF'}; % key events to display (don't show 'lever IN/OUT' or 'cueON' due to overlap with press

plot_ensemble_activity_trialavg(est_factors,currRat.spikes_conv,799,100,trial_indices,assembly_info,currRat.event_names,events_reshaped,events2show)

events_reshaped = reshape(full(currRat.event_matrix)',size(currRat.event_matrix,2),799,100);
which_trials = [1:100];
events_reshaped = events_reshaped(:,trial_indices,which_trials);
which_events = [1,4,5,6,7,8,9];

plot_ensemble_activity_trialwise(est_factors,trial_reshaped,events_reshaped,currRat.event_names,which_events,which_trials,assembly_info)

%% more ways to visualize assemblies

all_ensembles_raster(currRat.C_decon,assembly_info)

[within_ens_corr,between_ens_corr] = within_between_ens_correlations(currRat.spikes_conv,assembly_info,false);

stats_results =  compare_correlations(within_ens_corr,between_ens_corr,0.01);

viz_ktensor(est_factors, ... 
            'Plottype', {'bar', 'line', 'scatter'}, ...
            'Modetitles', {'neurons', 'time', 'trials'})









