%% set up paths and stuff
clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));
addpath(genpath('helper_func'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);

%% Initialize instance of SubjectObj for user-chosen rat and load data
tic
rat_id = 8;
behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';

% define a bunch of parameters used in the time-series processing
params.smooth_flag = false; % whether to smooth data
params.win_len = 2; % if smoothing, the length in frames of the window over which to smooth (to compute a locally-linear regression)
params.win_overlap = ceil(params.win_len/2); % overlap (in frames) of the windows
params.threshold_params.stat_type = 'mean'; % before spike denoising, parameter for how to threshold the deconvolved spikes array
params.threshold_params.thr_factor = 1.5; % before spike denoising, parameter for how to threshold the deconvolved spikes array
params.denoising_params.dFF_thr = 6; % spike denoising param controlling factor of dFF-increase, post-spike
params.denoising_params.look_back_f = 10; % spike denoising param controlling how much history to take into account when computing pre-spike baseline
params.denoising_params.look_ahead_f = 50; % spike denoising param controlling how much future to take into account when estimating post-spike fluorescence change
params.denoising_params.counts_flag = true; % whether to estimate spike-counts in final denoised_spk output
params.conv_kernel = normpdf(-3:3,0,1.5); % kernel used to convolve denoised-spike traces
params.quick_view = true; %flag to quickly view entire population activity, in case there's some whack trial that needs cutting out

SA1 = process_single_rat(rat_id,'SA1',behav_data_folder,neural_data_folder,timestamps_folder,params);
SA2 = process_single_rat(rat_id,'SA2',behav_data_folder,neural_data_folder,timestamps_folder,params);
Ext1 = process_single_rat(rat_id,'Ext1',behav_data_folder,neural_data_folder,timestamps_folder,params);
Ext2 = process_single_rat(rat_id,'Ext2',behav_data_folder,neural_data_folder,timestamps_folder,params);
Ext3 = process_single_rat(rat_id,'Ext3',behav_data_folder,neural_data_folder,timestamps_folder,params);
Ext4 = process_single_rat(rat_id,'Ext4',behav_data_folder,neural_data_folder,timestamps_folder,params);
Reinstatement = process_single_rat(rat_id,'Reinstatement',behav_data_folder,neural_data_folder,timestamps_folder,params);

fprintf('Time taken: %.2f minutes',toc/60)
%% Clean up cell-to-index map by excluding neurons that were shitty 

registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';
cell_map = load_cell_map(registration_folder,rat_id);

SA1_badneur = [];
SA2_badneur = 2;
Ext1_badneur = [1,2,3,4,117]; 
Ext2_badneur = [1,2,100,137]; 
Ext3_badneur = [2,3,4,11,125,135]; 
Ext4_badneur = [4,7,104,106,118,131];
Reinstatement_badneur = 107;

[cell_map,~] = calculate_global_excludeIDX(cell_map,{SA1_badneur,SA2_badneur,Ext1_badneur,Ext2_badneur,Ext3_badneur,Ext4_badneur,Reinstatement_badneur});

data_type = 'spikes_denoised';
across_sess_array = make_sessions_array({SA1,SA2,Ext1,Ext2,Ext3,Ext4,Reinstatement},data_type);
ext_day_idx = 3:6;

[population_vectors] = populationVectors_accum(across_sess_array,ext_day_idx,cell_map,'leverOUT',[20 20],'average');

corrMat = across_sess_correlations(across_sess_array,1:7,cell_map,'HLON',[10 1],'average',true,session_names);

corrMat = across_sess_correlations(across_sess_array,1:7,cell_map,'HLON',[10 10],'spike_time',true,session_names);


just_data = population_vectors(3:end,:)';

day_labels = population_vectors(1,:);
trial_labels = population_vectors(2,:);

[~,scores] = pca(just_data);

for i = 3:6
    scatter(scores(day_labels == i & trial_labels == 0,1),scores(day_labels == i & trial_labels == 0,2),'ro');
    hold on; scatter(scores(day_labels == i & trial_labels == 1,1),scores(day_labels == i & trial_labels == 1,2),'bo');
    pause;
end

active_inds = sum(just_data,2) > 3;
active_shit = just_data(active_inds,:);
tf_idf_normed = calcTFIDF(active_shit')';
day_labels_act = day_labels(active_inds);
trial_labels_act = trial_labels(active_inds);

[~,scores] = pca(tf_idf_normed);

for i = 3:6
    scatter(scores(day_labels_act == i & trial_labels_act == 0,1),scores(day_labels_act == i & trial_labels_act == 0,2),'ro');
    hold on; scatter(scores(day_labels_act == i & trial_labels_act == 1,1),scores(day_labels_act == i & trial_labels_act == 1,2),'bo');
    pause;
end


affinity_matrix = 1 - pdist2(tf_idf_normed,tf_idf_normed,'cosine');

mappedX_cos = tsne_p(affinity_matrix,trial_labels_act,2);
plot(mappedX_cos(trial_labels_act==1,1),mappedX_cos(trial_labels_act==1,2),'bo');
hold on;
plot(mappedX_cos(trial_labels_act==0,1),mappedX_cos(trial_labels_act==0,2),'ro');

plot(mappedX_cos(day_labels_act==3,1),mappedX_cos(day_labels_act==3,2),'bo');
hold on; plot(mappedX_cos(day_labels_act==4,1),mappedX_cos(day_labels_act==4,2),'co');
hold on; plot(mappedX_cos(day_labels_act==5,1),mappedX_cos(day_labels_act==5,2),'go');
hold on; plot(mappedX_cos(day_labels_act==6,1),mappedX_cos(day_labels_act==6,2),'ro');


common_neurons_ExtOnly = find(sum(relevant_days(:,2:end) > 0,2) == 4);
common_neurons = find(sum(relevant_days > 0,2) == 5);

common_SA2_Ext1 = find(sum(relevant_days(:,1:2) > 0,2) == 2);

common_neur_SA2 = SA2.spikes_conv(relevant_days(common_neurons,1),:);
correlations_SA2 = corrcoef(common_neur_SA2');

common_neur_SA2 = reshape(common_neur_SA2,length(common_neurons), 799, 100);

common_neur_Ext1 = Ext1.spikes_conv(relevant_days(common_neurons,2),:);
correlations_Ext1 = corrcoef(common_neur_Ext1');

common_neur_Ext1 = reshape(common_neur_Ext1,length(common_neurons), 799, 100);

new_neurons = find(relevant_days(:,1) == 0  & relevant_days(:,2) > 0);

new_neur_Ext1 = Ext1.spikes_conv(relevant_days(new_neurons,2),:);
new_neur_Ext1 = reshape(new_neur_Ext1,length(new_neurons),799,100);
Ext1_avg_new = mean(mean(new_neur_Ext1,3),1);

xt = (1:799) * (1./SA2.Fs);
plot(xt,SA2_avg,'DisplayName','Self-Administration Session 2')
hold on; plot(xt,Ext1_avg,'k','DisplayName','Extinction Session 1 - same neurons')
hold on; plot(xt,Ext1_avg_new,'r','DisplayName','Extinction Session 1 - new neurons')

reshaped_events = permute(reshape(full(Ext1.event_matrix),[],100,length(Ext1.event_names)),[1,3,2]);
failed_trial_idx = sum(squeeze(reshaped_events(:,3,:)),1) == 0;

Ext1_avg_succ = mean(mean(common_neur_Ext1(:,:,~failed_trial_idx),3),1);
Ext1_avg_fail = mean(mean(common_neur_Ext1(:,:,failed_trial_idx),3),1);
Ext1_avg_succ_new = mean(mean(new_neur_Ext1(:,:,~failed_trial_idx),3),1);
Ext1_avg_fail_new = mean(mean(new_neur_Ext1(:,:,failed_trial_idx),3),1);


plot(xt,Ext1_avg_succ,'r','DisplayName','Ext1 successful trials, persistent neurons')
hold on; plot(xt,Ext1_avg_fail,'r','LineWidth',1.5,'DisplayName','Ext1 failed trials, persistent neurons')
hold on; plot(xt,Ext1_avg_succ_new,'b','DisplayName','Ext1 successful trials, new neurons')
hold on; plot(xt,Ext1_avg_fail_new,'b','LineWidth',1.5,'DisplayName','Ext1 failed trials, new neurons')

legend('show')
xlim([1 75])

Ext1_trial_fail = common_neur_Ext1(:,:,failed_trial_idx);
Ext1_trialavg_fail = mean(common_neur_Ext1(:,:,failed_trial_idx),3);

correlated_w_avg = zeros(length(common_neurons),1);

for i = 1:length(common_neurons)
    correlated_w_avg(i) = corr(mean(Ext1_trialavg_fail,1)',Ext1_trialavg_fail(i,:)');
end

[sorted_r,srt] = sort(correlated_w_avg,'descend');
top_regressors = srt(1:30);

Ext1_avg_succ_top = mean(mean(common_neur_Ext1(top_regressors,:,~failed_trial_idx),3),1);
Ext1_avg_fail_top = mean(mean(common_neur_Ext1(top_regressors,:,failed_trial_idx),3),1);

SA2_top30 = common_neur_SA2(top_regressors,:,:);
SA2_avg_top = mean(mean(SA2_top30,3),1);

Ext1_top30 = common_neur_Ext1(top_regressors,:,:);
Ext1_avg_top = mean(mean(Ext1_top30,3),1);

plot(xt,SA2_avg_top,'b','LineWidth',1.5,'DisplayName','Top Extinction Responders on SA2, Average Activity');
hold on; plot(xt,Ext1_avg_top,'r','LineWidth',1.5,'DisplayName','Top Extinction Responders on Ext1, Average Activity');
legend('show')
xlim([1 75])

plot(xt,Ext1_avg_succ_top,'r','DisplayName','Ext1 successful trials, persistent neurons')
hold on; plot(xt,Ext1_avg_fail_top,'r','LineWidth',1.5,'DisplayName','Ext1 failed trials, persistent neurons')
legend('show')
xlim([1 75])

for trial = 1:100
    imagesc(common_neur_Ext1(top_regressors,:,trial));
    if failed_trial_idx(trial) == 1
        type = 'Failed';
    else
        type = 'Pressed';
    end
    title(sprintf('Trial Number %d, Type: %s',trial,type));
    pause; 
end

phase1 = 1:33;
phase2 = 34:66;
phase3 = 67:100;

phase1_activity_failed = common_neur_Ext1(top_regressors,:,failed_trial_idx(phase1));
phase1_activity_succ = common_neur_Ext1(top_regressors,:,~failed_trial_idx(phase1));

phase2_activity_failed = common_neur_Ext1(top_regressors,:,failed_trial_idx(phase2));
phase2_activity_succ = common_neur_Ext1(top_regressors,:,~failed_trial_idx(phase2));

phase3_activity_failed = common_neur_Ext1(top_regressors,:,failed_trial_idx(phase3));
phase3_activity_succ = common_neur_Ext1(top_regressors,:,~failed_trial_idx(phase3));

plot(xt,mean(mean(phase1_activity_failed,3),1),'LineWidth',0.75,'Color',[.6 0 0],'DisplayName','First 33 Failed Trials');
hold on; plot(xt,mean(mean(phase1_activity_succ,3),1),'LineWidth',0.75,'Color',[0 0 .6],'DisplayName','First 33 Successful Trials');

plot(xt,mean(mean(phase2_activity_failed,3),1),'LineWidth',1.25,'Color',[.8 0 0],'DisplayName','Second 33 Failed Trials');
hold on; plot(xt,mean(mean(phase2_activity_succ,3),1),'LineWidth',1.25,'Color',[0 0 .8],'DisplayName','Second 33 Successful Trials');

plot(xt,mean(mean(phase3_activity_failed,3),1),'LineWidth',1.75,'Color',[0.99 0 0],'DisplayName','Third 33 Failed Trials');
hold on; plot(xt,mean(mean(phase3_activity_succ,3),1),'LineWidth',1.75,'Color',[0 0 .99],'DisplayName','Third 33 Successful Trials');

legend('show')
xlim([1 75])

%% do some geometry 

reshaped_events = permute(reshape(full(Ext1.event_matrix),[],100,length(Ext1.event_names)),[1,3,2]);
failed_trial_idx = sum(squeeze(reshaped_events(:,3,:)),1) == 0;

SA2_top_reg = reshape(common_neur_SA2(top_regressors,:,:),length(top_regressors),799*100);
Ext1_top_reg = reshape(common_neur_Ext1(top_regressors,:,:),length(top_regressors),799*100);

correlations_SA2_top_reg = corrcoef(SA2_top_reg');
correlations_Ext1_top_reg = corrcoef(Ext1_top_reg');

scatter(get_lower_tri(correlations_SA2_top_reg),get_lower_tri(correlations_Ext1_top_reg))




relevant_chunk = common_neur_Ext1(top_regressors,250:300,:);
[num_neurons,chunk_size,numTrials] = size(relevant_chunk);

unwrapped = reshape(relevant_chunk,num_neurons,chunk_size*numTrials);

labels = repmat(double(failed_trial_idx)+1,chunk_size,1);
labels = reshape(labels,chunk_size*numTrials,1);

active_idx = sum(unwrapped,1) >= 4;
active_shit = unwrapped(:,active_idx);
labels = labels(active_idx);

[coeff,scores] = pca(active_shit');
scatter(scores(labels==1,1),scores(labels==1,2),15,'bo');
hold on;
scatter(scores(labels==2,1),scores(labels==2,2),15,'ro');


tf_idf_normed = calcTFIDF(active_shit);
tf_idf_normed(isnan(tf_idf_normed)) = 0;
[coeff,scores] = pca(tf_idf_normed');
plot3(scores(labels==1,1),scores(labels==1,2),scores(labels==1,3),'bo');
hold on;
plot3(scores(labels==2,1),scores(labels==2,2),scores(labels==2,3),'ro');


mappedX_euclid = tsne(tf_idf_normed',[],2,10,50);
scatter(mappedX_euclid(labels==1,1),mappedX_euclid(labels==1,2),15,'bo');
hold on;
scatter(mappedX_euclid(labels==2,1),mappedX_euclid(labels==2,2),15,'ro');


affinity_matrix = 1 - pdist2(tf_idf_normed',tf_idf_normed','cosine');
mappedX_cos = tsne_p(affinity_matrix,[],3);
plot3(mappedX_cos(labels==1,1),mappedX_cos(labels==1,2),mappedX_cos(labels==1,3),'bo');
hold on;
plot3(mappedX_cos(labels==2,1),mappedX_cos(labels==2,2),mappedX_cos(labels==2,3),'ro');


%%

reshaped_events = permute(reshape(full(Ext1.event_matrix),[],100,length(Ext1.event_names)),[1,3,2]);
failed_trial_idx = sum(squeeze(reshaped_events(:,3,:)),1) == 0;

relevant_chunk = common_neur_Ext1(top_regressors,260:280,:); % chunk of top regressors centered around the press-aligned bump
[num_neurons,chunk_size,numTrials] = size(relevant_chunk);

unwrapped = reshape(relevant_chunk,num_neurons,chunk_size*numTrials);

labels_which = reshape(repmat([1:100],chunk_size,1),chunk_size*numTrials,1);
labels_type = reshape(repmat(double(failed_trial_idx)+1,chunk_size,1),chunk_size*numTrials,1);

active_idx = sum(unwrapped,1) >= 4;
active_shit = unwrapped(:,active_idx);

labels_which = labels_which(active_idx);
labels_type = labels_type(active_idx);

succ_colors = [zeros(length(find(labels_type==1)),2),linspace(0.6,1,length(find(labels_type==1)))'];
failed_colors = [linspace(0.6,1,length(find(labels_type==2)))',zeros(length(find(labels_type==2)),2)];


time_colorz = winter(length(find(active_idx)));
[coeff,scores] = pca(active_shit');
scatter(scores(:,1),scores(:,2),20,time_colorz);

scatter(scores(labels_type==1,1),scores(labels_type==1,2),25,succ_colors);
hold on;
scatter(scores(labels_type==2,1),scores(labels_type==2,2),25,failed_colors);

tf_idf_normed = calcTFIDF(active_shit);
tf_idf_normed(isnan(tf_idf_normed)) = 0;
[coeff,scores] = pca(tf_idf_normed');
scatter3(scores(labels_type==1,1),scores(labels_type==1,2),scores(labels_type==1,3),20,succ_colors);
hold on;
scatter3(scores(labels_type==2,1),scores(labels_type==2,2),scores(labels_type==2,3),20,failed_colors);

figure;

for i = 1:length(find(active_idx))
    if labels_type(i) == 1
        scatter(scores(i,1),scores(i,2),20,'bo');
    elseif labels_type(i) == 2
        scatter(scores(i,1),scores(i,2),20,'ro');
    end
    xlim([min(scores(:,1)) - 0.2 max(scores(:,1)) + 0.2]);
    ylim([min(scores(:,2)) - 0.2 max(scores(:,2)) + 0.2]);
    title(sprintf('Trial Number %d, trial type: %d',labels_which(i),labels_type(i)));
    pause; hold on;
end

