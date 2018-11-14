%% Object-Oriented Script Version of RM036 Analysis

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

%% Process fluorescence traces and estimate spiking (smoothing, deconvolution, spike-thresholding, etc.)

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
currRat.denoiseSpikes(currRat.C,threshold_timeseries(currRat.spikes,'mean',2.5),dFF_thr,look_back_f,look_ahead_f,counts_flag)

% visualize results of spike-denoising

currRat.overlay_spikes_calc(currRat.C,currRat.spikes_denoised)

% convolve spikes

spike_data = currRat.spikes_denoised;
transpose_flag = 0; 
kernel = normpdf(-3:3,0,1.5);
currRat.convSpikes(spike_data,transpose_flag,kernel)


%% align data to behaviorally-relevant timestamps

before_align = 100;
after_align = 200;
currRat.lock2events('spikes_conv',before_align,after_align);

% test for selectivity of neurons to specific trial- or behavioral-events

num_shuffles = 1000;
alpha = 0.05;
surround_time = [100 100]; % time around which to average data for pre- and post-event distributions
currRat.compute_selectivity('permutation',surround_time,num_shuffles,alpha)

currRat.compute_fidelity();

%% Displays

if strcmp(currRat.meta_data.data_type_align,'spikes') || strcmp(currRat.meta_data.data_type_align,'spikes_denoised')
    for i = 1:length(currRat.event_names)
        currRat.displayRasters(currRat.event_names{i},[],1:100);
    end
else
    for i = 1:length(currRat.event_names)
        currRat.display_multiunit(currRat.event_names{i},surround_time,[],[],0)
    end
end

% some saving in case you want to show trial_average firing rates of significant neurons 
% for i = 1:length(currRat.event_names)
%     
%     event_idx = cellfun(@(x) ~isempty(x), strfind(currRat.event_names,currRat.event_names{i}));
%     temp = event_locked{event_idx};
%     
%     dat2display = temp(currRat.stats_results.Modulation_direction(:,event_idx) == 1, :, :);
%     dat2display = mean(dat2display,3);
%     
%     plot(dat2display');
%     title(['Average firing rate aligned to event: ',currRat.event_names{i}])
%     hold on; plot( [currRat.meta_data.event_tmsp,currRat.meta_data.event_tmsp],[min(dat2display(:)),max(dat2display(:))],'r--','LineWidth',1.5);
%     axis tight;
%     saveas(gcf,[currRat.event_names{i},'_lockedFR.png']);
%     hold off;
%     pause; close gcf
%     
% end

%% show different statistical summaries with different combinations of 

currRat.stats_summary()

figname = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/stats_summary_spikes_conv.png';
saveas(gcf,figname);


%% plot signal correlations

sig_corr = currRat.compute_signal_correlation();

min_corr = -1;
max_corr = 1;

figure;
caxis([min_corr max_corr])
for cond = 1:length(sig_corr)
    subplot(3,3,cond)
    imagesc(sig_corr{cond})
    title(currRat.event_names{cond})
end

figname = ['/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/sig_correlations_',currRat.meta_data.data_type_align,'.png'];
saveas(gcf,figname);

%% plot noise correlations

min_corr = -1;
max_corr = 1;

figure;
caxis([min_corr max_corr])
noise_corr = currRat.compute_noise_correlations();
imagesc(noise_corr);
title('Noise correlations')

figname = ['/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/noise_correlations_',currRat.meta_data.data_type_align,'.png'];
saveas(gcf,figname);

%% plot response fidelity distributions

currRat.fidelity_distribution()

figname = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/response_fidelity_spikes_conv.png';
saveas(gcf,figname);

%%  Test different windows over which to compute neuronal selectivity
alpha = 0.05;
num_shuffles = 1000;
time_vector = 5:5:100;
timescale_results = zeros(length(time_vector),length(currRat.event_names),3);

for ii = 1:length(time_vector)
    currRat.compute_selectivity('permutation',[time_vector(ii) time_vector(ii)],num_shuffles,alpha);
    timescale_results(ii,:,1) = sum(currRat.stats_results.SigMatrix,1);
    timescale_results(ii,:,2) = sum(currRat.stats_results.Modulation_direction == 1,1);
    timescale_results(ii,:,3) = sum(currRat.stats_results.Modulation_direction == -1,1);
end

% plot the results

figure;
for cond = 1:size(timescale_results,2)
    plot(time_vector,squeeze(timescale_results(:,cond,1)),'DisplayName',currRat.event_names{cond})
    hold on;
end
title('Effect of increasing window size on effect detection')
xlabel('Time window size')
ylabel('Number of significantly-modulated units')
legend('show')

figname = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/window_length_significance.png';
saveas(gcf,figname);

%% event locked population activity to certain epochs

event_names = currRat.event_names;
event_tmsp = currRat.meta_data.event_tmsp;
for i = 1:length(event_names)
    event_of_interest = event_names{i};
    event_idx = find(cellfun(@(x) ~isempty(x), strfind(currRat.event_names,event_of_interest)));
    traces = currRat.event_locked{event_idx};
    time_vect = 
    plot(mean(sum(traces,1),3)); axis tight;
    ylimz = ylim;
    hold on; h2 = plot([event_tmsp, event_tmsp],[ylimz(1),ylimz(2)],'r--','LineWidth',1.5);
    legend(h2,sprintf('Time of event: %s',event_of_interest))
    title(sprintf('Population activity locked to %s',event_of_interest))
    pause; hold off;
end

%% sequence embedding

events2examine = {'HLON','cueON','pellets','bbk','leverOUT'};

labels = [];
mean_activity_patterns = [];
sequence_patterns = [];

for event_i = 1:length(events2examine)
    
    if strcmp(events2examine{event_i},'leverOUT')
        temp = currRat.pattern_vectors(events2examine{event_i},'spikes_conv',[100 0]);
        mean_activity_patterns = [mean_activity_patterns; temp];
        
        temp = currRat.sequence_vectors(events2examine{event_i},'spikes_conv',[100 0]);
        sequence_patterns = [sequence_patterns; temp];
    else
        
        temp = currRat.pattern_vectors(events2examine{event_i},'spikes_conv',[0 100]);
        mean_activity_patterns = [mean_activity_patterns; temp];
        
        temp = currRat.sequence_vectors(events2examine{event_i},'spikes_conv',[0 100]);
        sequence_patterns = [sequence_patterns; temp];
    end
    
    labels = [labels; (event_i) * ones(size(temp,1),1)];
    clear temp;
    
end

[eigenvects_mean,scores_mean] = pca(mean_activity_patterns);
[eigenvects_seq,scores_seq] = pca(sequence_patterns);

reduced_space_mean = scores_mean(:,1:5);
reduced_space_seq = scores_seq(:,1:5);

unique_labels = unique(labels);

TSS_mean = sum ( sum((reduced_space_mean - mean(reduced_space_mean,1)).^2,2), 1);
TSS_seq = sum ( sum((reduced_space_seq - mean(reduced_space_seq,1)).^2,2), 1);

WCSS_mean = zeros(1,length(unique_labels));
WCSS_seq = zeros(1,length(unique_labels));

for lab = 1:length(unique_labels)
    
    temp = reduced_space_mean(labels==unique_labels(lab),:);
    WCSS_mean(lab) = sum( sum( (temp - mean(temp,1)).^2,2), 1);
    
    temp = reduced_space_seq(labels==unique_labels(lab),:);
    WCSS_seq(lab) = sum( sum( (temp - mean(temp,1)).^2,2), 1);
    
end

BSS_mean = TSS_mean - sum(WCSS_mean);
BSS_seq = TSS_seq - sum(WCSS_seq);


%% Results displaying why sequence embedding is better

figure(1);
bar([BSS_mean/TSS_mean;BSS_seq/TSS_seq],'k');
ylabel('Cluster separation (Between Stimulus SSE)')
legendz = {'Firing rate vectors', 'Sequence vectors'};
set(gca,'xticklabel',legendz)
title('Stimulus separability for different embeddings')

figure(2)

individual_WCSS_mean = WCSS_mean./TSS_mean;
individual_WCSS_seq = WCSS_seq./TSS_seq;

bar([individual_WCSS_mean;individual_WCSS_seq]','grouped') 
ylabel('Within stimulus SSE (normalized by total SSE)')
legend('Firing rate vectors','Sequence vectors')
set(gca,'xticklabel',events2examine)

%% population average activity with aligned behavioral events

CI_bounds = [2.5 97.5]; %under Gaussianity-assumption, +/- 1 SD
bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
[trial_average,pop_average,CIs] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,[]);

%% plot average activity over different phases of the session

first_third = [1:33];
CI_bounds = [16 84]; %under Gaussianity-assumption, +/- 1 SD
bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
[trial_average1,pop_average1,CIs1] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,first_third);

second_third = [34:66];
CI_bounds = [16 84]; %under Gaussianity-assumption, +/- 1 SD
bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
[trial_average2,pop_average2,CIs2] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,second_third);

third_third = [67:100];
CI_bounds = [16 84]; %under Gaussianity-assumption, +/- 1 SD
bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
[trial_average3,pop_average3,CIs3] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,third_third);


%% compute fidelity distributions with larger bins (20 trials rather than 10 trials)

currRat.lock2events('spikes_conv',100,200)
currRat.compute_selectivity('permutation',[100 200],1000,0.01);
currRat.stats_summary()
currRat.compute_fidelity()
[~] = currRat.fidelity_distribution([],0.2);

figname = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/response_fidelity_spikes_conv_20TRIALBINS.png';
saveas(gcf,figname);

% now do cumulative distributions with original settings (all event types,
% but 5 trial bins
histogram_data = currRat.fidelity_distribution();
currRat.fidelity_dist_cdf(histogram_data)

%% threshold signal correlation matrices and then make a new color-coded one with all behavioral events 

sig_corr = currRat.compute_signal_correlation([],100);
combined_matrix = zeros([size(sig_corr{1}),3]);
colors = parula(length(currRat.event_names));

for event_idx = 1:length(currRat.event_names)
    
    temp = get_lower_tri(sig_corr{event_idx});
    bins = -1:0.01:1;
    cd = histc(temp,bins);
    cd = cumsum(cd/sum(cd));
    
    temp = sig_corr{event_idx} > bins(find(cd>0.99,1)); % threshold with 99% percentile of correlation values
    temp(1:(size(temp,1)+1):end) = 0; % set diagonal to 0
    
    colorized = zeros(size(combined_matrix));
    
    [rows,cols] = find(temp);
    
    for id = 1:length(rows)
        colorized(rows(id),cols(id),:) = colors(event_idx,:);
    end
    
    combined_matrix = combined_matrix + colorized;

end

%% overlap signal correlations for subset of events

events2check = {'leverIN','press','pellets'};
event_idx = zeros(1,3);
for ii = 1:length(events2check)
    event_idx(ii) = find(cellfun(@(x) ~isempty(x), strfind(currRat.event_names,events2check{ii})));
end

combined_matrix = zeros([size(sig_corr{1}),3]);
% colors = parula(length(events2check));
colors = [ 1 0 0; 0 1 0; 0 0 1];

for ii = 1:length(events2check)
    
    temp = get_lower_tri(sig_corr{event_idx(ii)});
    bins = -1:0.01:1;
    cd = histc(temp,bins);
    cd = cumsum(cd/sum(cd));
    
    temp = sig_corr{event_idx(ii)} > bins(find(cd>0.99,1)); % threshold with 99% percentile of correlation values
    temp(1:(size(temp,1)+1):end) = 0; % set diagonal to 0
    
    colorized = zeros(size(combined_matrix));
    
    [rows,cols] = find(temp);
    
    for id = 1:length(rows)
        colorized(rows(id),cols(id),:) = colors(ii,:);
    end
    
    combined_matrix = combined_matrix + colorized;

end 

%% PCA on full-timeseries

[trial_average, CIs] = currRat.PCA_trialplot('spikes_conv',1,'bootstrap',[16 84]);
