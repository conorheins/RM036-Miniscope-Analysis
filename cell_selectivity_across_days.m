%% paths
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

%%

rat_id = 21;

registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';

cell_map = load_cell_map(registration_folder,rat_id);

behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};
sess_2align = [1:2];

neuron_ids = find(sum(cell_map(:,sess_2align) > 0,2) == length(sess_2align));

session_names_2align = session_names(sess_2align);

session_objects = cell(1,length(session_names_2align));

%% actually run the damn thang

for sess = 1:length(session_names_2align)
    
    session_name = session_names_2align{sess};

    currRat = SubjectObj;

    currRat.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);
    
    %% remove bad trials, if any
    
    if strcmp(session_name,'Ext3')
        currRat.remove_trials(40); % for instance, Rat21, Ext3 had a big artifact at Trial 40
    end
    
    % narrow down to only those neurons that persist across the
    % specified sessions
    currRat.filter_neurons('keep',cell_map(neuron_ids,sess_2align(sess)));
    
    %% Process fluorescence traces and estimate spiking (deconvolution, spike-thresholding, etc.)
    
    % deconvolve data
    
    currRat.deconvTrials()
    
    % denoise spikes estimated via deconvolution with a local/global
    % thresholding process
    
    dFF_thr = 6;
    look_back_f = 10;
    look_ahead_f = 50;
    currRat.denoiseSpikes(currRat.C,threshold_timeseries(currRat.spikes,'median',2.5),dFF_thr,look_back_f,look_ahead_f)
    
    % convolve spikes
    
    spike_data = currRat.spikes_denoised;
    transpose_flag = 0;
    kernel = normpdf(-3:3,0,1.5);
    currRat.convSpikes(spike_data,transpose_flag,kernel)
    
    % align to events for selectivity testing
    
    before_align = 100;
    after_align = 200;
    currRat.lock2events('spikes_conv',before_align,after_align);
    
    % test for selectivity of neurons to specific trial- or behavioral-events
    
    num_shuffles = 1000;
    alpha = 0.05;
    surround_time = [50 50]; % time around which to average data for pre- and post-event distributions
    currRat.compute_selectivity('permutation',surround_time,num_shuffles,alpha)
    
    currRat.compute_fidelity();

    session_objects{sess} = currRat;
    
end

%% compare single-neuron properties across days

[rate_array] = compare_rates_acrossdays(session_objects,session_names_2align,'spikes_denoised',true);

events_to_test = {'leverOUT','cueOFF','HLON','pellets'};
[distance_matrix,distance_matrix_shuffled] = compare_sig_units(session_objects,events_to_test);

SA1_threshold = stats_results{1}.pVals(stats_results{1}.SigMatrix == 1 & stats_results{2}.SigMatrix == 1);
SA2_threshold = stats_results{2}.pVals(stats_results{1}.SigMatrix == 1 & stats_results{2}.SigMatrix == 1);

for event_i = 1:9
    
%     pVals_SA1 = stats_results{1}.pVals(:,event_i);
%     pVals_SA2 = stats_results{2}.pVals(:,event_i);

    fidelity_SA1 = stats_results{1}.fidelity_matrix(:,event_i);
    fidelity_SA2 = stats_results{2}.fidelity_matrix(:,event_i);
    
    sig_SA1 = stats_results{1}.SigMatrix(:,event_i);
    sig_SA2 = stats_results{2}.SigMatrix(:,event_i);
    
    filter_vect = sig_SA1 == 1 & sig_SA2 == 1;
%     scatter(pVals_SA1(filter_vect),pVals_SA2(filter_vect));
    scatter(fidelity_SA1(filter_vect),fidelity_SA2(filter_vect));
    
    pause; 
end
    
    
