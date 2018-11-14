%% Run series of PCAs with bootstrapped confidence intervals on all of Rat21's sessions. Show trial averaged PCA


%% paths
clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);

%% set up parameters and extract principal components, session by session

rat_id = 21;
behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

numPCs = 10;

trial_indices = [10:780];

trialavg_FR = zeros(length(session_names),799);
CI_FR = zeros(length(session_names),2,799);

trialavg_PC = zeros(length(session_names),length(trial_indices),numPCs);
CI_PC = zeros(length(session_names),2,length(trial_indices),numPCs);

event_matrices = cell(1,length(session_names));

for sess = 1:length(session_names)
    
    session_name = session_names{sess};

    currRat = SubjectObj;

    currRat.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);
    
    %% remove bad trials, if any
    
    if strcmp(session_name,'Ext3')
        currRat.remove_trials(40); % for instance, Rat21, Ext3 had a big artifact at Trial 40
    end
    
    %% Process fluorescence traces and estimate spiking (deconvolution, spike-thresholding, etc.)
    
    % deconvolve data
    
    currRat.deconvTrials()
    
    % denoise spikes estimated via deconvolution with a local/global
    % thresholding process
    
    dFF_thr = 6;
    look_back_f = 10;
    look_ahead_f = 50;
    counts_flag = true;
    currRat.denoiseSpikes(currRat.C,threshold_timeseries(currRat.spikes,'mean',1.5),dFF_thr,look_back_f,look_ahead_f,counts_flag)

    % convolve spikes
    
    spike_data = sqrt(currRat.spikes_denoised);
    transpose_flag = 0;
    kernel = normpdf(-3:3,0,1.5);
    currRat.convSpikes(spike_data,transpose_flag,kernel)
    
    bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
    CI_bounds = [16 84]; %under Gaussianity-assumption, +/- 1 SD
    plot_flag = false;

    [~,trialavg_FR(sess,:),CI_FR(sess,:,:)] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,[],0);
    
    [trialavg_PC(sess,:,:), CI_PC(sess,:,:,:)] = currRat.PCA_trialplot('spikes_conv',numPCs,trial_indices,'bootstrap',CI_bounds,plot_flag);
    
    event_matrices{sess} = currRat.event_matrix;
    
end

%% save results

save('fullPCA_allsessions.mat','trialavg_PC','CI_PC','event_matrices','session_names')
save('avgFiringRate_allsessions.mat','trialavg_FR','CI_FR','event_matrices','session_names')