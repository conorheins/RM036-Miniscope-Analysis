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

%%
rat_id = 23;

registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';

cell_map = load_cell_map(registration_folder,rat_id);

behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};
sess_2align = [1:7];

neuron_ids = find(sum(cell_map(:,sess_2align) > 0,2) == length(sess_2align));

numPCs = 10;

session_names_2align = session_names(sess_2align);

trialavg_FR = zeros(length(session_names_2align),799);
CI_FR = zeros(length(session_names_2align),2,799);

trialavg_PC = zeros(length(session_names_2align),799,numPCs);
CI_PC = zeros(length(session_names_2align),2,799,numPCs);

event_matrices = cell(1,length(session_names_2align));

%% actually run the damn thang

for sess = 1:length(session_names_2align)
    
    session_name = session_names_2align{sess};

    currRat = SubjectObj;

    currRat.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);
    
    %% remove bad trials, if any
    
%     if strcmp(session_name,'Ext3')
%         currRat.remove_trials(40); % for instance, Rat21, Ext3 had a big artifact at Trial 40
%     end
    
    % narrow down to only those neurons that persist across the
    % specified sessions
    currRat.filter_neurons('keep',cell_map(neuron_ids,sess));
    
    %% Process fluorescence traces and estimate spiking (deconvolution, spike-thresholding, etc.)
    
    % deconvolve data
    
    currRat.deconvTrials()
    
    % denoise spikes estimated via deconvolution with a local/global
    % thresholding process
    
    dFF_thr = 3.5;
    look_back_f = 30;
    look_ahead_f = 30;
    currRat.denoiseSpikes(currRat.C,currRat.spikes,dFF_thr,look_back_f,look_ahead_f)
    
    % convolve spikes
    
    spike_data = currRat.spikes_denoised;
    transpose_flag = 0;
    kernel = normpdf(-3:3,0,1.5);
    currRat.convSpikes(spike_data,transpose_flag,kernel)
    
    bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic' 
    CI_bounds = [16 84]; %under Gaussianity-assumption, +/- 1 SD

    [~,trialavg_FR(sess,:),CI_FR(sess,:,:)] = currRat.population_activity_wholeTrial('spikes_conv',bound_type,CI_bounds,[],0);
    
    [trialavg_PC(sess,:,:), CI_PC(sess,:,:,:)] = currRat.PCA_trialplot('spikes_conv',numPCs,'bootstrap',CI_bounds,0);
    
    event_matrices{sess} = currRat.event_matrix;
    
end


save('fullPCA_allsessions_persist.mat','trialavg_PC','CI_PC','event_matrices','session_names_2align')
save('avgFiringRate_allsessions_persist.mat','trialavg_FR','CI_FR','event_matrices','session_names_2align')