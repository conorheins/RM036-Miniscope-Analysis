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

%% Pseudo-code workflow:
% For rat in RatID_List:
%     create rat-specific, across-session structure/cell array for that
%     rat's information
%     for session in session_List:
%       1. use 'process_single_rat' to load data from given rat, single
%           session into a SubjectObj class object
%       2. do firing rate average calculations for that session to pull out
%           two trial-averages -- the one for response trials, and the one for
%           no-response trials. Do this PER NEURON, not for all neurons yet
%       3. pull those trial-averages, and put them into the appropriate
%           cell of the trialavg_FRs array
% Now, outside main loop:
% For rat in RatID_List
%     use cell registration map to get indices of persistent neurons, relative to each session
%     for session in session_List:
%       1. create cross-neuron averages to get the average firing rates for
%           the persistent population 
%       2. plot dat shit


behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';
% session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4'};
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};
trial_indices = [10:780];

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

RatIDs = [8];

trialavg_FRs = cell(length(RatIDs),1);

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    temp_rat_data = cell(length(session_names),1);
    
    for jj = 1:length(session_names)
        
        sess_id = session_names{jj};
    
        sess_jj_obj = process_single_rat(rat_id,sess_id,behav_data_folder,neural_data_folder,timestamps_folder,params);
        
        data_type = 'spikes_conv'; % which form of the neural data to use (e.g. 'C','C_decon','spikes_denoised','spikes_conv')
        bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic'
        CI_bounds = [16 84]; % under Gaussianity-assumption, +/- 1 standard deviation of density
        plot_flag = false;
        
        reshaped_events = permute(reshape(full(sess_jj_obj.event_matrix),[],sess_jj_obj.num_trials,length(sess_jj_obj.event_names)),[1,3,2]);
        press_idx = find(cellfun(@(x) ~isempty(x),strfind(sess_jj_obj.event_names,'press')));
        noresponse_trial_idx = sum(squeeze(reshaped_events(:,press_idx,:)),1) == 0;
        
        sess_averages.R_trials_allneur = zeros(size(sess_jj_obj.spikes_denoised,1),length(trial_indices),3);
        sess_averages.R_trials_popavg = zeros(length(trial_indices),3);
        sess_averages.NR_trials_allneur = zeros(size(sess_jj_obj.spikes_denoised,1),length(trial_indices),3);
        sess_averages.NR_trials_popavg = zeros(length(trial_indices),3);
        
        [trial_average,trial_avg_CIs,pop_avg,pop_avg_CIs] = sess_jj_obj.population_activity_wholeTrial(data_type,bound_type,CI_bounds,[],~noresponse_trial_idx,plot_flag);
        sess_averages.R_trials_allneur(:,:,1) = trial_average(:,trial_indices);
        sess_averages.R_trials_allneur(:,:,2:3) = trial_avg_CIs(:,trial_indices,:);
        sess_averages.R_trials_popavg(:,1) = pop_avg(trial_indices)';
        sess_averages.R_trials_popavg(:,2:3) = pop_avg_CIs(:,trial_indices)';
        
        
        if length(find(noresponse_trial_idx)) < 5
            fprintf('Not enough no-response trials to create trial-wise average\n')
            sess_averages.NR_trials_allneur = [];
            sess_averages.NR_trials_popavg = [];
        else
            [trial_average,trial_avg_CIs,pop_avg,pop_avg_CIs] = sess_jj_obj.population_activity_wholeTrial(data_type,bound_type,CI_bounds,[],noresponse_trial_idx,plot_flag);
            sess_averages.NR_trials_allneur(:,:,1) = trial_average(:,trial_indices);
            sess_averages.NR_trials_allneur(:,:,2:3) = trial_avg_CIs(:,trial_indices,:);
            sess_averages.NR_trials_popavg(:,1) = pop_avg(trial_indices)';
            sess_averages.NR_trials_popavg(:,2:3) = pop_avg_CIs(:,trial_indices)';
        end
        
        temp_rat_data{jj} = sess_averages;
        
        Sess_object = sess_jj_obj; clear sess_jj_obj;
        rat_dir = fullfile('rats_processed',sprintf('Rat%d',rat_id));
        write_name = fullfile(rat_dir,sprintf('Rat%d_%s.mat',rat_id,sess_id));
        save(write_name,'Sess_object');
        
    end
    
    trialavg_FRs{ii} = temp_rat_data;
    
end
    

%% plot trial averages -- popoulation wide

trial_indices = [10:780];
Fs = 10.49;

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    rat_data = trialavg_FRs{ii};
    
    visualize_all_sess(rat_data,session_names,trial_indices,Fs)
    pause; close all;
    
end

%%  plot trial averages -- single neurons

trial_indices = [10:780];
Fs = 10.49;
leverOUT_idx = [250:300];
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};
registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';

sess_range = 3:6;

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    cell_map = cell_map(:,sess_range);
   
    all_sessions_badneur = cell(1,length(sess_range));
    
    for jj = 1:length(session_names(sess_range))
        
        sess_id = session_names{sess_range(jj)};
        fprintf('%s\n',sess_id)
        rat_dir = fullfile('rats_processed',sprintf('Rat%d',rat_id));
        fnam = fullfile(rat_dir,sprintf('Rat%d_%s.mat',rat_id,sess_id));
        load(fnam);
        
        all_sessions_badneur{jj} = Sess_object.bad_neurons;
        
    end
    
    [cell_map,~] = calculate_global_excludeIDX(cell_map,all_sessions_badneur);
    common_neurons = find(sum(cell_map > 0,2) == length(sess_range));
    
    rat_data = trialavg_FRs{ii};

    for jj = 1:length(sess_range)
        
        common_neuron_tempidx = cell_map(common_neurons,jj);
        
        if sess_range(jj) <= length(rat_data)
            
            
            R_trials = rat_data{sess_range(jj)}.R_trials_allneur(common_neuron_tempidx,:,:);
            trial_avg_jj = R_trials(:,:,1);
            trial_avg_jj = bsxfun(@rdivide,trial_avg_jj - min(trial_avg_jj,[],2),max(trial_avg_jj,[],2) - min(trial_avg_jj,[],2));
            if jj == 1
                [~,srt] = sort(mean(trial_avg_jj(:,leverOUT_idx),2),'descend');
            else
                srt = srt;
            end
            
            num_neurons = size(trial_avg_jj,1);
            plot_colors = winter(ceil(1.25*num_neurons));
            plot_colors = plot_colors(1:num_neurons,:);
            plot_ndarray(trial_avg_jj(srt,:),1,0,plot_colors);
            title(sprintf('Sorted firing rates for Rat %d, Sess %s',rat_id,session_names{sess_range(jj)}));
            pause; close gcf;
            
        end
    end
    
end

        


%% limit analysis to cells that persist across time

registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';
sess_range = 2:7;

persistent_trialavg_FRs = cell(length(RatIDs),1);

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    cell_map = cell_map(:,sess_range);
   
    all_sessions_badneur = cell(1,length(session_names));
    
    for jj = 1:length(session_names)
        
        sess_id = session_names{jj};
        
        rat_dir = fullfile('rats_processed',sprintf('Rat%d',rat_id));
        fnam = fullfile(rat_dir,sprintf('Rat%d_%s.mat',rat_id,sess_id));
        load(sprintf('Rat%s_%s.mat',num2str(rat_id),session_names{jj}));
        
        all_sessions_badneur{jj} = Sess_object.bad_neurons;
        
    end
    
    [cell_map,~] = calculate_global_excludeIDX(cell_map,all_sessions_badneur);
    common_neurons = find(sum(cell_map > 0,2) == length(sess_range));
    
    num_common = length(common_neurons);
    
    temp_rat_data = cell(length(session_names),1);

    for jj = 1:length(session_names)
        
        sess_id = session_names{jj};
        
        load(sprintf('Rat%s_%s.mat',num2str(rat_id),session_names{jj}));
        
        data_type = 'spikes_conv'; % which form of the neural data to use (e.g. 'C','C_decon','spikes_denoised','spikes_conv')
        bound_type = 'bootstrap'; % options: 'bootstrap','jackknife', or 'analytic'
        CI_bounds = [16 84]; % under Gaussianity-assumption, +/- 1 standard deviation of density
        plot_flag = false;
        
        common_neuron_tempidx = cell_map(common_neurons,jj);

        reshaped_events = permute(reshape(full(Sess_object.event_matrix),[],Sess_object.num_trials,length(Sess_object.event_names)),[1,3,2]);
        press_idx = find(cellfun(@(x) ~isempty(x),strfind(Sess_object.event_names,'press')));
        noresponse_trial_idx = sum(squeeze(reshaped_events(:,press_idx,:)),1) == 0;
        
        sess_averages.R_trials_allneur = zeros(num_common,length(trial_indices),3);
        sess_averages.R_trials_popavg = zeros(length(trial_indices),3);
        sess_averages.NR_trials_allneur = zeros(num_common,length(trial_indices),3);
        sess_averages.NR_trials_popavg = zeros(length(trial_indices),3);
        
        [trial_average,trial_avg_CIs,pop_avg,pop_avg_CIs] = Sess_object.population_activity_wholeTrial(data_type,bound_type,CI_bounds,common_neuron_tempidx,~noresponse_trial_idx,plot_flag);
        sess_averages.R_trials_allneur(:,:,1) = trial_average(:,trial_indices);
        sess_averages.R_trials_allneur(:,:,2:3) = trial_avg_CIs(:,trial_indices,:);
        sess_averages.R_trials_popavg(:,1) = pop_avg(trial_indices)';
        sess_averages.R_trials_popavg(:,2:3) = pop_avg_CIs(:,trial_indices)';
        
        if length(find(noresponse_trial_idx)) < 5
            fprintf('Not enough no-response trials to create trial-wise average\n')
            sess_averages.NR_trials_allneur = [];
            sess_averages.NR_trials_popavg = [];
        else
            [trial_average,trial_avg_CIs,pop_avg,pop_avg_CIs] = Sess_object.population_activity_wholeTrial(data_type,bound_type,CI_bounds,common_neuron_tempidx,noresponse_trial_idx,plot_flag);
            sess_averages.NR_trials_allneur(:,:,1) = trial_average(:,trial_indices);
            sess_averages.NR_trials_allneur(:,:,2:3) = trial_avg_CIs(:,trial_indices,:);
            sess_averages.NR_trials_popavg(:,1) = pop_avg(trial_indices)';
            sess_averages.NR_trials_popavg(:,2:3) = pop_avg_CIs(:,trial_indices)';
        end
        
        temp_rat_data{jj} = sess_averages; clear sess_averages;
        
    end
    
    persistent_trialavg_FRs{ii} = temp_rat_data;
    
end

        
trial_indices = [10:780];
Fs = 10.49;

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    rat_data = persistent_trialavg_FRs{ii};
    
    visualize_all_sess(rat_data,session_names,trial_indices,Fs)
    pause; close all;
    
end
    
    
    
    
   
    
    




    

