%% Script 12.10.2018

% Goal: group average activity into a few behavioral bins and look at the
% average activity in those bins across time, for different neurons
% 'shared' vs. 'unshared'. Also for response/non-response trials

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

registration_folder = fullfile(base_directory,'CellReg-1.3.4/data');


%% parameters for plotting/rats to use/etc.
RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

data_type = 'spikes_conv'; 

SA_sessions = {'SA1','SA2'};

events_of_interest = {'HLON','leverOUT','leverIN','pellets'}; 

temporal_bins = {[-5,5],[-5,5],[-5,5],[-5,5]};
Fs = 10.49; % sampling rate, in frames/second

ExtPhase1 = {'Ext1','Ext2'};
ExtPhase2 = {'Ext3','Ext4'};

%%
for ii = 1:length(RatIDs)

    %% set up params
    rat_id = RatIDs(ii);
    rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    exist_sessions = {};
    all_data_files = dir([rat_folder,filesep,'*.mat']);
    all_data_files = {all_data_files(:).name};
    for dat_i = 1:length(all_data_files)
        num_remove = length(sprintf('Rat%d_',rat_id))+1;
        exist_sessions{dat_i} = all_data_files{dat_i}(num_remove:end-4);
    end
    
    sess_ids = [];
    for kk = 1:length(session_names)
        if any(strcmp(exist_sessions,session_names{kk}))
            sess_ids = [sess_ids,kk];
        end
    end
    
    % find neurons that show up in the intersection of SA1/SA2
    
    sa_idx = [];
    for sa_i = 1:length(SA_sessions)
        sa_idx = [sa_idx,find(strcmp(SA_sessions{sa_i},session_names(sess_ids)))];
    end
    
    sa_neurons = find(sum(cell_map(:,sa_idx)>0,2) == length(sa_idx));
    
    accum_data_SA = cell(length(events_of_interest),1);
    
    %% accumulate data from self-administration sessions
    
    for kk = 1:length(sa_idx)
        
        sess_name = session_names{sess_ids(sa_idx(kk))};
        
        sess_specific_idx = cell_map(sa_neurons,sess_ids(sa_idx(kk)));
        
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        reshaped_events = reshape(full(Sess_object.event_matrix),799,Sess_object.num_trials,size(Sess_object.event_matrix,2));
        reshaped_neur = reshape(Sess_object.(data_type)(sess_specific_idx,1:(799*Sess_object.num_trials)),[length(sess_specific_idx), 799,Sess_object.num_trials]);

        for event_i = 1:length(events_of_interest)
            
            [event_tmsp,trial_id] = find(squeeze(reshaped_events(:,:,strcmp(events_of_interest{event_i},Sess_object.event_names))));
            [trial_id,temp] = unique(trial_id,'first');
            event_tmsp = event_tmsp(temp);
            event_i_data = zeros(length(sess_specific_idx),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(event_tmsp));
            
            for jj = 1:length(event_tmsp)
                temp = event_tmsp(jj);
                win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                event_i_data(:,:,jj) = reshaped_neur(:,win_edges(1):win_edges(2),trial_id(jj));
            end
            
            accum_data_SA{event_i} = cat(3,accum_data_SA{event_i},event_i_data);
            
        end
    end
    
    %% Ext Phase 1 (Ext1 - Ext2)    
    
    % find neurons that show up in the intersection of Ext1/Ext2
    % (so-called 'ExtPhase1')
    
    extph1_idx = [];
    for ext_i = 1:length(ExtPhase1)
        extph1_idx = [extph1_idx,find(strcmp(ExtPhase1{ext_i},session_names(sess_ids)))];
    end
    
    extph1_neurons = find(sum(cell_map(:,extph1_idx)>0,2) == length(extph1_idx));
    
    % now find sa_neurons that are also in extph1 neurons
    
    shared_extph1 = extph1_neurons(ismember(extph1_neurons,sa_neurons));
    unshared_extph1 = extph1_neurons(~ismember(extph1_neurons,sa_neurons));
    
    accum_data_ExtPh1_s = cell(length(events_of_interest),1);
    accum_data_ExtPh1_us = cell(length(events_of_interest),1);
    
    %% loop through both sessions of ExtPhase1, accumulating adjacent sessions into a single session
    for kk = 1:length(extph1_idx)
        
        sess_name = session_names{sess_ids(extph1_idx(kk))};
        
        sess_specific_idx_s = cell_map(shared_extph1,sess_ids(extph1_idx(kk)));
        sess_specific_idx_us = cell_map(unshared_extph1,sess_ids(extph1_idx(kk)));
        
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        reshaped_events = reshape(full(Sess_object.event_matrix),799,Sess_object.num_trials,size(Sess_object.event_matrix,2));
        
        R_trial_bool = logical(sum(squeeze(reshaped_events(:,:,strcmp('cueON',Sess_object.event_names))),1));
        NR_trial_bool = ~R_trial_bool;
        
        reshaped_neur_s = reshape(Sess_object.(data_type)(sess_specific_idx_s,1:(799*Sess_object.num_trials)),[length(sess_specific_idx_s), 799,Sess_object.num_trials]);
        reshaped_neur_us = reshape(Sess_object.(data_type)(sess_specific_idx_us,1:(799*Sess_object.num_trials)),[length(sess_specific_idx_us), 799,Sess_object.num_trials]);

        for event_i = 1:length(events_of_interest)
            
            [event_tmsp,trial_id] = find(squeeze(reshaped_events(:,:,strcmp(events_of_interest{event_i},Sess_object.event_names))));
            
            if ~isempty(event_tmsp)
                [trial_id,temp] = unique(trial_id,'first');
                event_tmsp = event_tmsp(temp);
                
                R_trials_tmsp = event_tmsp(R_trial_bool);
                R_trials_id = trial_id(R_trial_bool);
                
                NR_trials_tmsp = event_tmsp(NR_trial_bool);
                NR_trials_id = trial_id(NR_trial_bool);
                
                event_i_data_Rs = zeros(length(sess_specific_idx_s),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(R_trials_id));
                event_i_data_NRs = zeros(length(sess_specific_idx_s),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(NR_trials_id));
                
                event_i_data_Rus = zeros(length(sess_specific_idx_us),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(R_trials_id));
                event_i_data_NRus = zeros(length(sess_specific_idx_us),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(NR_trials_id));
                
                for jj = 1:length(R_trials_id)
                    temp = R_trials_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data_Rs(:,:,jj) = reshaped_neur_s(:,win_edges(1):win_edges(2),R_trials_id(jj));
                    event_i_data_Rus(:,:,jj) = reshaped_neur_us(:,win_edges(1):win_edges(2),R_trials_id(jj));
                end
                
                for jj = 1:length(NR_trials_id)
                    temp = NR_trials_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data_NRs(:,:,jj) = reshaped_neur_s(:,win_edges(1):win_edges(2),NR_trials_id(jj));
                    event_i_data_NRus(:,:,jj) = reshaped_neur_us(:,win_edges(1):win_edges(2),NR_trials_id(jj));
                end
                
                accum_data_ExtPh1_s{event_i} = cell(1,2);
                accum_data_ExtPh1_s{event_i}{1} = cat(3,accum_data_ExtPh1_s{event_i}{1},event_i_data_Rs);
                accum_data_ExtPh1_s{event_i}{2} = cat(3,accum_data_ExtPh1_s{event_i}{2},event_i_data_NRs);
                
                accum_data_ExtPh1_us{event_i} = cell(1,2);
                accum_data_ExtPh1_us{event_i}{1} = cat(3,accum_data_ExtPh1_us{event_i}{1},event_i_data_Rus);
                accum_data_ExtPh1_us{event_i}{2} = cat(3,accum_data_ExtPh1_us{event_i}{2},event_i_data_NRus);
            end
                
        end
    end
    
    %% Ext Phase 2 (Ext3 - Ext4)    
    
    % find neurons that show up in the intersection of Ext3/Ext4
    % (so-called 'ExtPhase2')
    
    extph2_idx = [];
    for ext_i = 1:length(ExtPhase2)
        extph2_idx = [extph2_idx,find(strcmp(ExtPhase2{ext_i},session_names(sess_ids)))];
    end
    
    extph2_neurons = find(sum(cell_map(:,extph2_idx)>0,2) == length(extph2_idx));
    
    % now find sa_neurons that are also in extph1 neurons
    
    shared_extph2 = extph2_neurons(ismember(extph2_neurons,sa_neurons));
    unshared_extph2 = extph2_neurons(~ismember(extph2_neurons,sa_neurons));
    
    accum_data_ExtPh2_s = cell(length(events_of_interest),1);
    accum_data_ExtPh2_us = cell(length(events_of_interest),1);
    
    %% loop through both sessions of ExtPhase2, accumulating adjacent sessions into a single session
   
    for kk = 1:length(extph2_idx)
        
        sess_name = session_names{sess_ids(extph2_idx(kk))};
        
        sess_specific_idx_s = cell_map(shared_extph2,sess_ids(extph2_idx(kk)));
        sess_specific_idx_us = cell_map(unshared_extph2,sess_ids(extph2_idx(kk)));
        
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        reshaped_events = reshape(full(Sess_object.event_matrix),799,Sess_object.num_trials,size(Sess_object.event_matrix,2));
        
        R_trial_bool = logical(sum(squeeze(reshaped_events(:,:,strcmp('cueON',Sess_object.event_names))),1));
        NR_trial_bool = ~R_trial_bool;
        
        reshaped_neur_s = reshape(Sess_object.(data_type)(sess_specific_idx_s,1:(799*Sess_object.num_trials)),[length(sess_specific_idx_s), 799,Sess_object.num_trials]);
        reshaped_neur_us = reshape(Sess_object.(data_type)(sess_specific_idx_us,1:(799*Sess_object.num_trials)),[length(sess_specific_idx_us), 799,Sess_object.num_trials]);

        for event_i = 1:length(events_of_interest)
            
            [event_tmsp,trial_id] = find(squeeze(reshaped_events(:,:,strcmp(events_of_interest{event_i},Sess_object.event_names))));
            
            if ~isempty(event_tmsp)
                [trial_id,temp] = unique(trial_id,'first');
                event_tmsp = event_tmsp(temp);
                
                R_trials_tmsp = event_tmsp(R_trial_bool);
                R_trials_id = trial_id(R_trial_bool);
                
                NR_trials_tmsp = event_tmsp(NR_trial_bool);
                NR_trials_id = trial_id(NR_trial_bool);
                
                event_i_data_Rs = zeros(length(sess_specific_idx_s),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(R_trials_id));
                event_i_data_NRs = zeros(length(sess_specific_idx_s),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(NR_trials_id));
                
                event_i_data_Rus = zeros(length(sess_specific_idx_us),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(R_trials_id));
                event_i_data_NRus = zeros(length(sess_specific_idx_us),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(NR_trials_id));
                
                for jj = 1:length(R_trials_id)
                    temp = R_trials_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data_Rs(:,:,jj) = reshaped_neur_s(:,win_edges(1):win_edges(2),R_trials_id(jj));
                    event_i_data_Rus(:,:,jj) = reshaped_neur_us(:,win_edges(1):win_edges(2),R_trials_id(jj));
                end
                
                for jj = 1:length(NR_trials_id)
                    temp = NR_trials_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data_NRs(:,:,jj) = reshaped_neur_s(:,win_edges(1):win_edges(2),NR_trials_id(jj));
                    event_i_data_NRus(:,:,jj) = reshaped_neur_us(:,win_edges(1):win_edges(2),NR_trials_id(jj));
                end
                
                accum_data_ExtPh2_s{event_i} = cell(1,2);
                accum_data_ExtPh2_s{event_i}{1} = cat(3,accum_data_ExtPh2_s{event_i}{1},event_i_data_Rs);
                accum_data_ExtPh2_s{event_i}{2} = cat(3,accum_data_ExtPh2_s{event_i}{2},event_i_data_NRs);
                
                accum_data_ExtPh2_us{event_i} = cell(1,2);
                accum_data_ExtPh2_us{event_i}{1} = cat(3,accum_data_ExtPh2_us{event_i}{1},event_i_data_Rus);
                accum_data_ExtPh2_us{event_i}{2} = cat(3,accum_data_ExtPh2_us{event_i}{2},event_i_data_NRus);
            end
                
        end
    end
    
  
    %% now plot all that shit
    
    % N.B. use this line to get boostrapped intervals (equivalent to density of +/- 1 SEM):
    % errors_SA = bootci(1000,{@mean,mean(temp,3)},'type','per','alpha',0.32);

    
    for event_i = 1:length(events_of_interest)
        
        t_axis = -5:1/Fs:5;

        
        %% RESPONSE TRIALS
        
        % start with self-admin trials
        
        temp = accum_data_SA{event_i};
        means_SA = mean(mean(temp,3),1);
        errors_SA = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        % now move onto Extinction Phase 1
        
        temp = accum_data_ExtPh1_s{event_i}{1};
        means_ExtP1_Rs = mean(mean(temp,3),1);
        errors_ExtP1_Rs = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        temp = accum_data_ExtPh1_us{event_i}{1};
        means_ExtP1_Rus = mean(mean(temp,3),1);
        errors_ExtP1_Rus = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        % now move onto Extinction Phase 2
        
        temp = accum_data_ExtPh2_s{event_i}{1};
        means_ExtP2_Rs = mean(mean(temp,3),1);
        errors_ExtP2_Rs = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        temp = accum_data_ExtPh2_us{event_i}{1};
        means_ExtP2_Rus = mean(mean(temp,3),1);
        errors_ExtP2_Rus = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        
        figure(1);
        
        subplot(311)
        mseb(t_axis,means_SA,errors_SA)
        subplot(312)
        mseb(t_axis,[means_ExtP1_Rs;means_ExtP1_Rus],[errors_ExtP1_Rs;errors_ExtP1_Rus])
        subplot(313)
        mseb(t_axis,[means_ExtP2_Rs;means_ExtP2_Rus],[errors_ExtP2_Rs;errors_ExtP2_Rus])
        
        
        %% NON RESPONSE TRIALS
                
        % start with self-admin trials
        
%         temp = accum_data_SA{event_i};
%         means_SA = mean(mean(temp,3),1);
%         errors_SA = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        % now move onto Extinction Phase 1
        
        temp = accum_data_ExtPh1_s{event_i}{2};
        means_ExtP1_NRs = mean(mean(temp,3),1);
        errors_ExtP1_NRs = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        temp = accum_data_ExtPh1_us{event_i}{2};
        means_ExtP1_NRus = mean(mean(temp,3),1);
        errors_ExtP1_NRus = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        % now move onto Extinction Phase 2
        
        temp = accum_data_ExtPh2_s{event_i}{2};
        means_ExtP2_NRs = mean(mean(temp,3),1);
        errors_ExtP2_NRs = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        temp = accum_data_ExtPh2_us{event_i}{2};
        means_ExtP2_NRus = mean(mean(temp,3),1);
        errors_ExtP2_NRus = std(mean(temp,3),0,1)./sqrt(size(temp,1));
        
        figure(2);
        subplot(311)
        subplot(312)
        mseb(t_axis,[means_ExtP1_NRs;means_ExtP1_NRus],[errors_ExtP1_NRs;errors_ExtP1_NRus])
        subplot(313)
        mseb(t_axis,[means_ExtP2_NRs;means_ExtP2_NRus],[errors_ExtP2_NRs;errors_ExtP2_NRus])
        
                
     
    end
    
end

            
            
            
        
        
        
    
    
    
    
    
    
    
    
     
    
