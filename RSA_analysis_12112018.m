%% Script 12.11.2018

% Goal: concatenate activity vectors across multiple bins and look at the
% similarity between adjacent activity in those bins across time, for
% neurons that persist across those bins, and also separate by
% response/non-response trials


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

sessions2track = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4'};

events_of_interest = {'HLON','leverOUT','leverIN','pellets'}; 

temporal_bins = {[-5,5],[-5,5],[-5,5],[-5,5]};
Fs = 10.49; % sampling rate, in frames/second


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
    
    % find neurons that show up in the intersection of Self-Admin through
    % end of Extinction
    
    
    sess2track_ids = [];
    for kk = 1:length(sessions2track)
        sess2track_ids = [sess2track_ids,find(strcmp(sessions2track{kk},session_names(sess_ids)))];
    end
    
    absolute_idx_neurons = find(sum(cell_map(:,sess2track_ids)>0,2) == length(sess2track_ids));
    
    accum_data = cell(length(events_of_interest),1);

    
    %% accumulate data from appropriate sessions
    
    all_trial_flags = [];
  
    for kk = 1:length(sess2track_ids)
        
        sess_name = session_names{sess_ids(sess2track_ids(kk))};
        
        sess_specific_idx = cell_map(absolute_idx_neurons,sess_ids(sess2track_ids(kk)));
        
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        reshaped_events = reshape(full(Sess_object.event_matrix),799,Sess_object.num_trials,size(Sess_object.event_matrix,2));
        
        R_trial_bool = logical(sum(squeeze(reshaped_events(:,:,strcmp('cueON',Sess_object.event_names))),1));
        
        all_trial_flags = [all_trial_flags;[sess_ids(sess2track_ids(kk))*ones(length(R_trial_bool),1), R_trial_bool']];
                
        reshaped_neur = reshape(Sess_object.(data_type)(sess_specific_idx,1:(799*Sess_object.num_trials)),[length(sess_specific_idx), 799,Sess_object.num_trials]);

        for event_i = 1:length(events_of_interest)
            
            [event_tmsp,trial_id] = find(squeeze(reshaped_events(:,:,strcmp(events_of_interest{event_i},Sess_object.event_names))));
            
            if ~isempty(event_tmsp)
                [trial_id,temp] = unique(trial_id,'first');
                event_tmsp = event_tmsp(temp);
                event_i_data = zeros(length(sess_specific_idx),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(trial_id));
                
                for jj = 1:length(event_tmsp)
                    temp = event_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data(:,:,jj) = reshaped_neur(:,win_edges(1):win_edges(2),trial_id(jj));
                end
                
                accum_data{event_i} = cat(3,accum_data{event_i},event_i_data);
                
            end
                
        end
    end
    
    %% look at decay of pattern similarity from self-administration average vector 
    
    time_win = 1:105;
    
    for event_i = 1:length(events_of_interest)
        
        all_event_i_data = squeeze(prctile(accum_data{event_i}(:,time_win,:),90,2));
        tf_idf_d = calcTFIDF(all_event_i_data);
        
        data2use = all_event_i_data;
        
        sa_ids = [];
        SA_names = {'SA1','SA2'};
        for kk = 1:length(SA_names)
            sa_ids = [sa_ids,find(strcmp(SA_names{kk},session_names(sess_ids)))];
        end      
        
        sa_filter_vector = ismember(all_trial_flags(:,1),sa_ids);
        mean_vector_self_admin = mean(data2use(:,and(sa_filter_vector,all_trial_flags(:,2)==1)),2);
        
        ext_ids = [];
        Ext_names = {'Ext1','Ext2','Ext3','Ext4'};
        for kk = 1:length(Ext_names)
            ext_ids = [ext_ids,find(strcmp(Ext_names{kk},session_names(sess_ids)))];
        end
        
        ext_filter_vector = ismember(all_trial_flags(:,1),ext_ids);
        ext_vectors_R = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==1));
        ext_vectors_NR = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==0));
        
        distances_R = 1-squareform(pdist([mean_vector_self_admin,ext_vectors_R]','spearman'));
        distances_NR = 1-squareform(pdist([mean_vector_self_admin,ext_vectors_NR]','spearman'));
        
        plot(distances_R(1,2:104));
        hold on; plot(distances_NR(1,2:104));
        
        pause; close gcf
        
    end
    
    %%
        

        
        
    
  
             
end
