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

events_of_interest = {'HLON','leverOUT','leverIN'}; 

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
        
        clear Sess_object
    end
    
    %% look at decay of pattern similarity from self-administration average vector 
    
    time_win = 45:70;
    
    num_trials_to_bin = 3;
    
    for event_i = 1:length(events_of_interest)
        
        %%
        
        close gcf;

        all_event_i_data = squeeze(mean(accum_data{event_i}(:,time_win,:),2));
            
%         tf_idf_d = calcTFIDF(all_event_i_data);
        
        data2use = all_event_i_data;
        
        sa_ids = [];
        SA_names = {'SA1','SA2'};
        for kk = 1:length(SA_names)
            sa_ids = [sa_ids,find(strcmp(SA_names{kk},session_names(sess_ids)))];
        end      
        
        sa_filter_vector = ismember(all_trial_flags(:,1),sa_ids);
        
        % get activity averages across groups of num_trials_to_bin trials
        sa_patterns_R = data2use(:,and(sa_filter_vector,all_trial_flags(:,2)==1));
        
        binned_patterns_Rsa = bin_patterns_overTime(sa_patterns_R,num_trials_to_bin);
            
%         mean_vector_self_admin = mean(data2use(:,and(sa_filter_vector,all_trial_flags(:,2)==1)),2);
        
        ext_ids = [];
        Ext_names = {'Ext1','Ext2','Ext3','Ext4'};
        for kk = 1:length(Ext_names)
            ext_ids = [ext_ids,find(strcmp(Ext_names{kk},session_names(sess_ids)))];
        end

        ext_filter_vector = ismember(all_trial_flags(:,1),ext_ids);
        ext_patterns_R = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==1));
        binned_patterns_Rext = bin_patterns_overTime(ext_patterns_R,num_trials_to_bin);
        
        distances_R_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_Rext]','spearman'));
        
        ext_patterns_NR = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==0));
        binned_patterns_NRext = bin_patterns_overTime(ext_patterns_NR,num_trials_to_bin);
        
        distances_NR_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_NRext]','spearman'));
               
        within_sa_sims = tril(distances_R_binned(1:size(binned_patterns_Rsa,2),1:size(binned_patterns_Rsa,2)),-1);
        within_sa_sims(within_sa_sims == 0) = [];
        [sa_hist_vals,sa_hist_bins] = hist(within_sa_sims,20);
        sa_hist_vals = sa_hist_vals./sum(sa_hist_vals);
        
        ext_sims_R = distances_R_binned((size(binned_patterns_Rsa,2)+1):size(distances_R_binned,1),1:size(binned_patterns_Rsa,2));
        [ext_hist_vals_R,ext_hist_bins_R] = hist(ext_sims_R(:),20);
        ext_hist_vals_R = ext_hist_vals_R./sum(ext_hist_vals_R);
        
        ext_sims_NR = distances_NR_binned((size(binned_patterns_Rsa,2)+1):size(distances_R_binned,1), 1:size(binned_patterns_Rsa,2));
        [ext_hist_vals_NR,ext_hist_bins_NR] = hist(ext_sims_NR(:),20);
        ext_hist_vals_NR = ext_hist_vals_NR./sum(ext_hist_vals_NR);
        

        figure(event_i);
        set(gcf,'Position',[100 300 800 600]);
        
        plot(sa_hist_bins,sa_hist_vals,'Color',[0 0.5 1],'LineWidth',1.5,'DisplayName','Within Self-Admin  Similarity');
        mean_sim = mean(within_sa_sims); 
        hold on; plot(mean_sim*ones(2,1),[0,max(sa_hist_vals)],'--','Color',[0 0.5 1],'LineWidth',1,'DisplayName','SA Similarity Mean');
        
        plot(ext_hist_bins_R,ext_hist_vals_R,'Color',[0 0.6 0.5],'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Response Trials');
        mean_sim = mean(ext_sims_R(:));
        hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_R)],'--','Color',[0 0.6 0.5],'LineWidth',1,'DisplayName','ExtR Similarity Mean');
        
        plot(ext_hist_bins_NR,ext_hist_vals_NR,'Color',[0.6 0.2 0.1],'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Non-response Trials'); mean_sim = mean(ext_sims_R(:));
        mean_sim = mean(ext_sims_NR(:));
        hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_NR)],'--','Color',[0.6 0.2 0.1],'LineWidth',1,'DisplayName','ExtNR Similarity Mean');
        
        xlabel('Pattern Similarity')
        ylabel('Probability')
        ylim([0 0.2])
        legend('show'); 
        set(gca,'FontSize',16)
        title(sprintf('Pattern Similarity, Rat %d, Event: %s',rat_id,events_of_interest{event_i}))
      
        saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12122018/',...
            sprintf('Rat%d/PatternCorr_LockedTo_%s_%dTrialBins',rat_id,events_of_interest{event_i},num_trials_to_bin)),'png');
    end
    
    %%
        

        
        
    
  
             
end
