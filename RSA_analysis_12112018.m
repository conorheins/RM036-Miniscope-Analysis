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

sessions2track = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

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
        distances_R_binned = 1-squareform(pdist(binned_patterns_Rsa','spearman'));
        
            
%         mean_vector_self_admin = mean(data2use(:,and(sa_filter_vector,all_trial_flags(:,2)==1)),2);
        
%         ext_ids = [];
%         Ext_names = {'Ext1','Ext2','Ext3','Ext4'};
%         for kk = 1:length(Ext_names)
%             ext_ids = [ext_ids,find(strcmp(Ext_names{kk},session_names(sess_ids)))];
%         end
        
        early_ext_ids = [];
        early_ExtNames = {'Ext1','Ext2'};
        for kk = 1:length(early_ExtNames)
            early_ext_ids = [early_ext_ids,find(strcmp(early_ExtNames{kk},session_names(sess_ids)))];
        end
        early_ext_bool = ismember(all_trial_flags(:,1),early_ext_ids);

        
        late_ext_ids = [];
        late_ExtNames = {'Ext3','Ext4'};
        for kk = 1:length(late_ExtNames)
            late_ext_ids = [late_ext_ids,find(strcmp(late_ExtNames{kk},session_names(sess_ids)))];
        end
        late_ext_bool = ismember(all_trial_flags(:,1),late_ext_ids);

        
        reinstatement_ids = [];
        Reinstatement_names = {'Reinstatement'};
        for kk = 1:length(Reinstatement_names)
            reinstatement_ids = [reinstatement_ids,find(strcmp(Reinstatement_names{kk},session_names(sess_ids)))];
        end
        
        reinst_bool = ismember(all_trial_flags(:,1),reinstatement_ids);
         
        % get pattern-to-pattern similarity for self-admin
        
        within_sa_sims = tril(distances_R_binned(1:size(binned_patterns_Rsa,2),1:size(binned_patterns_Rsa,2)),-1);
        within_sa_sims(within_sa_sims == 0) = [];
        [sa_hist_vals,sa_hist_bins] = hist(within_sa_sims,20);
        sa_hist_vals = sa_hist_vals./sum(sa_hist_vals);
        
        % deal with response trials first
        
        early_ext_patternsR = data2use(:,and(early_ext_bool,all_trial_flags(:,2)==1));
        binned_patterns_Rearly_ext = bin_patterns_overTime(early_ext_patternsR,num_trials_to_bin); 
        distances_RearlyExt_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_Rearly_ext]','spearman'));
           
        late_ext_patternsR = data2use(:,and(late_ext_bool,all_trial_flags(:,2)==1));
        binned_patterns_Rlate_ext = bin_patterns_overTime(late_ext_patternsR,num_trials_to_bin); 
        distances_RlateExt_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_Rlate_ext]','spearman'));
        
        reinst_patterns_R = data2use(:,and(reinst_bool,all_trial_flags(:,2)==1));
        binned_patterns_Rreinst = bin_patterns_overTime(reinst_patterns_R,num_trials_to_bin);  
        distances_Rreinst_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_Rreinst]','spearman'));
        
        % get similarity distributions
        
        early_ext_sims_R = distances_RearlyExt_binned((size(binned_patterns_Rsa,2)+1):size(distances_RearlyExt_binned,1),1:size(binned_patterns_Rsa,2));
        [early_ext_hist_vals_R,early_ext_hist_bins_R] = hist(early_ext_sims_R(:),20);
        early_ext_hist_vals_R = early_ext_hist_vals_R./sum(early_ext_hist_vals_R);
        
        late_ext_sims_R = distances_RlateExt_binned((size(binned_patterns_Rsa,2)+1):size(distances_RlateExt_binned,1),1:size(binned_patterns_Rsa,2));
        [late_ext_hist_vals_R,late_ext_hist_bins_R] = hist(late_ext_sims_R(:),20);
        late_ext_hist_vals_R = late_ext_hist_vals_R./sum(late_ext_hist_vals_R);
        
        reinst_sims_R = distances_Rreinst_binned((size(binned_patterns_Rsa,2)+1):size(distances_Rreinst_binned,1),1:size(binned_patterns_Rsa,2));
        [Reinst_hist_vals_R,Reinst_hist_bins_R] = hist(reinst_sims_R(:),20);
        Reinst_hist_vals_R = Reinst_hist_vals_R./sum(Reinst_hist_vals_R);
        
  
%         figure(event_i);
%         set(gcf,'Position',[100 300 800 600]);
%         
%         plot(sa_hist_bins,cumsum(sa_hist_vals),'Color',[0 0.8 1],'LineWidth',1.5,'DisplayName','Within Self-Admin  Similarity');
%         hold on; 
%         
%         plot(early_ext_hist_bins_R,cumsum(early_ext_hist_vals_R),'Color',[0 0.4 0.6],'LineWidth',1.5,'DisplayName',...
%             'Early Extinction to Self-Admin Similarity: Response Trials');
%         
%         plot(late_ext_hist_bins_R,cumsum(late_ext_hist_vals_R),'Color',[0 0.6 0.4],'LineWidth',1.5,'DisplayName',...
%             'Late Extinction to Self-Admin Similarity: Response Trials');
%         
%         plot(Reinst_hist_bins_R,cumsum(Reinst_hist_vals_R),'Color',[0.8 0.4 0.0],'LineWidth',1.5,'DisplayName',...
%             'Reinstatement to Self-Admin Similarity: Response Trials');
%         
%         legend('show')
        
       
%         plot(ext_hist_bins_R,ext_hist_vals_R,'Color',,'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Response Trials');
%         mean_sim = mean(ext_sims_R(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_R)],'--','Color',[0 0.6 0.5],'LineWidth',1,'DisplayName','ExtR Similarity Mean');
%         
%         plot(ext_hist_bins_NR,ext_hist_vals_NR,'Color',[0.6 0.2 0.1],'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Non-response Trials'); mean_sim = mean(ext_sims_R(:));
%         mean_sim = mean(ext_sims_NR(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_NR)],'--','Color',[0.6 0.2 0.1],'LineWidth',1,'DisplayName','ExtNR Similarity Mean');
% 
%         plot(Reinst_hist_bins_R,Reinst_hist_vals_R,'Color',[0.8 0.4 0.0],'LineWidth',1.5,'DisplayName','Reinstatement to Self-Admin Similarity: Response Trials'); 
%         mean_sim = mean(reinst_sims_R(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(Reinst_hist_vals_R)],'--','Color',[0.8 0.4 0.0],'LineWidth',1,'DisplayName','Reinst R Similarity Mean');
%         
%         plot(Reinst_hist_bins_NR,Reinst_hist_vals_NR,'Color',[1.0 0.7 0.0],'LineWidth',1.5,'DisplayName','Reinstatement to Self-Admin Similarity: Non-response Trials'); 
%         mean_sim = mean(reinst_sims_NR(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(Reinst_hist_vals_NR)],'--','Color',[1.0 0.7 0.0],'LineWidth',1,'DisplayName','Reinst NR Similarity Mean');
         
        
        
        % now do non-response trials
        early_ext_patternsNR = data2use(:,and(early_ext_bool,all_trial_flags(:,2)==0));
        binned_patterns_NRearly_ext = bin_patterns_overTime(early_ext_patternsNR,num_trials_to_bin);
        distances_NRearlyExt_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_NRearly_ext]','spearman'));
        
        late_ext_patternsNR = data2use(:,and(late_ext_bool,all_trial_flags(:,2)==0));
        binned_patterns_NRlate_ext = bin_patterns_overTime(late_ext_patternsNR,num_trials_to_bin); 
        distances_NRlateExt_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_NRlate_ext]','spearman'));
        
        reinst_patterns_NR = data2use(:,and(reinst_bool,all_trial_flags(:,2)==0));
        binned_patterns_NRreinst = bin_patterns_overTime(reinst_patterns_NR,num_trials_to_bin);
        distances_NRreinst_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_NRreinst]','spearman'));

       % get similarity distributions
        
        early_ext_sims_NR = distances_NRearlyExt_binned((size(binned_patterns_Rsa,2)+1):size(distances_NRearlyExt_binned,1),1:size(binned_patterns_Rsa,2));
        [early_ext_hist_vals_NR,early_ext_hist_bins_NR] = hist(early_ext_sims_NR(:),20);
        early_ext_hist_vals_NR = early_ext_hist_vals_NR./sum(early_ext_hist_vals_NR);
        
        late_ext_sims_NR = distances_NRlateExt_binned((size(binned_patterns_Rsa,2)+1):size(distances_NRlateExt_binned,1),1:size(binned_patterns_Rsa,2));
        [late_ext_hist_vals_NR,late_ext_hist_bins_NR] = hist(late_ext_sims_NR(:),20);
        late_ext_hist_vals_NR = late_ext_hist_vals_NR./sum(late_ext_hist_vals_NR);
        
        reinst_sims_NR = distances_NRreinst_binned((size(binned_patterns_Rsa,2)+1):size(distances_NRreinst_binned,1),1:size(binned_patterns_Rsa,2));
        [Reinst_hist_vals_NR,Reinst_hist_bins_NR] = hist(reinst_sims_NR(:),20);
        Reinst_hist_vals_NR = Reinst_hist_vals_NR./sum(Reinst_hist_vals_NR);
        
        bar_means = [   [mean(within_sa_sims(:));
            mean(early_ext_sims_R(:));
            mean(late_ext_sims_R(:));
            mean(reinst_sims_R(:))],...
            [NaN;
            mean(early_ext_sims_NR(:));
            mean(late_ext_sims_NR(:));
            mean(reinst_sims_NR(:))]    ];
        
        bar_errors =  [   [sem(within_sa_sims(:));
            sem(early_ext_sims_R(:));
            sem(late_ext_sims_R(:));
            sem(reinst_sims_R(:))],...
            [NaN;
            sem(early_ext_sims_NR(:));
            sem(late_ext_sims_NR(:));
            sem(reinst_sims_NR(:))]  ];
        
        
%         bar_errors =  [   [std(within_sa_sims(:))./sqrt(length(within_sa_sims(:)));
%             std(early_ext_sims_R(:))./sqrt(length(early_ext_sims_R(:)));
%             std(late_ext_sims_R(:))./sqrt(length(late_ext_sims_R(:)));
%             std(reinst_sims_R(:))./sqrt(length(reinst_sims_R(:)))],...
%             [std(within_sa_sims(:))./sqrt(length(within_sa_sims(:)));
%             std(early_ext_sims_NR(:))./sqrt(length(early_ext_sims_NR(:)));
%             std(late_ext_sims_NR(:))./sqrt(length(late_ext_sims_NR(:)));
%             std(reinst_sims_NR(:))./sqrt(length(reinst_sims_NR(:)))]    ];
%         
        figure(event_i);
        set(gcf,'Position',[100 300 800 600]);
        clr = [0 0.1 0.9;
               0.8 0.1 0];
        colormap(clr)
        barwitherr(2.5.*bar_errors,bar_means)
        xticklabels({'Training','E1-E2','E3-E4','Reinstatement'})
        
        title(sprintf('Similarity to Self-Administration Trials : locked to %s',events_of_interest{event_i}))
        ylabel('Spearman Correlation')
        
        set(gca,'FontSize',16)
        
        legend({'Response Trials','Non Response Trials'})
              
        
%         ext_filter_vector = ismember(all_trial_flags(:,1),ext_ids);
%         ext_patterns_R = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==1));
%         binned_patterns_Rext = bin_patterns_overTime(ext_patterns_R,num_trials_to_bin);
%         
%         distances_Rext_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_Rext]','spearman'));
%         
%         ext_patterns_NR = data2use(:,and(ext_filter_vector,all_trial_flags(:,2)==0));
%         binned_patterns_NRext = bin_patterns_overTime(ext_patterns_NR,num_trials_to_bin);
%         
%         distances_NRext_binned = 1-squareform(pdist([binned_patterns_Rsa,binned_patterns_NRext]','spearman'));
        
%         
%         within_sa_sims = tril(distances_R_binned(1:size(binned_patterns_Rsa,2),1:size(binned_patterns_Rsa,2)),-1);
%         within_sa_sims(within_sa_sims == 0) = [];
%         [sa_hist_vals,sa_hist_bins] = hist(within_sa_sims,20);
%         sa_hist_vals = sa_hist_vals./sum(sa_hist_vals);
%         
%         ext_sims_R = distances_Rext_binned((size(binned_patterns_Rsa,2)+1):size(distances_Rext_binned,1),1:size(binned_patterns_Rsa,2));
%         [ext_hist_vals_R,ext_hist_bins_R] = hist(ext_sims_R(:),20);
%         ext_hist_vals_R = ext_hist_vals_R./sum(ext_hist_vals_R);
%         
%         ext_sims_NR = distances_NRext_binned((size(binned_patterns_Rsa,2)+1):size(distances_NRext_binned,1), 1:size(binned_patterns_Rsa,2));
%         [ext_hist_vals_NR,ext_hist_bins_NR] = hist(ext_sims_NR(:),20);
%         ext_hist_vals_NR = ext_hist_vals_NR./sum(ext_hist_vals_NR);
%         
%         reinst_sims_R = distances_Rreinst_binned((size(binned_patterns_Rsa,2)+1):size(distances_Rreinst_binned,1),1:size(binned_patterns_Rsa,2));
%         [Reinst_hist_vals_R,Reinst_hist_bins_R] = hist(reinst_sims_R(:),20);
%         Reinst_hist_vals_R = Reinst_hist_vals_R./sum(Reinst_hist_vals_R);
%         
%         reinst_sims_NR = distances_NRreinst_binned((size(binned_patterns_Rsa,2)+1):size(distances_NRreinst_binned,1),1:size(binned_patterns_Rsa,2));
%         [Reinst_hist_vals_NR,Reinst_hist_bins_NR] = hist(reinst_sims_NR(:),20);
%         Reinst_hist_vals_NR = Reinst_hist_vals_NR./sum(Reinst_hist_vals_NR);

%         figure(event_i);
%         set(gcf,'Position',[100 300 800 600]);
%         
%         plot(sa_hist_bins,sa_hist_vals,'Color',[0 0.5 1],'LineWidth',1.5,'DisplayName','Within Self-Admin  Similarity');
%         mean_sim = mean(within_sa_sims); 
%         hold on; plot(mean_sim*ones(2,1),[0,max(sa_hist_vals)],'--','Color',[0 0.5 1],'LineWidth',1,'DisplayName','SA Similarity Mean');
%         
%         plot(ext_hist_bins_R,ext_hist_vals_R,'Color',[0 0.6 0.5],'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Response Trials');
%         mean_sim = mean(ext_sims_R(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_R)],'--','Color',[0 0.6 0.5],'LineWidth',1,'DisplayName','ExtR Similarity Mean');
%         
%         plot(ext_hist_bins_NR,ext_hist_vals_NR,'Color',[0.6 0.2 0.1],'LineWidth',1.5,'DisplayName','Extinction to Self-Admin Similarity: Non-response Trials'); mean_sim = mean(ext_sims_R(:));
%         mean_sim = mean(ext_sims_NR(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(ext_hist_vals_NR)],'--','Color',[0.6 0.2 0.1],'LineWidth',1,'DisplayName','ExtNR Similarity Mean');
% 
%         plot(Reinst_hist_bins_R,Reinst_hist_vals_R,'Color',[0.8 0.4 0.0],'LineWidth',1.5,'DisplayName','Reinstatement to Self-Admin Similarity: Response Trials'); 
%         mean_sim = mean(reinst_sims_R(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(Reinst_hist_vals_R)],'--','Color',[0.8 0.4 0.0],'LineWidth',1,'DisplayName','Reinst R Similarity Mean');
%         
%         plot(Reinst_hist_bins_NR,Reinst_hist_vals_NR,'Color',[1.0 0.7 0.0],'LineWidth',1.5,'DisplayName','Reinstatement to Self-Admin Similarity: Non-response Trials'); 
%         mean_sim = mean(reinst_sims_NR(:));
%         hold on; plot(mean_sim*ones(2,1),[0,max(Reinst_hist_vals_NR)],'--','Color',[1.0 0.7 0.0],'LineWidth',1,'DisplayName','Reinst NR Similarity Mean');
        
%         xlabel('Pattern Similarity')
%         ylabel('Probability')
%         ylim([0 0.2])
%         legend('show'); 
%         set(gca,'FontSize',16)
%         title(sprintf('Pattern Similarity, Rat %d, Event: %s',rat_id,events_of_interest{event_i}))
       
        saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12172018/',...
            sprintf('Rat%d/PatternCorr_LockedTo_%s_%dTrialBins',rat_id,events_of_interest{event_i},num_trials_to_bin)),'png');
    end
    
    %%
        

        
        
    
  
             
end
