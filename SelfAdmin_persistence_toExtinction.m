%% Script 12.10.2018

% Goal: see if the neurons that are detected on either SA1-SA2 and a given
% extinction session, are also shared among the other extinction sessions,
% and to what degree

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

%% initialize parameters

RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

self_admin_sess = {'SA1','SA2'};
ext_ref_sessions = {'Ext1','Ext2','Ext3','Ext4'};

%% loop through rats

results_array = zeros(length(RatIDs),length(ext_ref_sessions),length(ext_ref_sessions),2); 
% 4-D array -- 1. number of rats; 2. number of reference sessions; 3.
% number of sessions to look across; 4. three versions of the persistence
% counts --  normalized by # neurons shared between self-admin and given extinction day, normalized by
% total number found in self-administration, and unnormalized counts

for ii = 1:length(RatIDs)

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
    
    % find neurons that show up in either SA1 or SA2
    
    sa_idx = [];
    for sa_i = 1:length(self_admin_sess)
        sa_idx = [sa_idx,find(strcmp(self_admin_sess{sa_i},session_names(sess_ids)))];
    end
    
    sa_neurons = find(sum(cell_map(:,sa_idx),2) > 0);
    
    % fill out results array with neurons that are shared with a given
    % extinction session, and then measure those neurons' persistence across the
    % other extinction sessions
    
    ext_idx = [];
    for ext_i = 1:length(ext_ref_sessions)
        ext_idx = [ext_idx,find(strcmp(ext_ref_sessions{ext_i},session_names(sess_ids)))];
    end
    
    for ext_i = 1:length(ext_idx)
        
        temp_id = ext_idx(ext_i);
        
        shared = find(cell_map(sa_neurons,temp_id) > 0);
        
        results_array(ii,ext_i,ext_i,1) = length(shared);
        results_array(ii,ext_i,ext_i,2) = 1;
        results_array(ii,ext_i,ext_i,3) = length(shared)/length(sa_neurons);

        all_others = ext_idx(ext_idx ~= temp_id);
                
        for j = 1:length(all_others)
            
            other_id = all_others(j);
            
            ext_j = find(ext_idx == other_id);
            
            shared2 = find(cell_map(sa_neurons(shared),other_id) > 0);
            
            results_array(ii,ext_i,ext_j,1) = length(shared2);
            results_array(ii,ext_i,ext_j,2) = length(shared2)./length(shared);
            results_array(ii,ext_i,ext_j,3) = length(shared2)./length(sa_neurons);


            
        end
        
    end
    
end

%% plotting

% normalized, i.e. percentages of total number of neurons shared between SA and given
% extinction day

mean_percentages = squeeze(mean(results_array(:,:,:,2),1));
std_percentages = squeeze(std(results_array(:,:,:,2),0,1));
% std_percentages = squeeze(std(results_array(:,:,:,2),0,1))./length(RatIDs);

colors = cool(length(ext_idx));
lineProp.col = mat2cell(colors,ones(length(ext_idx),1),3); % turn colors into cell array for use in mseb function

figure(1);
mseb(1:length(ext_idx),mean_percentages,std_percentages,lineProp);
xticks([1:length(ext_idx)]);
xticklabels(ext_ref_sessions);
set(gca,'FontSize',16);

xlabel('Extinction session')
ylabel('Percentage of cell-overlap, all common to SA1-SA2')

legend({'Relative to Ext1','Relative to Ext2','Relative to Ext3','Relative to Ext4'})
title('Among cells found in self-administration, % persistence across extinction, relative to single extinction sessions')

saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12102018/',...
    'Persistence_Plot/SA_Extinction_Persistence_normalized'),'png')


% unnormalized, i.e. raw cell counts 
mean_counts = squeeze(mean(results_array(:,:,:,1),1));
std_counts = squeeze(std(results_array(:,:,:,1),0,1));

colors = cool(length(ext_idx));
lineProp.col = mat2cell(colors,ones(length(ext_idx),1),3); % turn colors into cell array for use in mseb function

figure(3);
mseb(1:length(ext_idx),mean_counts,std_counts,lineProp);
xticks([1:length(ext_idx)]);
xticklabels(ext_ref_sessions);
set(gca,'FontSize',16);

xlabel('Extinction session')
ylabel('Numbers of persistent cells, all common to SA1-SA2')

legend({'Relative to Ext1','Relative to Ext2','Relative to Ext3','Relative to Ext4'})
title('Among cells found in self-administration, shared cell counts across extinction, relative to single extinction sessions')

saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12102018/',...
    'Persistence_Plot/SA_Extinction_Persistence_unnormalized'),'png')
            

% normalized, by total number of neurons found on self-admin
mean_sa_relative_prct = squeeze(mean(results_array(:,:,:,3),1));
std_sa_relative_prct = squeeze(std(results_array(:,:,:,3),0,1));

colors = cool(length(ext_idx));
lineProp.col = mat2cell(colors,ones(length(ext_idx),1),3); % turn colors into cell array for use in mseb function

figure(3);
mseb(1:length(ext_idx),mean_sa_relative_prct,std_sa_relative_prct,lineProp);
xticks([1:length(ext_idx)]);
xticklabels(ext_ref_sessions);
set(gca,'FontSize',16);

xlabel('Extinction session')
ylabel('Percentage of persistent cells relative to SA1-SA2')

legend({'Relative to Ext1','Relative to Ext2','Relative to Ext3','Relative to Ext4'})
title('Persistence relative to single extinction sessions, normalized by total number of SA-common cells')

saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12102018/',...
    'Persistence_Plot/SA_Extinction_Persistence_SAnormalized'),'png')

            
        
    
    

