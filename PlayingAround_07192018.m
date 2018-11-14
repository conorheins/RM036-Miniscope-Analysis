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

%% Read in data (already saved as 'subject-objects' for each session), rat by rat

RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names
sess2align = [1:3]; % indices of which sessions to look at persistence across
surround_time = [52 52]; % ~5 seconds before and after
alpha_level = 0.05; % test at the 0.05 significance level

anchor_sess = 'SA1';

%% Option 1: 
% single neurons (single-neuron trial averages, evolution over sessions)

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    neuron_ids = find(sum(cell_map(:,sess2align) > 0,2) == length(sess2align));
        
    load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,anchor_sess)));
    
    Sess_object.lock2events('c',52,52)
    Sess_object.compute_selectivity('permutation',surround_time,1000,alpha_level)
    
    sess_specific_ids = cell_map(neuron_ids,sess2align == find(strcmp(anchor_sess,session_names)));

    lever_responders = sess_specific_ids(Sess_object.stats_results.pVals(sess_specific_ids,2) < 0.05);  % column 2 is lever IN
    all_activity = Sess_object.event_locked{2}(lever_responders,:,:);
    mean_act = mean(all_activity,3);
%     sem_act = std(all_activity,0,3)./sqrt(size(all_activity,3)); %quicker, using analytic standard error 
    
    % use bootstrap samples to compute 5% - 95% confidence intervals
    bstrp_sem = bootci_array(1000,@(x)mean(x),all_activity,3);

    cross_sess_ids = neuron_ids(Sess_object.stats_results.pVals(sess_specific_ids,2) < 0.05);
    
    activity_matrix = zeros([size(mean_act),3,length(sess2align)]);
    activity_matrix(:,:,1,strcmp(anchor_sess,session_names)) = mean_act;
    activity_matrix(:,:,2,strcmp(anchor_sess,session_names)) = bstrp_sem(:,:,1);
    activity_matrix(:,:,3,strcmp(anchor_sess,session_names)) = bstrp_sem(:,:,2);
    
    for jj = 1:length(sess2align)
        
        sess_id = session_names{sess2align(jj)};
        
        if ~strcmp(sess_id,anchor_sess)
        
            load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
            
            Sess_object.lock2events('c',52,52)
            all_activity = Sess_object.event_locked{2}(cross_sess_ids,:,:);
            mean_act = mean(all_activity,3);
%             sem_act = std(all_activity,0,3)./sqrt(size(all_activity,3)); %quicker, using analytic standard error 

            % use bootstrap samples to compute 5% - 95% confidence intervals
            bstrp_sem = bootci_array(1000,@(x)mean(x),all_activity,3);

            activity_matrix(:,:,1,jj) = mean_act;
            activity_matrix(:,:,2,jj) = bstrp_sem(:,:,1);
            activity_matrix(:,:,3,jj) = bstrp_sem(:,:,2);
            
        end
           
    end
    
    for neur = 1:size(activity_matrix,1)
        
        trajectories = cell(3,size(activity_matrix,4));
        for jj = 1:length(trajectories)
            trajectories{1,jj} = activity_matrix(neur,:,1,jj);
            trajectories{2,jj} = activity_matrix(neur,:,2,jj);
            trajectories{3,jj} = activity_matrix(neur,:,3,jj);
        end
       
        
        g = gramm('x',1:105,'y',trajectories(1,:),'ymin',trajectories(2,:),'ymax',trajectories(3,:),'color',session_names(sess2align));
        g.geom_interval('geom','area');
        figure('Position',[100 100 800 450]);
        g.draw();
        
        pause; 
        close gcf;
        
    end
   
end

%% Option 2:
% population averages (within-trial population-average, with error plotted
% across trials, evolution across sessions)

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    neuron_ids = find(sum(cell_map(:,sess2align) > 0,2) == length(sess2align));
        
    load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,anchor_sess)));
    
    Sess_object.lock2events('c',52,52)
    Sess_object.compute_selectivity('permutation',surround_time,1000,alpha_level)
    
    sess_specific_ids = cell_map(neuron_ids,sess2align == find(strcmp(anchor_sess,session_names)));

    lever_responders = sess_specific_ids(Sess_object.stats_results.pVals(sess_specific_ids,2) < 0.05);  % column 2 is lever IN
    all_activity = Sess_object.event_locked{2}(lever_responders,:,:);
    mean_act = squeeze(mean(all_activity,1));
    
    % use bootstrap samples to compute 5% - 95% confidence intervals
    bstrp_sem = bootci_array(1000,@(x)mean(x),squeeze(mean_act),2);

    cross_sess_ids = neuron_ids(Sess_object.stats_results.pVals(sess_specific_ids,2) < 0.05);
    
    activity_matrix = zeros([size(mean_act,1),3,length(sess2align)]);
    activity_matrix(:,1,strcmp(anchor_sess,session_names)) = mean(mean_act,2);
    activity_matrix(:,2,strcmp(anchor_sess,session_names)) = bstrp_sem(:,1);
    activity_matrix(:,3,strcmp(anchor_sess,session_names)) = bstrp_sem(:,2);
    
    for jj = 1:length(sess2align)
        
        sess_id = session_names{sess2align(jj)};
        
        if ~strcmp(sess_id,anchor_sess)
        
            load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
            
            Sess_object.lock2events('c',52,52)
            all_activity = Sess_object.event_locked{2}(cross_sess_ids,:,:);
            mean_act = squeeze(mean(all_activity,1));

            % use bootstrap samples to compute 5% - 95% confidence intervals
            bstrp_sem = bootci_array(1000,@(x)mean(x),squeeze(mean_act),2);

            activity_matrix(:,1,jj) = mean(mean_act,2);
            activity_matrix(:,2,jj) = bstrp_sem(:,1);
            activity_matrix(:,3,jj) = bstrp_sem(:,2);
            
        end
           
    end
    
    % single sessions (trial averages across neurons, with error across trials)
    
    trajectories = cell(3,size(activity_matrix,3));
    baseline_val = min(activity_matrix(51,1,:)); %normalize all trajectories to smallest average value at middle time-step
    for jj = 1:length(trajectories)
        trajectories{1,jj} = activity_matrix(:,1,jj) - (activity_matrix(51,1,jj)- baseline_val);
        trajectories{2,jj} = activity_matrix(:,2,jj) - (activity_matrix(51,2,jj)- baseline_val);
        trajectories{3,jj} = activity_matrix(:,3,jj) - (activity_matrix(51,3,jj)- baseline_val);
    end
    
    
    g = gramm('x',1:105,'y',trajectories(1,:),'ymin',trajectories(2,:),'ymax',trajectories(3,:),'color',session_names(sess2align));
   
    g.geom_interval('geom','area');
    figure('Position',[100 100 800 450]);
    g.draw();
            
 
end
    
    
        
        

        
    
    
    
   
                 
end