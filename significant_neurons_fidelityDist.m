%% Script September 19th, 2018


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

%% parameters for plotting/rats to use/etc.
RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names
surround_time = [52 52]; % ~5 seconds before and after
event_nam = 'leverOUT';

display_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_09112018';

%% main-loop through rats/neurons/sessions

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    sess_nums = 1:size(cell_map,2);
    
    all_fidelity_dists = cell(length(sess_nums),1);
    
     rat_folder = sprintf('Rat%d',rat_id);
    
    if ~isdir(fullfile(display_folder,rat_folder))
        mkdir(fullfile(display_folder,rat_folder))
    end
    
    
    for kk = 1:length(sess_nums)
        
        sess_id = session_names{sess_nums(kk)};

        load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
        
        event_idx = find(strcmp(event_nam,Sess_object.event_names));
        
        Sess_object.lock2events('spikes_conv',52,52);
        Sess_object.compute_selectivity('permutation',surround_time,1000,0.05);
        Sess_object.compute_fidelity();
        
        all_fidelity_dists{kk} = Sess_object.stats_results.fidelity_matrix(Sess_object.stats_results.SigMatrix(:,event_idx) == 1,1);
    end
    
    %% plotting
    
    figure('Position',[100 100 1000 700]);
    
    for kk = 1:length(all_fidelity_dists)
        [counts,bin_centers] = hist(all_fidelity_dists{kk},10);
        plot(bin_centers,counts,'LineWidth',2.5,'DisplayName',session_names{sess_nums(kk)})
        hold on;
    end
    
    xlabel('Percentage of events of robust responding','FontSize',14)
    ylabel('Number of neurons','FontSize',14)
    legend('show');
    axis tight;
    title(sprintf('Rat %d, firing fidelity of neurons that respond to %s',rat_id,event_nam),'FontSize',16);
    
    fig_name = [rat_folder,filesep,sprintf('Rat%d_FidelityHist_10bins_%s',rat_id,event_nam),'.png'];
    saveas(gcf,fullfile(display_folder,fig_name));
    close gcf; 
    
end
    

            
            
            
            
            
            

            
            
            
            
        
        
        
        
    
    
    

    