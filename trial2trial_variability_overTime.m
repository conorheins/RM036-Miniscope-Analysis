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

RatIDs = [8,9,10,11,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names
sessions2analyze = [1:7]; % indices of which sessions to look at persistence across
surround_time = [52 52]; % ~5 seconds before and after
alpha_level = 0.05; % test at the 0.05 significance level

all_rats_fidelities = cell(1,length(RatIDs));

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    within_rat_fidelities = cell(1,length(sessions2analyze));
        
    for jj = 1:length(sessions2analyze)
        
        sess_id = session_names{sessions2analyze(jj)};
        
        load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
        
        Sess_object.lock2events('spikes_denoised',52,52)
        Sess_object.compute_selectivity('permutation',surround_time,1000,alpha_level)
        Sess_object.compute_fidelity();
        
        within_rat_fidelities{jj} = Sess_object.stats_results.fidelity_matrix(:,2); % second column is lever extension
        
    end
    
    all_rats_fidelities{ii} = within_rat_fidelities;
                 
end

%% plots of results

fidelity_stats_matrix = zeros(length(RatIDs),length(sessions2analyze),2); % 3rd dimension carries actual numbers (mean in first slice, SD in second)

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    for jj = 1:length(sessions2analyze)
        
        sess_id = session_names{sessions2analyze(jj)};
        
        fidelity_stats_matrix(ii,jj,1) = mean(all_rats_fidelities{ii}{jj}(~isnan(all_rats_fidelities{ii}{jj})));
        fidelity_stats_matrix(ii,jj,2) = std(all_rats_fidelities{ii}{jj}(~isnan(all_rats_fidelities{ii}{jj})));

    end
    
end

figure;
for ii = 1:length(RatIDs)
    temp_line = plot(squeeze(fidelity_stats_matrix(ii,:,1)),'LineWidth',2);
    eval(sprintf('h%d = temp_line',ii));
    hold on;

    plot(squeeze(fidelity_stats_matrix(ii,:,1)) + squeeze(fidelity_stats_matrix(ii,:,2)),'Color',temp_line.Color,'LineWidth',0.25);
    plot(squeeze(fidelity_stats_matrix(ii,:,1)) - squeeze(fidelity_stats_matrix(ii,:,2)),'Color',temp_line.Color,'LineWidth',0.25);
end
    
    



