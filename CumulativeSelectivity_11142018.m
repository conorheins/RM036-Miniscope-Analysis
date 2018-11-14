%% Script 11.14.2018
% look at cumulative neuron map of all neurons (namely, across all sessions) responsive to event of
% interest

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

event_name = 'leverOUT'; 
data_type = 'spikes_conv'; 

%%  main-loop through rats/neurons/sessions

ii = 1; % for debugging purposes

% for ii = 1:length(RatIDs)

rat_id = RatIDs(ii);

cell_map = load_cell_map(registration_folder,rat_id);

rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));
exist_sessions = {};
all_data_files = dir([rat_folder,filesep,'*.mat']);
all_data_files = {all_data_files(:).name};
for dat_i = 1:length(all_data_files)
    num_remove = length(sprintf('Rat%d_',rat_id))+1;
    exist_sessions{dat_i} = all_data_files{dat_i}(num_remove:end-4);
end

sess_ids = [];
for kk = 1:length(session_names)
    if sum(strcmp(exist_sessions,session_names{kk})) > 0
        sess_ids = [sess_ids,kk];
    end
end

kk = 1; % for debugging purposes

for kk = 1:length(sess_ids)
             
        sess_name = session_names{sess_ids(kk)};
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        Sess_object.lock2events(data_type,52,52);
        Sess_object.compute_selectivity('permutation',[52 52],1000,0.05);
        
        sig_responders = Sess_object.stats_matrix
        
        

end    
    
    
 
% end
    