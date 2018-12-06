%% Full run script for RM036 (Ca2+ Imaging in PFC)
% Conor Heins, 10.22.2017, edited/cleaned up 03.06.2018

%% clear workspace & configure paths
clear; clc; close all;

fprintf('Choose Base Directory...\n')
master_directory = uigetdir();
cd(master_directory); clear master_directory;

setup_paths;

%% convert image sequence to writable matfile and (optionally) filter the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNCOMMENT CODE BELOW IF YOU WANT THE PROGRAM TO ASK YOU FOR YOUR DESIRED
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

already_loaded = input('Would you like to load an existing .mat file? (y/n, default n)\n','s');
if isempty(already_loaded)
    already_loaded = 'n';
end

trial_based = input('Are you analyzing trial-based data (y/n, default y)\n','s');
if isempty(trial_based)
    trial_based = 'y';
end

if strcmp(trial_based,'y')
    trial_length = input('What is the length of a single trial? (default, 800 frames)');
    if isempty(trial_length)
        trial_length = 800; % default for Tarun's experiment RM036, with 80-second trials
    end
else 
    trial_length = [];
end

filter_flag = input('Would you like to filter the data (for motion-correction, for instance)? (y/n, default y)\n','s');
if isempty(filter_flag)
    filter_flag = 'y';
end

run_MC = input('Would you like to motion-correct the data? (y/n, default y)\n','s');
if isempty(run_MC)
    run_MC = 'y';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNCOMMENT CODE BELOW IF YOU WANT TO MANUALLY PRE-SET OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% already_loaded = 'n';
% trial_based = 'y';
% trial_length = 800;
% filter_flag = 'y';
% run_MC = 'y';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRAME DELETION & FILTER PARAMETERS (SHOULDN'T HAVE TO CHANGE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delFrames = 0; % integer or list of integers that give indices of frames to throw out (if trial-based, this index relative to the start of every trial)

gaussFilt = create_filter(6,16); %function that creates gaussian filter for image filtering step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND MEMORY MAP DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

memory_map_data;

%% Run motion correction on filtered data and apply shifts to original (NoRMCorre)

if strcmp(run_MC,'y')
    
    
    MC_fnam = fullfile(fdir,[fnam(1:end-4),'_MCfilt.mat']);
    
    options_mc = NoRMCorreSetParms('d1',d1,'d2',d2,...
        'bin_width',50,'max_shift',20,'iter',1,...
        'use_parallel',true,'memmap',true,...
        'mem_filename',MC_fnam,...
        'output_type','memmap');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % register data and apply shifts to removed percentile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tic; [M1,shifts1,template1] = normcorre_batch_onFilt(data,options_mc); 
    fprintf('Time taken to motion-correct data: %.2f minutes\n',(toc/60)) % register filtered data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply shifts to original (unfiltered) dataset and overwrite the data with registered version    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    
    options_mc.memmap = false;
    options_mc.output_type = 'mat';
    
    chunkSize = 40000; %% change to 40,000 when running a full data-set (i.e. with 80k frames)
    
    data = apply_shifts_wrapper(data,shifts1,options_mc,chunkSize,Ysiz);
      
    fprintf('Time taken to apply shifts to original dataset: %.2f minutes\n',(toc/60))

end

%% Pre-process data 

options_preprocess = CNMFSetParms;
sig = 5;

preprocess_threshold;

%% correlation and STD images

tic
chunkSize = 7990; 
[Cn,STD_image] = corr_STD_imgs_chunks(data,chunkSize,Ysiz);  
fprintf('Time taken to generate correlation/STD image stack: %.2f minutes\n',(toc/60))

data.Cn = Cn;
data.STD_image = STD_image;
save(fullfile(fdir,[fnam(1:end-4),'_CorrImage.mat']),'Cn','STD_image');

%% seed ROI locations based on local maxima in the correlation and STD image stacks and do quick distance-based merge operation 

tic

options_morph.se = strel('disk',2);
options_morph.prctiles = [90 90];
display_flag = false; % set this to true if you want to click through different correlation images

allROIs = initialize_centers(Cn,STD_image,options_morph,display_flag);

dmin = 2; display_flag = false; % set this to true if you want to view initialized cell map
[allROIs] = distmerge(allROIs,dmin,display_flag,mean(Cn.*STD_image,3));

% Initialization of spatial components: draw box around centroids then clean up with morphological operations

options_init.medfilt_param = [3,3];
options_init.nrgthr = 0.5;
options_init.close_elem = strel('square',3);
rectBounds = [12 12];

[A_init,numROIs] = initialize_A(allROIs,mean(Cn.*STD_image,3),rectBounds,options_init);

fprintf('Time taken to seed centers, initialize spatial components, and threshold: %.2f minutes\n',(toc/60))


%% extract calcium traces

tic
[A,C,cell_coords] = estimate_AandC(data,A_init,rectBounds,Ysiz);
data.A = A; data.C = C;
fprintf('Time taken to extract cell shapes & calcium traces: %.2f minutes\n',(toc/60))

    
    
    
    
    


    
    




    
