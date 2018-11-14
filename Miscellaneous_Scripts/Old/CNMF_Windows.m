%% Full run script for RM036 (Ca2+ Imaging in PFC)

%% clear workspace
clear; clc; close all;

%% Configure paths

fprintf('Choose Base Directory...\n')
CNMF_path = uigetdir();
cd(CNMF_path); clear CNMF_path;

cnmfe_setup;

fprintf('Choose NoRMCorre Directory...\n')
addpath(genpath(uigetdir()));

fprintf('Add any helper folders as you see fit...\n')
addpath(genpath(uigetdir()));

%% load data

% point to data with .matfile
[fnam, fdir] = uigetfile('*.mat');
dataFile = fullfile(fdir,fnam);
data = matfile(dataFile,'Writable',true);

%get information about the data
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrames = Ysiz(3);    %total number of frames
frameIndices = data.sortedIndices; % trials and frames array

%% create Source2D class object for storing results and parameters
Fs = 10;             % frame rate
ssub = 1;           % spatial downsampling factor
tsub = 1;           % temporal downsampling factor
gSig = 6;           % width of the gaussian kernel, which can approximates the average neuron shape
gSiz = 16;          % maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron_full.Fs = Fs;         % frame rate

% with dendrites or not 
with_dendrites = false;
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 20;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end

%% options for running deconvolution 

deconv_flag = false;
neuron_full.options.deconv_flag = deconv_flag;

if deconv_flag
    neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
        'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
        'optimize_pars', true, ...  % optimize AR coefficients
        'optimize_b', false, ... % optimize the baseline
        'optimize_smin', true);  % optimize the threshold
end


%% loop over temporal batches and run CNMF-E on each batch

chunkSize = 500;
x_start = 1:chunkSize:numFrames;
x_end   = min(x_start + chunkSize - 1,numFrames);
chunks = mat2cell([x_start',x_end'],ones(length(x_start),1),2);
clear x_start x_end;

allChunks = cell(length(chunks),1);

neuron = neuron_full.copy();

for i = 1:length(chunks)
    
    tic;
    
    Y = double(data.Y(:, :, chunks{i}(1):chunks{i}(2)));
    [d1s,d2s, T] = size(Y);
    
    fprintf('\nThe data from Chunk %d has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', i, d1s, d2s, T, d1s*d2s*T*8/(2^30));
    fprintf('Time cost in loading data:     %.2f seconds\n', toc);
    
    Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix
    
    %% initialization of A, C
    % parameters
    if i == 1
        debug_on = false;   % visualize the initialization procedue.
        save_avi = false;   %save the initialization procedure as an avi movie.
        patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large
        K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
        
        neuron.options.seed_method = 'auto';
        min_corr = 0.65;     % minimum local correlation for a seeding pixel
        min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
        min_pixel = 5;      % minimum number of nonzero pixels for each neuron
        bd = 1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
        neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
            'min_pixel', min_pixel, 'bd', bd);
        neuron.options.nk = 1;  % number of knots for detrending
        
        
        tic;
        % greedy method for initialization
        [center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
        fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);
    
    else
        % initialize new temporal components by regressing each spatial
        % component against the data over time
        tic
        cnmfe_initialize_C; 
        fprintf('Time cost in initializing new temporal components:        %.2f seconds\n', toc);
    end
    
    %% iteratively update A, C and B
    % parameters, merge neurons
    display_merge = false;          % visually check the merged neurons
    view_neurons = false;           % view all neurons
    
    % parameters, estimate the background
    spatial_ds_factor = 4;      % spatial downsampling factor. it's for faster estimation
    thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
    
    bg_neuron_ratio = 2;  % spatial range / diameter of neurons
    
    % parameters, estimate the spatial components
    update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
    Nspatial = 5;       % this variable has different meanings:
    %1) udpate_spatial_method=='hals' or 'hals_thresh',
    %then Nspatial is the maximum iteration
    %2) update_spatial_method== 'nnls', it is the maximum
    %number of neurons overlapping at one pixel
    
    % parameters for running iterations
    nC = size(neuron.C, 1);    % number of neurons
    
    maxIter = 2;        % maximum number of iterations
    miter = 1;
    while miter <= maxIter
        %% merge neurons, order neurons and delete some low quality neurons
        if miter ==1
            merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
            % corresponding to {sptial overlaps, temporal correlation of C,
            %temporal correlation of S}
        else
            merge_thr = [0.6, 0.5, 0.1];
        end
        % merge neurons
        cnmfe_quick_merge_nosort;              % run neuron merges
        
        %% udpate background (cell 1, the following three blocks can be run iteratively)
        % estimate the background
        tic;
        cnmfe_update_BG;
        fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
        clear Ybg Ybg_weights;
        % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
        
        %% update spatial & temporal components
        tic;
        for m=1:2
            %temporal
            neuron.updateTemporal_endoscope(Ysignal);
            cnmfe_quick_merge_nosort;              % run neuron merges
            %spatial
            neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
            neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
            if isempty(merged_ROI)
                break;
            end
        end
        fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
        
        %% pick neurons from the residual (cell 4).
        if miter==1
            neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
            [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par); % method can be either 'auto' or 'manual'
        end
        
        %% stop the iteration
        temp = size(neuron.C, 1);
        if or(nC==temp, miter==maxIter)
            break;
        else
            miter = miter+1;
            nC = temp;
        end
    end
    
   allChunks{i} = neuron.copy();
   
end
    
    


