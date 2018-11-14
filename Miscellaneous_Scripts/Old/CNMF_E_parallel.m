% Object-oriented CNMF_E on chunks of data centered on each neuron

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

%% GUI for choosing neurons

fprintf('Load correlation image...\n')
[fnam, fdir] = uigetfile('*.mat');
load(fullfile(fdir,fnam));

cMin = 0.8; %minimum correlation for thresholding
thresholded = Cn; thresholded(bsxfun(@lt,Cn,cMin))=0;

% stack of averages of thresholded matrices (first, second, and third
% chunks)
aggregate_Cn = cat(3,mean(thresholded(:,:,1:33),3),mean(thresholded(:,:,34:67),3),mean(thresholded(:,:,68:end),3));

%initialize center coordinates
center_coords = [];

for slice = 1:size(aggregate_Cn,3)
    
    figure;
    imagesc(aggregate_Cn(:,:,slice));
    
    if ~isempty(center_coords)
        for i = 1:size(center_coords,1)
            rectangle('Position',[center_coords(i,1)-2, center_coords(i,2)-2, 4, 4],'Curvature',[1,1],'LineWidth',1,'EdgeColor','r');
        end
    end
             
    xrange = get(gca, 'xlim');
    yrange = get(gca, 'ylim');
    title('click pixels out of the field to end');
    
    while true
        [tmp_x, tmp_y] = ginput(1);
        if tmp_x<xrange(1) || tmp_x>xrange(2) || (tmp_y<yrange(1)) || (tmp_y>yrange(2))
            break;
        else
            rectangle('Position',[tmp_x-2, tmp_y-2, 4, 4],'Curvature',[1,1],'LineWidth',1,'EdgeColor','r');
            center_coords = [center_coords;[tmp_x,tmp_y]];  
        end
    end
    
    sprintf('Press space-bar to advance to next slice')
    pause;
    close(gcf);
    
end

%% pengcheng method for extracting A and C

% select its neighbours for estimation of ai and ci, the box size is
%[2*rectBounds(1)+1, 2*rectBounds(2)+1]
 
writeName = ['Rat9_',fnam(1:end-4),'_initAC.mat'];

tic

A = [];
C = [];

rectBounds = [7 7];

for i = 1:length(center_coords)
    
    fprintf('Extracting cell no. %d of %d\n',i,length(center_coords));
    
    r = round(center_coords(i,2));
    c = round(center_coords(i,1));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    
    % routine to read from contiguous pixels from memory mapped file, even
    % in case of linear indices
    
    breakpts = [find(diff(ind_nhood)>1);length(ind_nhood)];
    
%     Y_box = data.Y(rsub,csub,:);
%     Ythr_box = data.Ythresh(rsub,csub,:);
    
    Ybox = zeros(nr,nc,numFrames);
    Ythr_box = zeros(nr,nc,numFrames);
    
    for col = 1:length(breakpts)
        colInds = ind_nhood((breakpts(col)-nr+1):breakpts(col));
        Ybox(:,col,:) = data.Y(colInds,:);
        Ythr_box(:,col,:) = data.Ythresh(colInds,:);
    end
    Ybox = reshape(Ybox,nr*nc,[]);
    Ythr_box = reshape(Ythr_box,nr*nc,[]);
    
    
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
    sz = [nr, nc];
    [atemp, ci_raw, ind_success] =  extract_ac(cast(Ythr_box,'double'), cast(Ybox,'double'), ind_ctr, sz);
    ai = zeros(d1*d2,1); ai(ind_nhood) = atemp;
    A = [A,ai];
    C = [C;ci_raw];
    
    clear Ythr_box; clear Y_box;
    
end

fprintf('Time to extract initial spatial & temporal components:   %.2f seconds\n',toc)
save(writeName,'A','C','center_coords');

%% merge neurons based on overlap of spatial and temporal correlations

%use intersection of minimum correlation of spatial and temporal components to determine merge
%threshold

temp = bsxfun(@times, A>0, 1./sqrt(sum(A>0))); % binarize spatial components and normalize by Frobenius energy

spatial_thr = 0.4; temporal_thr = 0.3; 

A_overlap = temp'*temp; clear temp;

C_overlap = corrcoef(C');

flag_merge = (A_overlap >= spatial_thr & C_overlap >= temporal_thr);

flag_merge = flag_merge-diag(diag(flag_merge));

%use Bron-Kerbosch algorithm to find cliques in adjacency matrix 'flag_merge'

cliques = maximalCliques(sparse(flag_merge),'v2');
toMerge = cliques(:,sum(cliques,1)>1); 

[nr, n2merge] = size(toMerge);
merged_ROIs = cell(n2merge,1);
ind_del = false(nr, 1 );    % indicator of deleting corresponding neurons

% optional visualization of found clusters
for i = 1:n2merge
    members = find(toMerge(:,i));
    figure;
    for j = 1:length(members)
        subplot(121)
        imagesc(reshape(sum(A(:,members(1:j)),2),d1,d2));
        title(['Clique No. ', num2str(i)]);
        subplot(122)
        plot(C(members(j),2000:4000)); hold on;
        pause; 
    end
    close(gcf)
end

A_init = A; 
C_init = C; 

for i = 1:n2merge
    IDs = find(toMerge(:, i));   % IDs of neurons within this cluster
    merged_ROIs{i} = IDs;
    
    % determine searching area
    active_pixel = (sum(A(:,IDs), 2)>0);
    
    % update spatial/temporal components of the merged neuron
    recon = A(active_pixel, IDs)*C(IDs, :);
    ci = C(IDs(1), :);
    for miter=1:10
        ai = recon*ci'/(ci*ci');
        ci = ai'*recon/(ai'*ai);
    end
    
    % normalize ci
    sn = GetSn(ci);
    A_init(active_pixel, IDs(1)) = ai*sn;
    C_init(IDs(1), :) = ci/sn;
    
    ind_del(IDs(2:end))=true;
end

A_init(:,ind_del) = [];
C_init(ind_del,:) = [];

K = size(A_init,2); % new number of neurons after merging (for subsequent steps)

clear center_coords spatial_thr temporal_thr flag_merge cliques toMerge nr n2merge merged_ROIs IDs
clear A C A_overlap C_overlap active_pixel recon ai ci miter sn ind_del;

%% plotting spatial components overlaid with cell-coordinates (estimated either with maximum or with center-of-mass)

% center-of-mass calculation to find cell-coordinates
cm = com(A_init,d1,d2); 
x = cm(:,2); y = cm(:,1);

% normalize all components and add to create full-cell mask, and overlay
% coordinates
Anorm = bsxfun(@rdivide,A_init,max(A_init,[],1));
imagesc(reshape(sum(Anorm,2),d1,d2));
hold on; scatter(x,y,'ro');

clear Anorm cm;

%% CNMF-E centered on chunks drawn around initialized cell-centroids

ssub = 1; tsub = 1; Fs = 10; 
neuron_full = Sources2D();
neuron_full.Fs = Fs; %frame rate

%Image dimensions and subsampling
options.d1 = d1;
options.d2 = d2;
options.ssub = ssub; %spatial downsampling factor
options.tsub = tsub; %temporal downsampling factor
options.init_method = 'greedy'; % initialization method, either done greedily or with a sparse NMF ('greedy' (default) vs. 'sparse_NMF')

% greedy_corr parameters (greedyROI_corr.m)
options.min_corr = 0.3; %minimum local correlation for initializing a neuron (default: 0.3)
options.nk       = 1;   %number of knots in B-Spline basis for detrending
options.center_psf = 1; % center (mean-subtract) the filter used to accentuate cells (Default: 1, true);
 
% greedyROI parameters (greedyROI.m)
options.gSig        = 6;             % half size of neurons to be found (default: [5,5]) -- also used to compute correlation image
options.gSiz        = 16;            % half size of bounding box for each neuron (default: 2*gSig+1) -- also used to compute correlation image
options.nb          = 1;           % number of background components (default: 1)
options.nIter       = 5;           % maximum number of rank-1 NMF iterations during refining (default: 5)
options.med_app     = 1;           % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
options.save_memory = 0;           % process data sequentially to save memory (default: 0)
options.chunkSiz    = 100;         % filter this number of timesteps each time (default: 100)
options.windowSiz   = [32,32];     % size of window over which is computed sequentially (default: 32 x 32)

% sparse_NMF parameters (sparse_NMF_initialization.m)
options.snmf_max_iter = 100;     % max # of sparse NMF iterations (default: 100)
options.err_thr       = 1e-4;    % relative change threshold for stopping sparse_NMF (default: 1e-4)
options.eta           = 1;       % frobenious norm factor: eta*max(Y(:))^2 (default: 1)
options.beta          = 0.8;     % sparsity factor (default: 0.5)

% % HALS parameters (HALS_2d.m)
options.bSiz    = 3;             % expand kernel for HALS growing (default: 3)
options.maxIter = 5;             % maximum number of HALS iterations (default: 5)

% Noise and AR coefficients calculation (preprocess_data.m)
options.noise_range    = [0.25 0.5]; % frequency range over which to estimate the noise (default: [0.25,0.5])
options.noise_method   = 'logmexp';  % method for which to estimate the noise level (default: 'logmexp')
options.flag_g         = false;      % compute global AR coefficients (default: false)
options.lags           = 5;          % number of extra lags when computing the AR coefficients (default: 5)
options.include_noise  = true;          % include early lags when computing AR coefs (default: 0)
options.pixels         = [];         % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
options.split_data     = 0;          % split data into patches for memory reasons (default: 0)
options.block_size     = [64,64];    % block size for estimating noise std in patches (default: [64,64])
options.cluster_pixels = true;       % cluster pixels to active/inactive based on the PSD density (default: true)

% UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
options.search_method  = 'ellipse';      % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'ellipse')
options.use_parallel   = 0;           % update pixels in parallel (default: 1 if present)
 
% determine_search_location.m
options.min_size = 3;                           % minimum size of ellipse axis (default: 3)
options.max_size = 6;                           % maximum size of ellipse axis (default: 8)
options.dist     = 2;                           % expansion factor of ellipse (default: 3)
options.se       = strel('disk',3,0);           % morphological element for dilation (default: strel('disk',4,0))
 
% threshold_components.m
options.nrgthr  = 0.99;                  % energy threshold (default: 0.99)
options.clos_op = strel('square',3);     % morphological element for closing (default: strel('square',3))
options.medw    = [3,3];                 % size of median filter (default: [3,3])

% UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
options.deconv_method = 'constrained_foopsi';     % method for spike deconvolution (default: 'constrained_foopsi')
options.restimate_g   = 1;                        % flag for updating the time constants for each component (default: 1)
options.temporal_iter = 2;                        % number of block-coordinate descent iterations (default: 2)
options.temporal_parallel = false;                 % flag for parallel updating of temporal components (default: true if present)

% CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
options.method       = 'cvx';   % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
options.bas_nonneg   =     1;   % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 0)
options.fudge_factor =  0.99;   % scaling constant to reduce bias in the time constant estimation (default 0.99 - no scaling)
options.resparse     =     1;   % number of times that the solution is resparsened (default: 0)

% MERGING (merge_ROIs.m)
options.merge_thr  = 0.85;      % merging threshold (default: 0.85)
options.fast_merge = 1;        % flag for using fast merging (default 1)

%  DF/F (extract_DF_F.m)
options.df_prctile = 25;          % percentile to be defined as baseline (default 50, median)
options.df_window  = 100;         % length of running window (default [], no window)

% CONTOUR PLOTS (plot_contours.m)
options.cont_threshold  = 0.9;             %the value above which a specified fraction of energy is explained (default 90%)  

% VIDEO (make_patch_video.m)
options.ind             = [1 2 3 4];       % indeces of components to be shown (deafult: 1:4)
options.skip_frame      = 1;          % skip frames when showing the video (default: 1 (no skipping))
options.sx              = 18;    % half size of representative patches (default: 16)
options.make_avi        = 0;    % flag for saving avi video (default: 0)
options.show_background = 1;    % flag for displaying the background in the denoised panel (default: 1)
options.show_contours   = 1;    % flag for showing the contour plots of the patches in the FoV (default: 0)
options.cmap            = 'default';    % colormap for plotting (default: 'default')
options.name            =['video_',datestr(now,30),'.avi'];    % name of saved video file (default: based on current date)

% PLOT COMPONENTS (view_patches.m)
options.plot_df         = 1;    % flag for displaying DF/F estimates (default: 1)
options.make_gif        = 0;    % save animation (default: 0)
options.save_avi        = 0;    % save video (default: 0)
options.pause_time      = Inf;  % time to pause between each component (default: Inf, user has to click)

% CLASSIFY COMPONENTS (classify components.m)
options.cl_thr          = 0.8;  % overlap threshold for energy for a component to be classified as true (default: 0.8)

% ORDER COMPONENTS (order_components.m)
options.nsd             = 3;    % number of standard deviations (default: 3)
options.nfr             = 5;    % number of consecutive frames (default: 5)

% parameters for microendoscope
options.min_pnr         =  10;    
options.seed_method     = 'auto';  
options.min_pixel       =  8;    % minimum number of nonzero pixels for a neuron
options.bd              =  0;    % number of pixels to be ignored in the boundary
options.deconv_flag     =  1;    % perform deconvolution or not

%Options for running deconvolution
options.deconv_options  = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold 
options.smin            = 5;     % mimimum spike size

neuron_full.options = options;


%% Create convolution kernel to model the shape of calcium transients, and update neuron object accordingly
tau_decay = 1;  
tau_rise = 0.1;
nframe_decay = ceil(10*tau_decay*neuron_full.Fs);  % number of frames in decaying period
bound_pars = false;     % bound tau_decay/tau_rise or not
neuron_full.kernel = create_kernel('exp2', [tau_decay, tau_rise]*neuron_full.Fs, nframe_decay, [], [], bound_pars); % add kernel to Sources2D parameters

%% Load and analyze the data 
sframe=1;						% user input: first frame to read (optional, default:1)
num2read=numFrames;             % user input: how many frames to read (optional, default: until the end)

K = size(A_init,2);
% how much spatial edge to have bounding neuron i's spatial component
edge_width = 3;

neuron_full.A = zeros(size(A_init));
neuron_full.C = zeros(size(C_init));
neuron_full.C = zeros(size(C_init));
neuron_full.S = zeros(size(C_init));

for i = 1:K
    
    % initialize search location for CNMF_E by starting with active pixels
    % for neuron i
    active_pix_i = find(A_init(:,i)>0);
    
    % create indices for rectangular neighborhood around neuron i
    [r, c]  = ind2sub([d1, d2], active_pix_i);
    rsub = max(1, min(r)-edge_width):min(d1, max(r)+edge_width);
    csub = max(1, min(c)-edge_width):min(d2, max(c)+edge_width);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    
    % create image-indices for pulling out data from memory-mapped data file
    [i_rows, i_cols] = ind2sub([d1,d2],ind_nhood);
    i_rows = min(i_rows):max(i_rows);
    i_cols = min(i_cols):max(i_cols);

    d1s = length(i_rows); d2s = length(i_cols);
    
    % create neuron_slave object copied from master 'neuron_full' object
    neuron_slave = neuron_full.copy();
    sframe_patch = sframe; num2read_patch = num2read;
    
    % load data and increase bit-depth to double-precision
    Y = double(data.Y(i_rows,i_cols,sframe_patch+(1:num2read_patch)-1)); 
    
    neuron_slave.options.d1 = d1s; neuron_slave.options.d2 = d2s;
    Y = neuron_slave.reshape(Y, 1);       % convert a 3D video into a 2D matrix
    
    %% Initialization of A,C -- parameters
    patch_par = [1,1]*1;
    K = 1;
    
    %% iteratively update A, C and B -- parameters
    % parameters, merge neurons
    
    % parameters, estimate the background
    spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
    thresh = 10;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
    
    bg_neuron_ratio = 1;  % spatial range / diameter of neurons
    
    % parameters, estimate the spatial components
    update_spatial_method = 'hals_thresh';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
    Nspatial = 5;       % this variable has different meanings:
    %1) udpate_spatial_method=='hals' or 'hals_thresh',
    %then Nspatial is the maximum iteration
    %2) update_spatial_method== 'nnls', it is the maximum
    %number of neurons overlapping at one pixel
    
     %% Initialize spatial and temporal components
        
    neuron_slave.initComponents_endoscope(Y, K, patch_par, [], []); % set debug_on and save_avi to false/empty
    
    %% Continue with full CNMF if neurons are found in patch
    
    [~, srt] = sort(max(neuron_slave.C, [], 2), 'descend');
    neuron_slave.orderROIs(srt);
    
    neuron_slave.options.maxIter = 3;   % iterations to update C
    
    % parameters for running iteratiosn
    nC = size(neuron_slave.C, 1);    % number of neurons
    
    maxIter = 2;        % maximum number of iterations
    miter = 1;
    while miter <= maxIter
        %% merge neurons, order neurons and delete some low quality neurons
        if miter ==1
            merge_thr = [1e-1, 0.8, .1];     % thresholds for merging neurons
            % corresponding to {sptial overlaps, temporal correlation of C,
            %temporal correlation of S}
        else
            merge_thr = [0.4, 0.5, 0.1];
        end
        
        if nC <= 1
            break
        end
        
        % merge neurons
        neuron_slave.quickMerge(merge_thr); % run neuron merges
        %sort neurons
        [~,srt] = sort(max(neuron_slave.C,[],2).*max(neuron_slave.A,[],1)','descend');
        neuron_slave.orderROIs(srt);
        
        %% udpate background (cell 1, the following three blocks can be run iteratively)
        
        Ybg = Y-neuron_slave.A*neuron_slave.C_raw;
        rr = ceil(neuron_slave.options.gSiz * bg_neuron_ratio);
        active_px = []; %(sum(IND, 2)>0); %If some missing neurons are not covered by active_px, use [] to replace IND
        [Ybg, Ybg_weights] = neuron_slave.localBG(Ybg,spatial_ds_factor,rr,active_px,neuron_slave.P.sn,thresh); %estimate local background
        
        %subtract background from the raw data to obtain signal for
        %subsequent CNMF
        Ysignal = Y - Ybg;
        
        nC = size(neuron_slave.C, 1);    % number of neurons
        
        % estimate noise
        if ~isfield(neuron_slave.P,'sn') || isempty(neuron_slave.P.sn)
            % estimate the noise for all pixels
            b0 = zeros(size(Ysignal,1), 1);
            sn = b0;
            for m = 1:size(neuron_slave.A,1)
                [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m,:));
            end
            Ysignal = bsxfun(@minus,Ysignal,b0);
        end
        
        %% update spatial & temporal components
        tic;
        for m=1:2
            %temporal
            neuron_slave.updateTemporal_endoscope(Ysignal);
            
            if nC > 1
                % merge neurons
                [merged_ROI, ~] = neuron_slave.quickMerge(merge_thr); % run neuron merges
                %sort neurons
                [~,srt] = sort(max(neuron_slave.C,[],2).*max(neuron_slave.A,[],1)','descend');
                neuron_slave.orderROIs(srt);
            end
            
            %spatial
            neuron_slave.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
            neuron_slave.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
            
            if ~exist('merged_ROI','var') || isempty(merged_ROI)
                break;
            end
        end
        fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
        
        %% pick neurons from the residual (cell 4).
        if miter==1
            
            %reset thresholds for picking neurons from residual
            min_corr = 0.8;     % minimum local correlation for a seeding pixel
            min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel
            min_pixel = 8;      % minimum number of nonzero pixels for each neuron
            bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
            neuron_slave.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
                'min_pixel', min_pixel, 'bd', bd, 'deconv_flag', true);
            
            neuron_slave.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
            [center_new, Cn_res, pnr_res] = neuron_slave.pickNeurons(Ysignal - neuron_slave.A*neuron_slave.C, patch_par, 'auto'); % method can be either 'auto' or 'manual'
        end
        
        %% stop the iteration
        temp = size(neuron_slave.C, 1);
        if or(nC==temp, miter==maxIter)
            break;
        else
            miter = miter+1;
            nC = temp;
        end
    end
    
    neuron_full.A(ind_nhood,i) = neuron_slave.A(:,overlappers==i);
    neuron_full.C(i,:) = neuron_slave.C(overlappers==i,:);
    neuron_full.C_raw(i,:) = neuron_slave.C_raw(overlappers==i,:);
    neuron_full.S(i,:) = neuron_slave.S(overlappers==i,:);
    
end
    
    
            


