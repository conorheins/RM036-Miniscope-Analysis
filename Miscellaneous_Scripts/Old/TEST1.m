% CNMF_E on small chunks of data centered on neural coordinates

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

%% merge neurons based on overlap of spatial and temporal correlations

%use intersection of minimum correlation of spatial and temporal components to determine merge
%threshold

temp = bsxfun(@times, A>0, 1./sqrt(sum(A>0))); % binarize spatial components and normalize by Frobenius energy

spatial_thr = 0.3; temporal_thr = 0.3; 

A_overlap = temp'*temp; clear temp;

C_overlap = corrcoef(C');

flag_merge = (A_overlap >= spatial_thr & C_overlap >= temporal_thr);

% Old method: Using inter-centroid distance to determine merging criterion

% find max pixel indices for each trimmed spatial component

% [~, tempInd] = max(A);
% [y,x] = ind2sub([d1,d2],tempInd);
% new_coords = [x',y']; clear x y;
% dmin = 5;
% distMatrix = sqrt(bsxfun(@minus, new_coords(:,1), new_coords(:,1)').^2 + bsxfun(@minus, new_coords(:,2), new_coords(:,2)').^2);
% flag_merge = (distMatrix <= dmin);

flag_merge = flag_merge-eye(size(A,2));

%use Bron-Kerbosch algorithm to find cliques in adjacency matrix 'flag_merge'

cliques = maximalCliques(sparse(flag_merge),'v2');
toMerge = cliques(:,sum(cliques,1)>1); 

[nr, n2merge] = size(toMerge);
merged_ROIs = cell(n2merge,1);
ind_del = false(nr, 1 );    % indicator of deleting corresponding neurons

% optional visualization of found clusters
% for i = 1:n2merge
%     members = find(toMerge(:,i));
%     figure;
%     for j = 1:length(members)
%         subplot(121)
%         imagesc(reshape(sum(A(:,members(1:j)),2),d1,d2));
%         title(['Clique No. ', num2str(i)]);
%         subplot(122)
%         plot(C(members(j),2000:4000)); hold on;
%         pause; 
%     end
%     close(gcf)
% end

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

% [~, tempInd] = max(A_final);
% [y,x] = ind2sub([d1,d2],tempInd);
% new_coords = [x',y']; clear x y;

% center-of-mass calculation to find cell-coordinates
cm = com(A_init,d1,d2); 
x = cm(:,2); y = cm(:,1);

% normalize all components and add to create full-cell mask, and overlay
% coordinates
Anorm = bsxfun(@rdivide,A_init,max(A_init,[],1));
imagesc(reshape(sum(Anorm,2),d1,d2));
% hold on; scatter(new_coords(:,1),new_coords(:,2),'ro');
hold on; scatter(x,y,'ro');

clear Anorm cm;

% optional drawing of a rectangle around cells that are putative merges
% that the merging algorithm missed

% box = getrect(gcf);
% tempInd = find(x>box(1) & x < box(1)+box(3) ...
%         & y > box(2) & y < box(2)+box(4));

%%  parameters, estimate the background
gSiz = 16;
spatial_ds_factor = 1;      % spatial downsampling factor. it's for faster estimation
thresh = 3;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
bg_neuron_ratio = 2;  % spatial range / diameter of neurons
rr = ceil(gSiz * bg_neuron_ratio);
active_px = [];

%% options for running deconvolution 

deconv_flag = true;

deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'constrained', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold

%% parameters for temporal and spatial updates

numIters = 5; %number of iterations for temporal HALS

update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 5;       % this variable has different meanings:
%1) udpate_spatial_method=='hals' or 'hals_thresh',
%then Nspatial is the maximum iteration
%2) update_spatial_method== 'nnls', it is the maximum
%number of neurons overlapping at one pixel

% search method for 'determine_search_location' within
% 'updateSpatial_endoscope_noObj'
search_method = 'ellipse';

% trimming parameters
thr = 0.02;
sz = 2;
minpixel = 5;

% how much spatial edge to have bounding neuron i's spatial component
edge_width = 3;

%% Main loop -- loop through K CNMF_Es centered on each neuron

% save background results into memory-mapped data file
% data.Ybg = zeros(d1*d2,2,'uint16');

C_decon_final = zeros(size(C_init));
C_raw_final = zeros(size(C_init));
S_final = zeros(size(C_init));

A_final = A_init;
C_final = zeros(size(C_init));

ind_small = cell(K,1);

tic
for i = 1:K
    
%     tic
    
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

    % find overlapping neurons in the neighborhood and initialize 
    % low-rank matrix A with N + 1 columns (where N is number of
    % neurons overlapping with neuron i)
    A_temp = A_init(ind_nhood,:);
    overlappers = find(sum(A_temp)>0);
    A = A_temp(:,overlappers);
    i_ind = find(overlappers==i);
    
    % initialize C with extracted temporal components of the overlapping
    % neurons
    C = C_init(overlappers,:);
    
    % create image-indices for pulling out data from memory-mapped data file
    [i_rows, i_cols] = ind2sub([d1,d2],ind_nhood);
    i_rows = min(i_rows):max(i_rows);
    i_cols = min(i_cols):max(i_cols);

    d1s = length(i_rows); d2s = length(i_cols);
    
    % load data and increase bit-depth
    Y = double(data.Y(i_rows,i_cols,:)); 
    Y = reshape(Y,d1s*d2s,numFrames);

%%     Regression to solve for C
    A_ind = (A>0); 
    Amean = sum(A)./sum(A_ind); %average of non-zero pixels
    A_centered = bsxfun(@minus, A, Amean); %subtract component-wise averages
    A_centered(~A_ind) = 0; %re-set zero pixels to zero after subtraction
    C = (A_centered'*A)\(A_centered'*Y); %regression to solve for C
    clear A_ind Amean A_centered; 
    
    C(bsxfun(@lt,C,0))=0;
    
    C_final(i,:) = C(overlappers==i,:);
    
end

toc
    
    

%% BG Update

%     Ybg = Y-A*C; % reconstruct data and subtract from raw to initialize local background
%     active_px = []; % If some missing neurons are not covered by active_px, use [] to replace active_px
% 
%     [Ybg, ~] = local_background(reshape(Ybg,d1s,d2s,[]), spatial_ds_factor, rr, active_px, [], thresh); %estimate local background
%     subtract the background from the raw data.
%     Ybg = reshape(Ybg,d1s*d2s,[]);
%     Ysignal = Y - Ybg;
%     
%     data.Ybg(ind_nhood,:) = cast(Ybg,'uint16'); % reduce bit-depth of Ybg & write it to memory-mapped file 
%     clear Ybg 
%     
%     get noise using fft for short time series
%     sn = GetSn(Ysignal);
%     
%      % estimate the noise for all pixels
%     b0 =zeros(size(Ysignal,1), 1);
%     sn = b0;
%     parfor m=1:size(A,1)
%         [b0(m), sn(m)] = estimate_baseline_noise(Ysignal(m, :));
%     end     
%     Ysignal = bsxfun(@minus, Ysignal, b0);
%     % update spatial & temporal components
%     
%     temporal components update --  if no deconvolution, only A, C_raw, and sn are required outputs
%     [A,C,C_raw,S,kernel_pars,smin,sn] = updateTemporal_endoscope_noObj(Ysignal,A,C,numIters,deconv_flag,deconv_options);
%     
%     spatial
%     spatialSz = [d1s;d2s];
%     A = updateSpatial_endoscope_noObj(Ysignal,A,C,sn,Nspatial,update_spatial_method,search_method,spatialSz);
%     
%     trim spatial components after HALS update
%     [A,temp_ind_small] = trimSpatial_noObj(A,spatialSz,thr,sz,minpixel);
%     
%     put cell i's spatial component into final A matrix
%     A_final(ind_nhood,i) = A(:,overlappers==i);
%     
%     ind_small{i} = overlappers(temp_ind_small);
%     
%     C_decon_final(i,:) = C(overlappers==i,:);
%     C_raw_final(i,:) = C_raw(overlappers==i,:);
%     S_final(i,:) = S(overlappers==i,:);
%     
%     clear Y Ysignal Ybg sn
%     
%     fprintf('Time taken to analyze neuron %d:    %.2f seconds\n',i,toc)
    
end
fprintf('Total run time:     %.2f seconds\n',toc)

%% OLD METHOD FOR FINDING OVERLAPPING PIXELS -- SAVE FOR LATER MAYBE

% temp = real(A_init>0); % binarize spatial components
% A_overlap = temp'*temp; clear temp;
% A_overlap = A_overlap - diag(diag(A_overlap));
% [cells_1, cells_2] = ind2sub([K,K],find(A_overlap>0));

    
%     if ismember(i,cells_1)
%         overlappers = cells_2(cells_1==i);
%         
%         A_sub = [A_init(ind_nhood,i), zeros(length(ind_nhood,length(overlappers));
%         
%         for j = 1:length(overlappers)
%             active_pix_j = find(A_init(:,overlappers(j))>0);
%             
% %             [j_cols, j_rows] = ind2sub([d1,d2],active_pix_j);
%             
%             common_pix = intersect(ind_nhood,active_pix_j);
%             
% %             [comm_cols, comm_rows] = ind2sub([d1,d2],common_pix);
        
        
        
    
    
    
    