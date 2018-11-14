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

tic

A = [];
C = [];

rectBounds = [7 7];

for i = 1:length(center_coords)
    
    r = round(center_coords(i,2));
    c = round(center_coords(i,1));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    
    Y_box = dat1.Y(rsub,csub,:); Y_box = double(reshape(Y_box,nr*nc,[]));
    Ythr_box = dat2.Ythresh(rsub,csub,:); Ythr_box = double(reshape(Ythr_box,nr*nc,[]));

    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
    sz = [nr, nc];
    [atemp, ci_raw, ind_success] =  extract_ac(cast(Ythr_box,'double'), cast(Y_box,'double'), ind_ctr, sz);
    ai = zeros(d1*d2,1); ai(ind_nhood) = atemp;
    A = [A,ai];
    C = [C;ci_raw];
    
    clear Ythr_box; clear Y_box;
    
end

fprintf('Time to extract initial spatial & temporal components:   %.2f seconds\n',toc)
save('Rat9_SA2_extractAC.mat','A','C','center_coords');

%% Conor method for extracting A & C: 

% 1. Run a bunch of rank-1 NMFs on each chosen neuron center

% A = [];
% C = [];
% 
% rectBounds = [7 7];
% 
% tic
% for i = 1:length(center_coords)
%     x = round(center_coords(i,1)); y = round(center_coords(i,2));
%     currCellData = double(data.Ythresh((x-rectBounds(1)):(x+rectBounds(1)),(y-rectBounds(2)):(y+rectBounds(2)),:));
%     tempSz = size(currCellData);
%     [a_temp,c_temp] = nnmf(reshape(currCellData,prod(tempSz(1:2)),[]),1);
%     C = [C;c_temp];
%     mask = zeros(d1,d2); mask((x-rectBounds(1)):(x+rectBounds(1)),(y-rectBounds(2)):(y+rectBounds(2)),:)=reshape(a_temp',tempSz(1),tempSz(2));
%     A = [A,reshape(mask,d1*d2,[])];
%     A = [A, a_temp];
% end
% toc


% 2. Center each spatial component by average of its non-zero pixels and regress against raw data to re-initialize C

% Y = double(data.Y(:,:,1:800));
% Y = reshape(Y,prod(Ysiz(1:2)),[]);
% 
% A_expand = zeros(d1,d2,size(Atrim,2));
% for i = 1:size(Atrim,2)
%     x = round(center_coords(i,1)); y = round(center_coords(i,2));
%     A_expand((x-rectBounds(1)):(x+rectBounds(1)),(y-rectBounds(2)):(y+rectBounds(2)),i) = reshape(Atrim(:,i),maskSz(1),maskSz(2));
% end
% A_expand = reshape(A_expand,d1*d2,[]);
% A_ind = (A_expand>0); 
% Amean = sum(A_expand)./sum(A_ind); %average of non-zero pixels
% A_centered = bsxfun(@minus, A_expand, Amean); %subtract component-wise averages
% A_centered(~A_ind) = 0; %re-set zero pixels to zero after subtraction
% C = (A_centered'*A_expand)\(A_centered'*Y); %regression to solve for C
% clear A_ind Amean A_centered; 
% 
% C = C + abs(min(C,[],2)); % make non-negative by positively-shifting each
% temporal trace by its minimum value


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

%%  local background estimation followed by CNMF, given initial A & C

% Workflow:
% 1) Set parameters
% For each temporal chunk:
% 2) estimate background & save results to disk,
% 3) subtract background from raw and run CNMF on residual with A_init and
% C_init, taking turn with HALS for spatial & temporal components
% 4) trim spatial components, use spatial components as new initialization
% for next temporal chunk
% 5) run deconvolution on final C traces for each cell and run merges

%%  parameters, estimate the background
gSiz = 16;
spatial_ds_factor = 3;      % spatial downsampling factor. it's for faster estimation
thresh = 3;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)
bg_neuron_ratio = 2;  % spatial range / diameter of neurons
rr = ceil(gSiz * bg_neuron_ratio);
active_px = [];

%% options for running deconvolution 

deconv_flag = false;

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
search_method = 'pos_pix';

% trimming parameters
thr = 0.02;
sz = 2;
minpixel = 5;


%% Main loop

% determine chunk size based on RAM/SNR considerations, etc.
chunkSize = 799;
x_start = 1:chunkSize:numFrames;
x_end   = min(x_start + chunkSize - 1,numFrames);
chunks = mat2cell([x_start',x_end'],ones(length(x_start),1),2);
clear x_start x_end;

% save cell array for delete indices for each chunk
ind_small = cell(length(chunks),1);

% save background results into memory-mapped data file
data.Ybg = zeros(d1*d2,2,'uint16');

C_decon_final = zeros(size(C_init));
C_raw_final = zeros(size(C_init));

A = A_init;

tic
for win_iter = 1:length(chunks)
    
    fprintf('Running CNMF_E on temporal chunk %d of %d\n', win_iter,length(chunks));
    
    tic
    chunk = chunks{win_iter}(1):chunks{win_iter}(2);

    Y = double(data.Y(:,:,chunk)); % load data and increase bit-depth
    Y = reshape(Y,d1*d2,[]);  % vectorize
    
    A_ind = (A>0); 
    Amean = sum(A)./sum(A_ind); %average of non-zero pixels
    A_centered = bsxfun(@minus, A, Amean); %subtract component-wise averages
    A_centered(~A_ind) = 0; %re-set zero pixels to zero after subtraction
    C = (A_centered'*A)\(A_centered'*Y); %regression to solve for C
    clear A_ind Amean A_centered; 
    
    C(bsxfun(@lt,C,0))=0;
        
    Ybg = Y-A*C; % reconstruct data and subtract from raw to initialize background
    active_px = []; % If some missing neurons are not covered by active_px, use [] to replace active_px

    [Ybg, ~] = local_background(reshape(Ybg,d1,d2,[]), spatial_ds_factor, rr, active_px, [], thresh); %estimate local background
    % subtract the background from the raw data.
    Ybg = reshape(Ybg,d1*d2,[]);
    Ysignal = Y - Ybg;
 
    data.Ybg(:,chunk) = cast(Ybg,'uint16'); % reduce bit-depth of Ybg & write it to memory-mapped file 
    clear Ybg 
    
    % get noise using fft for short time series
    sn = GetSn(Ysignal);

    %% update spatial & temporal components
    
    %temporal components update --  if no deconvolution, only A, C_raw, and sn are required outputs
    [A,C_update,C_raw,S,kernel_pars,smin,sn] = updateTemporal_endoscope_noObj(Ysignal,A,C,numIters,deconv_flag,deconv_options);
    
    %spatial
    A = updateSpatial_endoscope_noObj(Ysignal,A,C_update,sn,Nspatial,update_spatial_method,search_method);
    
    % trim spatial components
    [A,temp_ind_small] = trimSpatial_noObj(A,Ysiz,thr,sz,minpixel);
    
    ind_small{win_iter} = temp_ind_small;
    C_decon_final(:,chunk) = C_update;
    C_raw_final(:,chunk) = C_raw;
    
    clear Y sn
    
    fprintf('Time taken to analyze chunk %d:    %.2f seconds\n',win_iter,toc)
    
end
fprintf('Total run time:     %.2f seconds\n',toc)




 

        
    