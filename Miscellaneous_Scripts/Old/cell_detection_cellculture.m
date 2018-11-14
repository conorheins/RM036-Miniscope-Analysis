%% Configure paths

fprintf('Choose Base Directory...\n')
CNMF_path = uigetdir();
cd(CNMF_path); clear CNMF_path;

cnmfe_setup;

fprintf('Choose NoRMCorre Directory...\n')
addpath(genpath(uigetdir()));

fprintf('Add any helper folders as you see fit...\n')
addpath(genpath(uigetdir()));


%% load existing data-file to analyze or add to

existing_file = 'y';
filter_data = 'n';

if strcmp(existing_file,'y')
    [fnam,fdir] = uigetfile('*.mat');
    dataFile = fullfile(fdir,fnam);
    data = matfile(dataFile,'Writable',true);
else
    data_dir = uigetdir();
    fnam = uigetfile([data_dir,filesep,'*.BTF']);
    
    all_files = dir([data_dir,filesep,'*.BTF']);
    all_files = {all_files.name};
    numFiles = length(all_files);
    
    order = zeros(numFiles,1);
    for i = 1:numFiles
        sequence_ind = regexp(all_files{i},'\d*.BTF');
        if isnan(str2double(all_files{i}(sequence_ind)))
            order(i) = 1;
        else
            order(i) = str2double(all_files{i}(sequence_ind));
        end
    end
    
    [~,srt] = sort(order);
    all_files = all_files(srt);
    
    Y = smod_bigread2(fullfile(data_dir,all_files{1}),1,1); Yfilt = zeros(size(Y),'uint16');
    dataFile = fullfile(data_dir,[all_files{1}(1:end-4),'.mat']);
    save(dataFile,'Y','Yfilt','-v7.3')
    
    tot_frames = 0;
    
    for i = 1:numFiles
        
        Y = smod_bigread2(fullfile(data_dir,all_files{i}));
        local_frames = size(Y,3);
        
        data.Y(:,:,(tot_frames+1):(tot_frames+local_frames)) = Y;
        
        tot_frames = tot_frames+local_frames;
        
        clear Y;
    end
end


if strcmp('filter_data','y')
    %% Filter memory-mapped file
    tic
    [data] = filter_memmap(data,create_filter(5,12),10000); % gaussian filter with std 5 and 12 pixel support, chunks of 10k frames
    fprintf('Time taken to filter the data: %.2f minutes\n',(toc/60))
end

Ysiz = size(data,'Y');
d1 = Ysiz(1); d2 = Ysiz(2); numFrames = Ysiz(3);


%% Threshold filtered data

patch_size = [100 100];

x_start = 1:patch_size(1):Ysiz(1);
x_end = min(x_start+patch_size(1)-1,Ysiz(1));
y_start = 1:patch_size(2):Ysiz(2);
y_end = min(y_start + patch_size(2)-1,Ysiz(2));

[X1,Y1]=meshgrid(x_start,y_start); [X2,Y2] = meshgrid(x_end,y_end);
patches = mat2cell([X1(:),X2(:),Y1(:),Y2(:)],ones(numel(X1),1),2*length(Ysiz(1:end-1)));

clear X1 X2 x_end x_start Y1 Y2 y_end y_start;

Ythresh = zeros(d1,d2,2,'uint16');
data.Ythresh = Ythresh;

options = CNMFSetParms;
sig = 5; %SNR ratio that needs to be reached for signal to be kept

tic

for i = 1:length(patches)
    
    % get patch's data
    
    Y_patch = data.Yfilt(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
    [d1p,d2p,~] = size(Y_patch);
    Y_patch = reshape(Y_patch,d1p*d2p,numFrames);
    
    % subtract each pixel's temporal median
    Y_patch = bsxfun(@minus,Y_patch,median(Y_patch,2));
        
    % get spectral estimate of noise for each pixel
    Ynoise_patch = get_noise_fft(Y_patch,options);
    
    %threshold patch with cut-off noise*sig -- anything below threshold is
    %set to 0
    Y_patch(bsxfun(@lt, Y_patch, Ynoise_patch*sig)) = 0;
    
    data.Ythresh(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),1:numFrames) = reshape(Y_patch,d1p,d2p,[]);
   
end

clear temp Y_patch Ynoise_patch Ythresh 

fprintf('Time taken to threshold filtered data: %.2f minutes',(toc/60));

%% correlation image & std image

chunkSize = 3298;
x_start = 1:chunkSize:numFrames;
x_end   = min(x_start + chunkSize - 1,numFrames);
chunks = mat2cell([x_start',x_end'],ones(length(x_start),1),2);
clear x_start x_end;

tic

Cn = zeros(d1,d2,length(chunks));
STD_image = zeros(d1,d2,length(chunks));
sz = [1 2];

for i = 1:length(chunks)
    
    chunkData = data.Ythresh(:,:,chunks{i}(1):chunks{i}(2));    
    Cn(:,:,i) = correlation_image(cast(chunkData,'double'),sz,d1,d2);
    STD_image(:,:,i) = std(cast(chunkData,'double'),[],3);
    
end
    
toc

clear chunkData

data.Cn = Cn;
data.STD_image = STD_image;


%% blob detection on statistical summary images of the standard deviation & correlation images

Cn = data.Cn;
STD_image = data.STD_image;

se = strel('disk',1); %structuring element for opening

combined = max(Cn,[],3).*max(STD_image,[],3);
thresholded = combined;
thresholded(bsxfun(@lt,combined,prctile(combined(:),97.5)))=0;
opened = imopen(thresholded,se);

temp = imregionalmax(opened); % find local maxima of opened image
ROIs = bwconncomp(temp); %find ROIs with 8-pixel connectivity on the regional maxima image
centers = zeros(ROIs.NumObjects,2);
for i = 1:ROIs.NumObjects
    [y,x] = ind2sub([d1 d2],ROIs.PixelIdxList{i}(3));
    centers(i,:) = [x, y];
end

figure;
imagesc(max(Cn,[],3)); hold on; scatter(centers(:,1),centers(:,2),'ro');

K = length(centers);
meanCn = mean(Cn,3);

%% 1st method of initialization: draw box around mean correlation image then clean up with cleanup_footprints & circular_constraints

options.medfilt_param = [3,3];
options.nrgthr = 0.8;
options.close_elem = strel('square',3);

data_thr = zeros(d1,d2,nr);
rectBounds = [12 12];

for i = 1:length(centers)
   
    r = round(centers(i,2));
    c = round(centers(i,1));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    relative_r = r - rsub(1) + 1;
    relative_c = c - csub(1) + 1;
    
    center_ind = sub2ind([length(rsub),length(csub)],relative_r,relative_c);
    
    data_temp = meanCn(rsub,csub);
    data_thr_temp = cleanup_footprints(data_temp,center_ind,options);
    data_thr(rsub,csub,i) = circular_constraints(data_thr_temp,center_ind);
    
end

A_init = reshape(data_thr,d1*d2,nr);
K = size(A_init,2);

%% 2nd method of initialization: do pengcheng-style A/C extraction on each detected cell-center

tic

A_init = [];
C_init = [];

rectBounds = [12 12];
trim_options.medfilt_param = [3,3];
trim_options.nrgthr = 0.8;
trim_options.close_elem = strel('square',3);


for i = 142:length(centers)
    
    fprintf('Extracting cell no. %d of %d\n',i,length(centers));
    
    r = round(centers(i,2));
    c = round(centers(i,1));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
       
    Y_box = data.Y(rsub,csub,:);
    Ythr_box = data.Ythresh(rsub,csub,:);
    
    Y_box = reshape(Y_box,nr*nc,[]);
    Ythr_box = reshape(Ythr_box,nr*nc,[]);
    
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
    sz = [nr, nc];
    [atemp, ci_raw, ind_success] =  extract_ac(cast(Ythr_box,'double'), cast(Y_box,'double'), ind_ctr, sz);
    temp = zeros(d1*d2,1); temp(ind_nhood) = atemp;
    temp = cleanup_footprints(reshape(temp,[d1,d2]),sub2ind([d1,d2],r,c),trim_options);
    temp = circular_constraints(temp,sub2ind([d1,d2],r,c));
    ai = reshape(temp,d1*d2,1);
    A_init = [A_init,ai];
    C_init = [C_init;ci_raw];
    
    clear Ythr_box Y_box ind_nhood ai ci_raw temp;
    
end

fprintf('Time to extract initial spatial & temporal components:   %.2f minutes\n',(toc/60))

K = size(A_init,2);


%% run CNMF on patches 

A_final = zeros(d1*d2,K);
C_final = zeros(K,numFrames);
deconv_flag = 1;
numIters = 5;
deconv_options  = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'thresholded', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', false, ... % optimize the baseline
    'optimize_smin', true);  % optimize the threshold 

trim_options.medfilt_param = [3,3];
trim_options.nrgthr = 0.8;
trim_options.close_elem = strel('square',3);
trim_options.min_pixel = 5; 

patch_size = [128 128];

x_start = 1:patch_size(1):Ysiz(1);
x_end = min(x_start+patch_size(1)-1,Ysiz(1));
y_start = 1:patch_size(2):Ysiz(2);
y_end = min(y_start + patch_size(2)-1,Ysiz(2));

[X1,Y1]=meshgrid(x_start,y_start); [X2,Y2] = meshgrid(x_end,y_end);
patches = mat2cell([X1(:),X2(:),Y1(:),Y2(:)],ones(numel(X1),1),2*length(Ysiz(1:end-1)));

clear X1 X2 x_end x_start Y1 Y2 y_end y_start;

for patch = 1:length(patches)
    
    d1p = length(patches{patch}(1):patches{patch}(2)); d2p = length(patches{patch}(3):patches{patch}(4));
    [X1,Y1] = meshgrid(patches{patch}(1):patches{patch}(2),patches{patch}(3):patches{patch}(4));
    patch_inds = sub2ind([d1,d2],Y1(:),X1(:));
    A_patch = zeros(length(patch_inds),K);
    
    neurons_in_patch = [];
    
    % find neurons in current patch whose spatial components comprise at
    % least 50% of the respective component's total energy. if less than
    % 50%, then that neuron will be initialized in a different patch.
    % Neurons whose energy is split by  > 3 patches will have to be dealt
    % with some other way...think about this
    
    for i = 1:K
        A_patch(:,i) = A_init(patch_inds,i);
        [temp,~] = sort(A_init(:,i).^2,'ascend');
        all_energy = cumsum(temp);
        patch_energy = cumsum(A_init(patch_inds,i).^2);
        if patch_energy(end) > 0.5 * all_energy(end) 
            neurons_in_patch = [neurons_in_patch,i];
        end
    end
    
    A_patch = A_patch(:,neurons_in_patch);
    C_patch = C_init(neurons_in_patch,:);
%     % center each (patched) spatial component by average of its non-zero pixels
%     Atemp = A_patch;
%     A_ind = (Atemp>0);
%     Amean = sum(Atemp)./sum(A_ind); %average of non-zero pixels
%     A_centered_patch = bsxfun(@minus, Atemp, Amean); %subtract component-wise averages
%     A_centered_patch(~A_ind) = 0; %re-set zero pixels to zero after subtraction
%     clear A_ind Amean Atemp;
    
    tic
    patch_data = double(data.Y(patches{patch}(1):patches{patch}(2),patches{patch}(3):patches{patch}(4),:));
    fprintf('Time taken to load data from patch no. %d : %.2f seconds\n',patch,toc)
    patchSz = size(patch_data);
    patch_data = reshape(patch_data,prod(patchSz(1:2)),[]);
    
%     C_patch = (A_centered_patch'*A_patch)\(A_centered_patch'*patch_data);
%     C_patch = max(C_patch,0);
    
    %% CNMF updates
    
    for iter = 1:2
        [~,C_patch,C_raw,S,kernel_pars,smin,sn ] = updateTemporal_endoscope_noObj(patch_data,A_patch,C_patch,numIters,deconv_flag,deconv_options);
        [A_patch] = updateSpatial_endoscope_noObj(patch_data,A_patch,C_patch,[],5,'hals','ellipse',patchSz);
        for i = 1:size(A_patch,2)   
            [~,center_ind] = max(A_patch(:,i));
            temp1 = cleanup_footprints(reshape(A_patch(:,i),d1p,d2p),center_ind,trim_options);
            [~,center_ind] = max(reshape(temp1,d1p*d2p,[]));
            temp2 = circular_constraints(temp1,center_ind);
            A_patch(:,i) = reshape(temp2,d1p*d2p,1);
            clear temp1 temp2 center_ind
        end
    end
    
    A_final(patch_inds,neurons_in_patch) = A_patch;
    C_final(neurons_in_patch,:) = C_patch;
    
end
    



%% merge neurons based on overlap of spatial and temporal correlations

%use intersection of minimum correlation of spatial and temporal components to determine merge
%threshold

temp = bsxfun(@times, A_final>0, 1./sqrt(sum(A_final>0))); % binarize spatial components and normalize by Frobenius energy

spatial_thr = 0.3;

A_overlap = temp'*temp; clear temp;

flag_merge = (A_overlap >= spatial_thr);

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
        imagesc(reshape(sum(A_final(:,members(1:j)),2),d1,d2));
        title(['Clique No. ', num2str(i)]);
        subplot(122)
        plot(C_final(members(j),4000:6000)); hold on;
        pause; 
    end
    close(gcf)
end

A_merge = A_final; 
C_merge = C_final; 

for i = 1:n2merge
    IDs = find(toMerge(:, i));   % IDs of neurons within this cluster
    merged_ROIs{i} = IDs;
    
    % determine searching area
    active_pixel = (sum(A_final(:,IDs), 2)>0);
    
    % update spatial/temporal components of the merged neuron
    recon = A_final(active_pixel, IDs)*C_final(IDs, :);
    ci = C_final(IDs(1), :);
    for miter=1:10
        ai = recon*ci'/(ci*ci');
        ci = ai'*recon/(ai'*ai);
    end
    
    % normalize ci
    sn = GetSn(ci);
    A_merge(active_pixel, IDs(1)) = ai*sn;
    C_merge(IDs(1), :) = ci/sn;
    
    ind_del(IDs(2:end))=true;
end

A_merge(:,ind_del) = [];
C_merge(ind_del,:) = [];

K = size(A_merge,2); % new number of neurons after merging (for subsequent steps)

clear spatial_thr temporal_thr flag_merge cliques toMerge nr n2merge merged_ROIs IDs
clear  A_overlap C_overlap active_pixel recon ai ci miter sn ind_del;

%% plotting spatial components overlaid with cell-coordinates (estimated either with maximum or with center-of-mass)

% center-of-mass calculation to find cell-coordinates
cm = com(A_merge,d1,d2); 
x = cm(:,2); y = cm(:,1);

% normalize all components and add to create full-cell mask, and overlay
% coordinates
Anorm = bsxfun(@rdivide,A_merge,max(A_merge,[],1));
imagesc(reshape(sum(Anorm,2),d1,d2));
hold on; scatter(x,y,'ro');

clear Anorm cm;


