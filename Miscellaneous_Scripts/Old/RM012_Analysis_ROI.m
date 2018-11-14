%% Full run script for RM036 (Ca2+ Imaging in PFC)

% last edit, 08202017, incorporated reshaped matrices (2D instead of 3D) into memory-mapping & motion-correcting

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

%% write image sequence to matfile and optionally filter the data

% set options
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
        trial_length = 800;
    end
else 
    trial_length = [];
end

filter_flag = input('Would you like to filter the data (for motion-correction, for instance)? (y/n, default n)\n','s');
if isempty(filter_flag)
    filter_flag = 'n';
end

delFrames = 0; % integer or list of integers that give indices of frames to throw out (if trial-based, this index relative to the start of every trial)

gaussFilt = create_filter(6,16); %function that creates gaussian filter for image filtering step

% load data

if strcmp(already_loaded,'y')
    fprintf('Choose an existing .mat file to analyze or add data to...\n');
    [fnam, fdir] = uigetfile('*.mat');
    dataFile = fullfile(fdir,fnam);
    data = matfile(dataFile,'Writable',true);
    add_data = input(sprintf('Would you like to add more data to existing file %s ? (y/n, default n)\n',fnam),'s');
    if strcmp(add_data,'y')
        add_data = 1;
        [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,add_data,trial_based,trial_length,delFrames);
        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origIndices = data.sortedIndices;
        sortedIndices_new = [origIndices;frameIndices];
        data.sortedIndices = sortedIndices_new; clear origIndices sortedIndices_new;
    else
        add_data = 0;
        if strcmp(filter_flag,'y')
            chunkSize = 10000;
            data = filter_memmap(data,gaussFilt,chunkSize);
        end       
    end
elseif strcmp(already_loaded,'n')
    
    fprintf('Choose a folder to determine name of data file...\n');
    dataFile = [uigetdir(),'.mat'];
    [fdir, fnam, temp] = fileparts(dataFile);
    fnam = strcat(fnam,temp); clear temp;
    
    add_data = input(sprintf('Would you like to add more data after loading new file %s ? (y/n, default n)\n',fnam),'s');
    
    Y = [];
    if strcmp(filter_flag,'y')
        Yfilt = []; save(dataFile,'Y','Yfilt','-v7.3')
    else
        save(dataFile,'Y','-v7.3');
    end
    
    data = matfile(dataFile,'Writable',true);
    
    [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,0,trial_based,trial_length,delFrames);
    data.Ysiz = Ysiz; data.sortedIndices = frameIndices;
    
    
    if strcmp(add_data,'y')
        add_data = 1;
        [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,add_data,trial_based,trial_length,delFrames);
        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origIndices = data.sortedIndices;
        sortedIndices_new = [origIndices;frameIndices];
        data.sortedIndices = sortedIndices_new; clear origIndices sortedIndices_new;
    else
        add_data = 0; 
    end
end

%get information about the data
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrames = Ysiz(3);    %total number of frames
frameIndices = data.sortedIndices; % trials and frames array 

if strcmp(trial_based,'y')
    fprintf('\nThe data has been mapped to a hard disk. It has %d X %d pixels X %d frames, and is comprised of %d total trials. \nLoading all data requires %.2f GB RAM\n\n', ...
        d1, d2, numFrames,length(unique(frameIndices(:,1))),prod(Ysiz)*8/(2^30));
else
     fprintf('\nThe data has been mapped to a hard disk. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', ...
        d1, d2, numFrames,prod(Ysiz)*8/(2^30));
end

%% Run motion correction (NoRMCorre), if desired

run_MC = input('Would you like to motion-correct the data? (y/n, default y)\n','s');
if isempty(run_MC)
    run_MC = 'y';
end

fprintf('Where would you like to save the filtered & motion-corrected data...\n')
% MC_fnam = fullfile(uigetdir(),[fnam(1:end-4),'_MCfilt.mat']);
MC_fnam = fullfile(fdir,[fnam(1:end-4),'_MCfilt.mat']);

if strcmp(run_MC,'y')
    %% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
    options_r = NoRMCorreSetParms('d1',d1,'d2',d2,...
        'bin_width',50,'max_shift',20,'iter',1,...
        'use_parallel',true,'memmap',true,...
        'mem_filename',MC_fnam,...
        'output_type','memmap');
    
    %% register data and apply shifts to removed percentile
    tic; [M1,shifts1,template1] = normcorre_batch_onFilt(data,options_r); toc % register filtered data
    
    %% apply shifts to original (unfiltered) dataset and overwrite the data with registered version
    
    % read data into workspace in chunks, apply shifts in workspace, then write back to
    % matfile
    
    options_r.memmap = false;
    options_r.output_type = 'mat';
    
    chunkSize = 5000;
    
    if isempty(chunkSize)
        
        %default to 4GB chunkSize
        details = whos(data,'Y');
        data_type = details.class;
        
        % 16-, 32- and 64-bit data types supported
        if strcmp(data_type,'uint16')
            chunkSize = min(numFrames, round(2*(2^30/d1/d2)));
        elseif strcmp(data_type,'single')
            chunkSize = min(numFrames, round(4*(2^30/d1/d2)));
        elseif strcmp(data_type,'double')
            chunkSize = min(numFrames, round(8*(2^30/d1/d2)));
        end
        
    end
    
    %% after motion correction, create video of results (filtered uncorrected vs. filtered corrected) and let user choose either:
    % a) accept motion correction & overwrite original filtered
    % (uncorrected) data with filtered (corrected), but save raw corrected
    % to another Matfile (e.g. M1)
    % b) reject motion correction results & start over
    
    tic
    for i = 1:chunkSize:numFrames
        
        sframemsg = ['Registering frame ',num2str(i),' to frame ',num2str(i+min(chunkSize,numFrames-i+1)-1),' of ',num2str(numFrames), ' total frames'];
        disp(sframemsg)
        
        if length(size(data,'Y')) == 3
            Ytemp = data.Y(:,:,i:(i+min(chunkSize,numFrames-i+1)-1));
            Ytemp = apply_shifts(Ytemp,shifts1(i:(i+min(chunkSize,numFrames-i+1)-1)),options_r);
            data.Y(:,:,i:(i+min(chunkSize,numFrames-i+1)-1)) = Ytemp;
        elseif length(size(data,'Y')) == 2
            Ytemp = data.Y(:,i:(i+min(chunkSize,numFrames-i+1)-1));
            Ytemp = reshape(Ytemp,d1,d2,[]);
            Ytemp = apply_shifts(Ytemp,shifts1(i:(i+min(chunkSize,numFrames-i+1)-1)),options_r);
            data.Y(:,i:(i+min(chunkSize,numFrames-i+1)-1)) = reshape(Ytemp,d1*d2,[]);     
        end    
        
    end
    toc
    clear Ytemp;

end

%% Pre-process data in patches to reduce computational load

% construct patches

if length(size(M1,'Yfilt')) == 3
    patch_size = [100 100];

    x_start = 1:patch_size(1):Ysiz(1);
    x_end = min(x_start+patch_size(1)-1,Ysiz(1));
    y_start = 1:patch_size(2):Ysiz(2);
    y_end = min(y_start + patch_size(2)-1,Ysiz(2));
    
    [X1,Y1]=meshgrid(x_start,y_start); [X2,Y2] = meshgrid(x_end,y_end);
    patches = mat2cell([X1(:),X2(:),Y1(:),Y2(:)],ones(numel(X1),1),2*length(Ysiz(1:end-1)));
    
    clear X1 X2 x_end x_start Y1 Y2 y_end y_start;
    
    Ythresh = zeros(d1,d2,2,'uint16');
    % M1.Ythresh = Ythresh;
    data.Ythresh = Ythresh;
    
elseif length(size(M1,'Yfilt')) == 2
    
    patch_size = 10000; % number of pixels in patch
    patches = 1:patch_size:(d1*d2);
    Ythresh = zeros(d1*d2,2,'uint16'); % Conor edit 08202017 -- make Ythresh a reshaped matrix
    % M1.Ythresh = Ythresh;
    data.Ythresh = Ythresh;
    
end

% do patch-wise processing
 
options = CNMFSetParms;
sig = 5;

tic

for i = 1:length(patches)
    
    % get patch's data
    
    if length(size(M1,'Yfilt')) == 3 
        Y_patch = M1.Yfilt(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
        %     Y_patch = data.Yfilt(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
        [d1p,d2p,~] = size(Y_patch);
        Y_patch = reshape(Y_patch,d1p*d2p,Ysiz(3));
    elseif length(size(M1,'Yfilt')) == 2
        Y_patch = M1.Yfilt(temp,:);
    end

    temp = patches(i):(patches(i)+min(patch_size,(d1*d2)-patches(i)+1)-1);

    % subtract each pixel's temporal median
    Y_patch = bsxfun(@minus,Y_patch,median(Y_patch,2));
        
    % get spectral estimate of noise for each pixel
    Ynoise_patch = get_noise_fft(Y_patch,options);
    
    %threshold patch with cut-off noise*sig -- anything below threshold is
    %set to 0
    Y_patch(bsxfun(@lt, Y_patch, Ynoise_patch*sig)) = 0;
    
    if length(size(data,'Ythresh')) == 3 
%         M1.Ythresh(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),1:numFrames) = Y_patch;
        data.Ythresh(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),1:numFrames) = Y_patch;
    elseif length(size(data,'Ythresh')) == 2 
        data.Ythresh(temp,1:numFrames) = Y_patch;
    end
   
end
clear temp Y_patch Ynoise_patch Ythresh 

toc

%% correlation image 

chunkSize = 799;
x_start = 1:chunkSize:numFrames;
x_end   = min(x_start + chunkSize - 1,numFrames);
chunks = mat2cell([x_start',x_end'],ones(length(x_start),1),2);
clear x_start x_end;

tic

Cn = zeros(d1,d2,length(chunks));
sz = [1 2];

for i = 1:length(chunks)
    
%     chunkData = M1.Ythresh(:,:,chunks{i}(1):chunks{i}(2));
%     chunkData = data.Ythresh(:,:,chunks{i}(1):chunks{i}(2));
    chunkData = data.Ythresh(:,chunks{i}(1):chunks{i}(2));
    chunkData = reshape(chunkData,d1,d2,[]);
    
    Cn(:,:,i) = correlation_image(cast(chunkData,'double'),sz,d1,d2);
    
end
    
toc

clear chunkData

save(fullfile(fdir,[fnam(1:end-4),'_CorrImage.mat']),'Cn');

%% morphology on the multiple correlation images

se = strel('disk',2); %structuring element for opening

cMin = 0.8; %minimum correlation for thresholding

thresholded = Cn; thresholded(bsxfun(@lt,Cn,cMin))=0;
opened = zeros(size(thresholded));
for slice = 1:size(thresholded,3)
    opened(:,:,slice) = imopen(thresholded(:,:,slice),se); %morpho-opening with strel to remove outlier pixels
end

centers = cell(1,size(thresholded,3));
for slice = 1:size(thresholded,3)
    temp = imregionalmax(opened(:,:,slice)); % find local maxima of opened image
    ROIs = bwconncomp(temp); %find ROIs with 8-pixel connectivity on the regional maxima image
    centers{slice} = zeros(ROIs.NumObjects,2);
    for roi = 1:ROIs.NumObjects
        [y,x] = ind2sub([d1 d2],ROIs.PixelIdxList{roi}(3));
        centers{slice}(roi,:) = [x, y];
    end
end


%verify ROI detection by overlaying on correlation image
for slice = 1:size(thresholded,3)
    imagesc(thresholded(:,:,slice)); hold on; scatter(centers{slice}(:,1),centers{slice}(:,2),'ro');
    pause; hold off;
end


%% merge redundant cells found multiply in different correlation images


dmin = 2;

allROIs = cell2mat(centers');

distMatrix = sqrt(bsxfun(@minus, allROIs(:,1), allROIs(:,1)').^2 + bsxfun(@minus, allROIs(:,2), allROIs(:,2)').^2);
flag_merge = (distMatrix<=dmin);
[l,c] = graph_connected_comp(sparse(flag_merge)); 
MC = bsxfun(@eq, reshape(l, [],1), 1:c);

nr = size(MC,1);
total = size(MC,2);
nonMergers = find(sum(MC,1)==1);
mergers = find(sum(MC,1)>1);

newCenters = zeros(nr,2);
ind_del = [];

for m = 1:total
    
    IDs = find(MC(:,m));
    
    if length(IDs)==1
        newCenters(IDs,:) = allROIs(IDs,:);
    else
        newCenters(IDs(1),:) = mean(allROIs(IDs,:));
        ind_del = [ind_del;IDs(2:end)];
    end
    
end

newCenters(ind_del,:) = [];

%view neurons overlaid on correlation image
imagesc(mean(thresholded,3)); hold on; scatter(newCenters(:,1),newCenters(:,2),'ro');

%% manually exclude neurons by drawing a rectangle around regions of the FOV
trim_flag = input('Do you want to select any neurons for deletion? (y/n, default n)\n','s');


imagesc(mean(thresholded,3)); hold on; scatter(newCenters(:,1),newCenters(:,2),'ro');
while strcmp(trim_flag,'y')
    excludeBox = getrect(gcf);
    ind_del = find(newCenters(:,1)>excludeBox(1) & newCenters(:,1) < excludeBox(1)+excludeBox(3) ...
        & newCenters(:,2) > excludeBox(2) & newCenters(:,2) < excludeBox(2)+excludeBox(4));
    newCenters(ind_del,:) = [];
    trim_flag = input('Do you want to select any neurons for deletion? (y/n, default n)\n','s');
    if strcmp(trim_flag,'n')
        break
    end
end

close gcf;

nr = length(newCenters);

imagesc(mean(thresholded,3)); hold on; scatter(newCenters(:,1),newCenters(:,2),'ro');


%% do a series of low-rank nNMFs on chunks of data centered on the cell locations identified in the above steps to get temporal/spatial components

% A = zeros(nr,d1*d2);
% C = zeros(numFrames,nr);

A = [];
C = [];

rectBounds = [6 6];

% already_done = [];

tic
for i = 1:nr
    temp = round(newCenters(i,:));
    x = temp(1); y = temp(2);
    neighbors = find(newCenters(:,1) > x - rectBounds(1) & newCenters(:,1) < x + rectBounds(1) & newCenters(:,2) > y - rectBounds(2) & newCenters(:,2) < y + rectBounds(2));
%     temp = ismember(neighbors,already_done);
    currCellData = double(data.Ythresh((x-rectBounds(1)):(x+rectBounds(1)),(y-rectBounds(2)):(y+rectBounds(2)),:));
    tempSz = size(currCellData);
    rank = length(neighbors);
    [w,h] = nnmf(reshape(currCellData,prod(tempSz(1:2)),[])',rank);
%     C(:,neighbors(~temp)) = w;
    C = [C,w];
    mask = zeros(d1,d2,rank); mask((x-rectBounds(1)):(x+rectBounds(1)),(y-rectBounds(2)):(y+rectBounds(2)),:)=reshape(h',tempSz(1),tempSz(2),rank);
%     A(neighbors(~temp),:) = reshape(mask,d1*d2,[]);
    A = [A,reshape(mask,d1*d2,[])];
%     already_done = [already_done;neighbors(~temp)];
end
toc

Atrim = A;

%spatially trim and remove disconnected components 

se = strel('disk', 2);
thr = 0.02;

for i=1:size(Atrim,2)
    ai = A(:,i);
    ai_open = imopen(reshape(ai,d1,d2), se);
    
    temp = ai_open>max(ai)*thr;
    l = bwlabel(temp, 8);   % remove disconnected components
    [~, ind_max] = max(ai_open(:));
    
    ai(l(:)~=l(ind_max)) = 0;
    ai(ai<max(ai)*0.3)=0;
    
    ai(ai>0) = A(ai>0,i);
    
    Atrim(:,i) = ai(:);
end

C = C';


%visualize spatial components before and after trimming to check parameters
figure;
for i = 1:nr
    subplot(121);
    imagesc(reshape(A(i,:),d1,d2));
    subplot(122);
    imagesc(reshape(Atrim(i,:),d1,d2));
    pause;
end

maxNorm = bsxfun(@times,Atrim,(1./max(Atrim,[],2)));

Csmooth = zeros(size(C));
for i = 1:nr
    Csmooth(:,i) = locsmooth(C(:,i),10,3);
end
    
reconstruction = Csmooth*maxNorm;

%create video of reconstruction
for i = 1:numFrames
    imagesc(reshape(reconstruction(i,:),d1,d2)); 
    pause(0.05);
end


%% create video 
% v = VideoWriter('example_spontaneous_activity.avi');
% open(v)
% 
% imagesc(reshape(reconstruction(1,:),d1,d2));
% set(gca,'nextplot','replacechildren');
% 
% for i = 2:numFrames
%     imagesc(reshape(reconstruction(i,:),d1,d2));
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% 
% close(v)
    
    
    
    
    


%% extract traces from original video

tic

nr = size(newCenters,1);

traces = zeros(nr,numFrames);
backgrounds = zeros(nr,numFrames);

roi_r = 5;
annulus_r = 8;
tol = 2;

for i = 1:length(chunks)
    
    sframemsg = ['Extracting traces from frame ',num2str(chunks{i}(1)),' to frame ',num2str(chunks{i}(2)),' of ',num2str(numFrames), ' total frames'];
    disp(sframemsg)
    
    chunkData = data.Y(:,:,chunks{i}(1):chunks{i}(2));
    
%     chunkData = M1.Yfilt(:,:,chunks{i}(1):chunks{i}(2));
    chunkData = cast(reshape(chunkData,d1*d2,[]),'double');
    
    for j= 1:nr
        
        center = round(newCenters(j,:));
        [x,y]=meshgrid(-(center(1)-1):(d1-center(1)),-(center(2)-1):(d2-center(2)));
        
        roi_mask=((x.^2+y.^2)<=roi_r^2);
        
        annulus_mask = (x.^2 + y.^2)<annulus_r^2 & (x.^2 + y.^2)>(annulus_r-tol)^2;
        
        roi_mask = reshape(roi_mask,d1*d2,1);
        annulus_mask = reshape(annulus_mask,d1*d2,1);
        
        traces(j,chunks{i}(1):chunks{i}(2)) = mean(chunkData(roi_mask,:));
        backgrounds(j,chunks{i}(1):chunks{i}(2)) = median(chunkData(annulus_mask,:));
        
    end
   
end

toc

save(fullfile(fdir,[fnam(1:end-4),'_traces.mat']),'traces','backgrounds','newCenters');


%%

smooth_bg = zeros(size(backgrounds));

for i = 1:nr
    smooth_bg(i,:) = locsmooth(backgrounds(i,:),10,10,5);
end
    

figure;
for i = 1:nr
    subplot(121); plot(traces(i,:)); hold on; pause; plot(backgrounds(i,:)); hold off;
    subplot(122); plot((traces(i,:)-backgrounds(i,:))./locsmooth(backgrounds(i,:),10,10,5)); pause;  
end


DFF = (traces-backgrounds)./smooth_bg;
DFF_smooth = DFF;
for i = 1:nr
    DFF_smooth(i,:) = locsmooth(DFF(i,:),10,3);
end

DFF_smooth = DFF_smooth';
deriv = diff(DFF_smooth);
noise_est = mean(deriv)+2.5*std(deriv);
events = bsxfun(@gt,deriv,noise_est);
nonZero = events(sum(events,2)>1,:);

    


chunkSize = 5000;

tic
for i = 1:chunkSize:numFrames
        
        sframemsg = ['Transferring registered Yfilt frame ',num2str(i),' to frame ',num2str(i+min(chunkSize,numFrames-i+1)-1),' of ',num2str(numFrames), ' total frames'];
        disp(sframemsg)
        
        Ytemp = M1.Yfilt(:,:,i:(i+min(chunkSize,numFrames-i+1)-1));
        
        data.Yfilt(:,:,i:(i+min(chunkSize,numFrames-i+1)-1)) = Ytemp;
        
end
toc    
    
    




    
