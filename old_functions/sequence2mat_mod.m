%% Reads in sequence of .tif files and memory maps them to a .mat file. Assumes images are in NeuralMapper output format:
% (i.e. unsigned 16-bit .tifs with naming format: 'image_XXX_XXX.tif', where XXX is  an integer
% R. Conor Heins, 05022017
% edit 1.1, Conor and Tarun Madangopal, 05072017
% edit 1.2, Conor Heins, 06092017 -- allows sequences of .tifs to be
%   filtered during read-in, also explicitly accomodates continuous (non-trial-based) recordings 
% edit 1.3, Conor Heins, 08192017 -- reshape Y data to 2D for easier handling
% edit 1.4, Conor heins, 10172017 -- undo edit 1.3, uses 3D tensors for all
% operations. Makes things easier for cell detection later on

function [data,Ysiz,sortedIndices] = sequence2mat_mod(data,filter_flag,psf,add_data_flag,trial_based,trial_length,exInd)

% ask user what the first trial index is
firstInd = input('Enter number of the first trial in this stack (default = 1):\n');

%choose image directory (folder with image sequence) and get 'img' struct
%of file names and info
fprintf('Choose a folder to load new image data from...\n')
stack_dir = uigetdir();
imgs = dir([stack_dir, filesep,'*.tif']);

%choose which frames to throw away (e.g. first frame of every trial, which is often very low fluorescence)
if ~exist('exInd','var') || nargin < 7
    frames2exclude = input('Are there any frames you would like to exclude?\n(enter single digit or vector, or press enter if you want to keep everything):\n');
else
    frames2exclude = exInd;
end

% determine first trial number of stack, in case of trial-based analysis,
% and create corresponding trial-and-frame array for sorting/excluding frames later on
if strcmp(trial_based,'y')
    if isempty(firstInd)
        firstInd = 1;
    end
    trialandframe = zeros(length(imgs),2);
else %otherwise, create frame array for sorting/excluding frames later on
    frames = zeros(length(imgs),1); 
end

% fill out (trial and) frame list for sorting/excluding frames
for file = 1:length(imgs)
    fnam = imgs(file).name;
    IDs = regexp(fnam,'\d*\d*','match'); % use a regular expression to find trial index and frame index within each file name
    if strcmp(trial_based,'y')
        trialandframe(file,1)=str2num(IDs{1})+firstInd-1;
        trialandframe(file,2)=str2num(IDs{2});
    else
        frames(file) = str2num(IDs{2});
    end
end

% throw out trials with less than trial_length frames 
if strcmp(trial_based,'y')
    firstCol = trialandframe(:,1);
    allTrials = unique(firstCol);
    badTrials = zeros(size(firstCol));
    for i = 1:length(allTrials)
        if length(trialandframe(firstCol==allTrials(i),:)) < trial_length
            badTrials = badTrials + firstCol==allTrials(i);
        end
    end
    trialandframe(logical(badTrials),:) = [];
    imgs(logical(badTrials))=[];
end

% find indices of frames to throw out (i.e. the first 5 frames of every
% trial, etc.)
excludeIDs = [];
if ~isempty(frames2exclude)
    for i = 1:length(frames2exclude)
        if strcmp(trial_based,'y')
            excludeIDs = [excludeIDs, find(trialandframe(:,2)==frames2exclude(i))];
        else
            excludeIDs = [excludeIDs, find(frames == frames2exclude(i))];
        end
    end
end

% throw out bad frames and sort frames and images
if strcmp(trial_based,'y')
    trialandframe(excludeIDs,:) = [];
    [sortedIndices, sortInds] = sortrows(trialandframe,[1 2]);
else
    frames(excludeIDs) = [];
    [sortedIndices,sortInds] = sort(frames,'ascend');
end

imgs(excludeIDs)=[]; 
sortedImgs = imgs(sortInds);

tic;

if strcmp(filter_flag,'y')
    fprintf('converting the selected image sequence to *.mat version and filtering...\n');
else
    fprintf('converting the selected image sequence to *.mat version...\n');
end

info = imfinfo(fullfile(stack_dir,sortedImgs(1).name)); %read first image from sequence to get image info
d1 = info.Height;   % height of the image 
d2 = info.Width;    % width of the image 
T = length(sortedImgs);
Ysiz = [d1, d2, T]';


if add_data_flag
%     firstFrame = size(data.Y,2)+1; % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
    firstFrame = size(data.Y,3)+1; % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on 
else
%     data.Y = uint16(zeros(d1*d2,2)); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
%     data.Yfilt = uint16(zeros(d1*d2,2)); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
    data.Y = uint16(zeros(d1,d2,2)); % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on
    data.Yfilt = uint16(zeros(d1,d2,2)); % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on
    firstFrame = 1;
end

Tchunk = min(T, round(2^29/d1/d2)); %each chunk uses at most 4GB

firstChunk = sequence_bigread2(stack_dir,sortedImgs,1,Tchunk);
% data.Y(:,firstFrame:(firstFrame+Tchunk-1)) = reshape(firstChunk,d1*d2,[]); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
data.Y(:,:,firstFrame:(firstFrame+Tchunk-1)) = firstChunk; % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on

if strcmp(filter_flag,'y')
    temp = imfilter(firstChunk,psf,'replicate');
%     data.Yfilt(:,firstFrame:(firstFrame+Tchunk-1)) = reshape(temp,d1*d2,[]); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
    data.Y(:,:,firstFrame:(firstFrame+Tchunk-1)) = temp; % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on
end

if Tchunk==T
    return; 
else
    t0 = Tchunk+1; 
    while t0<=T
        num2read = min(t0+Tchunk-1, T) - t0 + 1; 
        tmpY = sequence_bigread2(stack_dir,sortedImgs, t0, num2read); 
%         data.Y(:,(firstFrame-1)+(1:num2read)+t0-1) = reshape(tmpY,d1*d2,[]); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
        data.Y(:,:,(firstFrame-1)+(1:num2read)+t0-1) = tmpY; % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on
        if strcmp(filter_flag,'y')
            tmpY_f = imfilter(tmpY,psf,'replicate'); 
%             data.Yfilt(:,(firstFrame-1)+(1:num2read)+t0-1) = reshape(tmpY_f,d1*d2,[]); % edit 1.3 08192017 -- initialize Y to 2D shape, makes calling it easier
            data.Yfilt(:,:,(firstFrame-1)+(1:num2read)+t0-1) = tmpY_f; % edit 1.4 10172017 -- back to old version, initialie Y to 3D shape, makes cell detection easier later on
        end
        t0 = t0 + num2read; 
    end 
end

if strcmp(filter_flag,'y')
    fprintf('Time cost in converting data to *.mat file and filtering:     %.2f seconds\n', toc); % edit 1.4 10172017 -- changed display message
else
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
end

end