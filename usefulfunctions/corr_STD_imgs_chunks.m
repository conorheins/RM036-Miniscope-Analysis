function [ Cn,STD_image ] = corr_STD_imgs_chunks(data,chunkSize,Ysiz)
% generate sequence of correlation and standard deviation images from
% Ythresh (filtered & pre-processed) data, reading it in in chunks
%   If no chunk size provided, defaults to making one image for each trial

if ~exist('Ysiz','var') || isempty(Ysiz)
    Ysiz = size(data,'Y');
    d1 = Ysiz(1); d2 = Ysiz(2); numFrames = Ysiz(3);
else
    d1 = Ysiz(1); d2 = Ysiz(2); numFrames = Ysiz(3);
end

if ~exist('chunkSize','var') || isempty(chunkSize)
    
    frameIndices = data.sortedIndices;
    trials_ids = unique(frameIndices(:,1));
    chunkSize = length(find((frameIndices(:,1) == trial_ids(1)))); % default chunkSize to number of frames in a trial, using first trial as a reference for finding frame length
    
end

x_start = 1:chunkSize:numFrames;
x_end   = min(x_start + chunkSize - 1,numFrames);
chunks = mat2cell([x_start',x_end'],ones(length(x_start),1),2);

Cn = zeros(d1,d2,length(chunks));
STD_image = zeros(d1,d2,length(chunks));
sz = [1 2];

for i = 1:length(chunks)
    
    chunkData = data.Ythresh(:,chunks{i}(1):chunks{i}(2));
    chunkData = reshape(chunkData,d1,d2,length(chunks{i}(1):chunks{i}(2)));
    
    Cn(:,:,i) = correlation_image(cast(chunkData,'double'),sz,d1,d2);
    STD_image(:,:,i) = std(cast(chunkData,'double'),0,3);
    
end

end

