function [ data ] = apply_shifts_wrapper(data,shifts1,options,chunkSize,Ysiz)
% Apply shifts learned from NoRMCorre to original dataset in
% chunks
% read data into workspace in chunks, apply shifts in workspace, then write back to
% matfile

if ~exist('Ysiz','var') || isempty(Ysiz)
    Ysiz = size(data,'Y');
    d1 = Ysiz(1); d2 = Ysiz(2); numFrames = Ysiz(3);
else
    d1 = Ysiz(1); d2 = Ysiz(2); numFrames = Ysiz(3);
end

if ~exist('chunkSize','var') || isempty(chunkSize)
    
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
    

for i = 1:chunkSize:numFrames
    
    sframemsg = ['Registering frame ',num2str(i),' to frame ',num2str(i+min(chunkSize,numFrames-i+1)-1),' of ',num2str(numFrames), ' total frames'];
    disp(sframemsg)
    
    if length(size(data,'Y')) == 3
        Ytemp = data.Y(:,:,i:(i+min(chunkSize,numFrames-i+1)-1));
        Ytemp = apply_shifts(Ytemp,shifts1(i:(i+min(chunkSize,numFrames-i+1)-1)),options);
        data.Y(:,:,i:(i+min(chunkSize,numFrames-i+1)-1)) = Ytemp;
    elseif length(size(data,'Y')) == 2
        Ytemp = data.Y(:,i:(i+min(chunkSize,numFrames-i+1)-1));
        Ytemp = reshape(Ytemp,d1,d2,[]);
        Ytemp = apply_shifts(Ytemp,shifts1(i:(i+min(chunkSize,numFrames-i+1)-1)),options);
        data.Y(:,i:(i+min(chunkSize,numFrames-i+1)-1)) = reshape(Ytemp,d1*d2,[]);
    end
    
end


end

