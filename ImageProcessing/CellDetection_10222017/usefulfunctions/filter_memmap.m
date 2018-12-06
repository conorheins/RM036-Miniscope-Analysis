function [ memm_file,chunkSize ] = filter_memmap(memm_file,psf,chunkSize)
%FILTER_MEMMAP Looks for filtered version of data in object, and if absent,
% creates one
%   Initializes Yfilt data if there is none in the memmap object, optimizes
%   chunk size if none is provided


details = whos(memm_file,'Y');
size = details.size;
data_type = details.class;

d1 = size(1); d2 = size(2); T = size(3);

if ~exist('chunkSize','var') || isempty(chunkSize) || nargin < 3  
    % default to 4GB chunk size
    if strcmp(data_type,'uint16')
        chunkSize = min(T, round(2*(2^30/d1/d2))); 
    elseif strcmp(data_type,'single')
        chunkSize = min(T, round(4*(2^30/d1/d2))); 
    elseif strcmp(data_type,'double')
        chunkSize = min(T, round(8*(2^30/d1/d2))); 
    end
end
        
if isempty(whos(memm_file,'Yfilt')) || ismatrix(memm_file.Yfilt)
    Yfilt = zeros(d1,d2,chunkSize,data_type);
    memm_file.Yfilt = Yfilt;
end

for i = 1:chunkSize:T
    
    sframemsg = ['Filtering frame ',num2str(i),' to frame ',num2str(min(T,i+chunkSize-1)),' of ',num2str(T), ' total frames'];
    disp(sframemsg)
    
    Ytemp = memm_file.Y(:,:,i:(i+min(chunkSize,T-i+1)-1));
    
    memm_file.Yfilt(:,:,i:(i+min(chunkSize,T-i+1)-1)) = imfilter(Ytemp,psf,'replicate');
    
end

end

