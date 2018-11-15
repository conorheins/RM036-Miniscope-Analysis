function [ A_init,numROIs ] = initialize_A(allROIs,data_img,rectBounds,options)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('options','var')
    options.medfilt_param = [3,3];
    options.nrgthr = 0.5;
    options.close_elem = strel('square',3);
end
if ~exist('rectBounds','var')
    rectBounds = [12 12];
end

[d1,d2] = size(data_img);
numROIs = size(allROIs,1);

data_thr = zeros(d1,d2,numROIs);
numPixels = zeros(1,numROIs);

for i = 1:numROIs
   
    r = round(allROIs(i,2));
    c = round(allROIs(i,1));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    relative_r = r - rsub(1) + 1;
    relative_c = c - csub(1) + 1;
    
    center_ind = sub2ind([length(rsub),length(csub)],relative_r,relative_c);
    
    data_temp = data_img(rsub,csub);
    data_thr_temp = cleanup_footprints(data_temp,center_ind,options);
    data_thr(rsub,csub,i) = circular_constraints(data_thr_temp,center_ind);
    numPixels(i) = nnz(data_thr_temp);
    
end

data_thr(:,:,numPixels<=5) = [];
numROIs = size(data_thr,3);

A_init = reshape(data_thr,d1*d2,numROIs);

end

