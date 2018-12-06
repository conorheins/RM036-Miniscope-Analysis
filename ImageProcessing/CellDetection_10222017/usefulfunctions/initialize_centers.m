function [ allROIs ] = initialize_centers(Cn,STD_image,options,display_flag)
% INITIALIZE_CENTERS, seed ROI centers based on product of summary
% statistics (correlation image(s) and standard-deviation image(s))
%   Steps of morphological operations for finding ROI candidates:
%   1. threshold both images by top n-th percentile of intensity-values
%   (determined by values in <options.thresholds> )
%   2. Take product of both image-stacks as summary statistic
%   3. morphologically open summary statistic stacks using user-defined structuring element (given by
%   <options.se>)
%   4. seed ROI-candidates at local maxima of the opened images (one set of
%   ROIs for each slice in the stack).
%   5. Concatenate all ROIs found across all slices into 'allROIs' --
%   simple merging should be used after this to eliminate multiply-detected
%   ROIs in the same location

if ~isfield(options,'se'); se = strel('disk',2); else; se = options.se; end %structuring element for opening
if ~isfield(options,'thresholds'); thresholds = [90 90]; else; thresholds = options.thresholds; end %structuring element for opening
if ~exist('display_flag','var'); display_flag = true; end % flag for displaying detected centers overlaid on each statistic-slice

[d1,d2,numSlices] = size(Cn);

Cn_thr_image = reshape(Cn,d1*d2,[]);
STD_thr_image = reshape(STD_image,d1*d2,[]);
Cn_thr =  prctile(Cn_thr_image,thresholds(1),1);
STD_thr = prctile(STD_thr_image,thresholds(2),1);
Cn_thr_image(bsxfun(@lt,Cn_thr_image,Cn_thr))=0;
STD_thr_image(bsxfun(@lt,STD_thr_image,STD_thr))=0;

combined = reshape( (Cn_thr_image .* STD_thr_image), d1, d2, numSlices);

opened = zeros(size(combined));
for slice = 1:numSlices
    opened(:,:,slice) = imopen(combined(:,:,slice),se); %morpho-opening with strel to remove outlier pixels
end

centers = cell(1,numSlices);
for slice = 1:size(opened,3)
    temp = imregionalmax(opened(:,:,slice)); % find local maxima of opened image-slice
    ROIs = bwconncomp(temp); %find ROIs with 8-pixel connectivity on the regional maxima image
    centers{slice} = zeros(ROIs.NumObjects,2);
    for roi = 1:ROIs.NumObjects
        [y,x] = ind2sub([d1 d2],ROIs.PixelIdxList{roi}(3));
        centers{slice}(roi,:) = [x, y];
    end
end

if display_flag
    figure;
    %  verify ROI detection by overlaying on correlation image
    for slice = 1:numSlices
        imagesc(combined(:,:,slice)); hold on; scatter(centers{slice}(:,1),centers{slice}(:,2),'ro');
        pause; hold off;
    end
end

allROIs = cell2mat(centers');


end

