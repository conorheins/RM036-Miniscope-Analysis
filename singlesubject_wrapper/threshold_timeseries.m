function [ thresholded,thresholds ] = threshold_timeseries( timeseries,type,factor )
%threshold_timeseries: Applies row-wise thresholds to a multidimensional timeseries (assuming rows are individual timeseries) 
%   to denoise timeseries, using either standard deviations from the mean
%   or median absolute deviations from the median

if size(timeseries,1) > size(timeseries,2)
    timeseries = timeseries';
end

if ~exist('type','var') || isempty(type)
    type = 'mean';
end

if strcmp(type,'mean')
    thresholds = mean(timeseries,2) + factor*std(timeseries,0,2);
elseif strcmp(type,'median')
    thresholds = median(timeseries,2) + factor*mad(timeseries,1,2);
end

thresholded = timeseries; thresholded(bsxfun(@gt,timeseries,thresholds)) = 1;

end

