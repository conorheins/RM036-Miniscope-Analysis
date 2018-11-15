function [ convolved_array ] = convolve_timeseries(data_array,transpose_flag,kernel)
%convolve_timeseries Quick way to do row- or column-wise convolution (e.g.
%of point-processes with transiently-increasing rates.
% assumes data is formatted as N x T, where N is dimensions/units, T is
% time

if ~exist('transpose_flag','var') || ~transpose_flag
    data = data_array;
elseif transpose_flag
    data = data_array';
end

if ~exist('kernel','var') || isempty(kernel)
    kernel = normpdf(-3:3,0,1.5);
end

convolved_array = zeros(size(data));

for n = 1:size(data,1)
    convolved_array(n,:) = conv(data(n,:),kernel,'same');
end

if ~exist('transpose_flag','var') || ~transpose_flag
    convolved_array = convolved_array;
elseif transpose_flag
    convolved_array = convolved_array';
end

end

