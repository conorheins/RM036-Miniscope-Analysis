function viz_array = plot_ndarray(X,transpose_flag,max_norm,color)
% PLOT_NDARRAY Plots multidimensional array of data (usually time-series)
% by stacking them vertically on top of one another, optional normalization
% and color options
%   INPUTS: 
%           X -- data matrix, of either (T x N) or (N x T), where T is
%           samples/observations, N is dimensionality/number of lines to be plotted
%
%           tranpose_flag -- flag to tranpose data, to achieve dimensions [T x N]. If [T x N], then
%           tranpose_flag can be left empty or set to 0. If 1, then time is
%           assumed to be in the second dimension, N in the first, and the
%           data matrix is accordingly transposed.
%
%           max_norm -- flag (True or False, 1 or 0) to normalize data to
%           max value
%           
%           colors -- N x 3 matrix of colors or single string to plot each
%           dimension/unit n with

% check if data needs to be transposed
if ~exist('transpose_flag','var') || ~transpose_flag
    data = X;
elseif transpose_flag
    data = X';
end

% max normalize
if exist('max_norm','var')
    if max_norm
        data = bsxfun(@rdivide,data,max(data,[],1));
        data(isnan(data)) = 0;
    end
end

% use top 90% of max variance of the population to set interval between
% traces

interval = 0.99 * max(range(data));

% uses interval and number of dimensions to create shift_matrix to translate all traces by for plotting 
[num_samples,num_dims] = size(data);
shift_matrix = repmat([1:num_dims]*interval,num_samples,1);
% viz_array = data + shift_matrix;
viz_array = data - shift_matrix; % edit on May 5 2018

figure;

if exist('color','var')
    if ischar(color)
        plot(viz_array,color); axis tight;
        return
    end
    if size(color,2) == 3
        all_lines = plot(viz_array); axis tight;
        if length(color) < num_dims
            error('Must have at least as many colors as there are traces!')
        elseif length(color) > num_dims
            color = color(randperm(length(color),num_dims),:);
        end
        cc = mat2cell(color,ones(num_dims,1),repmat(3,1,1));
        set(all_lines,{'color'},cc);
    elseif size(color,2) ~= 3
        error('Color matrix must be (at-least-N) x 3!')
    end
    
elseif ~exist('color','var') || isempty('color')
    all_lines = plot(viz_array); axis tight;
    color = jet(3*num_dims);
    color = color(randperm(length(color),num_dims),:);
    cc = mat2cell(color,ones(num_dims,1),repmat(3,1,1));
    set(all_lines,{'color'},cc);
end

end

