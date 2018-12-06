function [A_update] = updateSpatial_endoscope_noObj(Y,A,C,sn,num,method,search_method,sz)
%% udpate spatial components

%% inputs:
%   Y: d X T matrix, data
%   A: d x K matrix, spatial components
%   C: K x T matrix, temporal components
%   sn: 1 x K vector of noise values
%   num: scalar. If method=='hals', then num is the number of iterations to
%       update A; If method=='nnls', then num is the maximum number of neurons
%       overlapping at one pixel
%   method: method for updating the spatial components {'hals', 'nnls'}.
%       default: 'nnls'
%   search_method: method ('ellipse', 'dilate', or []) for determing pixel indices for 
%   updating each spatial component
%   size argument to initialize two fields of 'params' struct

%% Author: Pengcheng Zhou, Carnegie Mellon University.

%% input parameters number of iterations
if ~exist('method', 'var')||isempty(method)
    method = 'nnls';
end
if ~exist('num', 'var')||isempty(num)
    if strcmpi(method, 'nnls')
        num=5;
    else
        num = 10;
    end
end

%% determine the search locations

params.d1 = sz(1); params.d2 = sz(2);
IND = logical(determine_search_location(A, search_method,params));

%% estimate the noise
if and(strcmpi(method, 'hals_thresh') || strcmpi(method, 'nnls_thresh'), isempty(sn))
    %% estimate the noise for all pixels
    b0 =zeros(size(A,1), 1);
    sn = b0;
    parfor m=1:size(A,1)
        [b0(m), sn(m)] = estimate_baseline_noise(Y(m, :));
    end
    Y = bsxfun(@minus, Y, b0); 
end

%% update spatial components
if strcmpi(method, 'hals')
    A_update = HALS_spatial(Y, A, C, IND, num);
elseif strcmpi(method, 'hals_thresh')
    A_update = HALS_spatial_threshold(Y, A, C, IND, num, sn); 
elseif strcmpi(method, 'lars')
     [A_update, ~] = update_spatial_components_nb(Y,C,A); 
elseif strcmpi(method, 'nnls_thresh')&&(~isempty(IND_thresh)) 
    A_update = nnls_spatial_thresh(Y, A, C, IND, num, sn); 
else
    A_update = nnls_spatial(Y, A, C, IND, num);
end

end
