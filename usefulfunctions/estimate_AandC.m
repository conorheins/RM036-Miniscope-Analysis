function [A,C,cell_coords] = estimate_AandC(data,A_init,rectBounds,Ysiz)
% Estimate A and C with regression (for A) rank-1 NMF for C (plus
% background normalization)
%   Detailed explanation goes here

if ~exist('rectBounds','var')
    rectBounds = [12 12];
end

A = [];
C = [];

d1 = Ysiz(1); d2 = Ysiz(2); 
numROIs = size(A_init,2);

tic
Y = data.Y;
Ythresh = data.Ythresh;
frameIndices = data.sortedIndices;
fprintf('Time taken to load data into RAM: %.2f minutes\n',(toc/60))

cell_coords = com(A_init,d1,d2);

for i = 1:numROIs
    
    fprintf('Extracting neuron (%d / %d)\n',i,numROIs)
    r = round(cell_coords(i,1));
    c = round(cell_coords(i,2));
    rsub = max(1, -rectBounds(1)+r):min(d1, rectBounds(2)+r);
    csub = max(1, -rectBounds(1)+c):min(d2, rectBounds(2)+c);
    [cind, rind] = meshgrid(csub,rsub);
    [d1p,d2p] = size(cind);
    ind_nhood = sub2ind([d1,d2],rind(:),cind(:));
    
    relative_r = r - rsub(1) + 1;
    relative_c = c - csub(1) + 1;
    
    center_ind = sub2ind([length(rsub),length(csub)],relative_r,relative_c);
    
    HY_box = cast(Ythresh(ind_nhood,:),'double');
    Y_box = cast(reshape(Y(rsub,csub,:),d1p*d2p,[]),'double');
    ai = reshape(A_init(ind_nhood,i),d1p,d2p);
    
    y0 = HY_box(center_ind,:);
    tmp_corr = reshape(corr(y0', HY_box'), d1p, d2p);
    cell_pixels = find(tmp_corr > 0.3 & ai > 0);
    bg_pixels = find(tmp_corr < 0.3 & ai == 0);
    cell_data = Y_box(cell_pixels, :); % equivalent of 'F-raw'
    bg = median(Y_box(bg_pixels,:),1); % equivalent of 'F0'
    bg_smooth = zeros(size(bg)); % equivalent of 'F0_smoothed'
    trial_ids = unique(frameIndices(:,1));
    for trial = 1:length(trial_ids)
        thisTrial = find(frameIndices(:,1)==trial_ids(trial));
        bg_smooth(thisTrial) = locsmooth(bg(thisTrial),10,10,5);
    end
    
    ci = (mean(cell_data,1) - bg)./bg_smooth; % temporal component calculated as (Fraw-F0)/F0-style ratio
    
    
    %% estimate ai by fitting a linear regression 
    % with the raw data patch (pixel values over time) as regressands and
    % the background & temporal component as regressors
    T = length(ci);
    X = [ones(T,1), bg', ci'];
    temp = (X'*X)\(X'*Y_box');
    atemp = max(0, temp(3,:)');
    
    %% threshold the spatial shape and remove outliers
    % remove outliers
    temp =  full(atemp>quantile(atemp(:), 0.5));
    l = bwlabel(reshape(temp, d1p, d2p), 4);
    temp(l~=l(center_ind)) = false;
    atemp(~temp(:)) = 0;
    
    ai = zeros(d1*d2,1);
    ai(ind_nhood) = atemp;
    
    A = [A,ai];
    C = [C;ci];
    
end

cell_coords = com(A,d1,d2);

end

