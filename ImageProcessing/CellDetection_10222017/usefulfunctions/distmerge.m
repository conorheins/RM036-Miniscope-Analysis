function [ newCenters,numROIs ] = distmerge( centers,dmin,display_flag,display_img )
% Determines every inter-ROIs euclidean distance, thresholds based on
% dmin, then does conncomps to find groups of centroids that probably are
% the same neuron

if ~exist('dmin','var')
    dmin = 2;
end

distMatrix = sqrt(bsxfun(@minus, centers(:,1), centers(:,1)').^2 + bsxfun(@minus, centers(:,2), centers(:,2)').^2);
flag_merge = (distMatrix<=dmin);
[l,c] = graph_connected_comp(sparse(flag_merge)); 
MC = bsxfun(@eq, reshape(l, [],1), 1:c);

nr = size(MC,1);
total = size(MC,2);

newCenters = zeros(nr,2);
ind_del = [];

for m = 1:total
    
    IDs = find(MC(:,m));
    
    if length(IDs)==1
        newCenters(IDs,:) = centers(IDs,:);
    else
        newCenters(IDs(1),:) = mean(centers(IDs,:));
        ind_del = [ind_del;IDs(2:end)];
    end
    
end

newCenters(ind_del,:) = [];

if display_flag
    %view neurons overlaid on correlation image
    imagesc(display_img); hold on; scatter(newCenters(:,1),newCenters(:,2),'ro');
end

numROIs = length(newCenters);

end

