function [ A_merge, C_merge, K_merge ] = conor_neuron_merge( A, C, sz, merge_type, thresholds, display_flag )
%CONOR_NEURON_MERGE Merges neurons based on (either OR or AND) spatiotemporal correlations 
%   Inputs: A -- spatial components (sparse or full matrix)
%           C -- temporal components
%           sz -- spatial size in format [d1, d2]
%           merge_type -- whether to use A correlations ('spatial'), C
%           correlations ('temporal') or both ('both')
%           thresholds -- thresholds for merging, in format [A_threshold,
%           C_threshold]). If merge_type is not 'both', will use the first
%           or second element of thresholds as the respectively A- or C-correlation threshold
%           display_flag -- allows user to manually inspect each merge
%           candidate
%   Outputs: A_Merge -- new spatial components
%            C_merge -- new temporal components
%            K_merge -- new total neuron number, post merging 

if ~exist('merge_type','var')
    merge_type = 'both';
end
if ~exist('display_flag','var')
    display_flag = false;
end

A_thr = thresholds(1); C_thr = thresholds(2);

if ~ismember(merge_type,{'spatial','temporal','both'})
    error('Invalid input for merge_type -- must be "spatial", "temporal", or "both"')
    return;
end

if strcmp(merge_type,'both')
    
    temp = bsxfun(@times, A>0, 1./sqrt(sum(A>0))); % binarize spatial components and normalize by Frobenius energy
    A_overlap = temp'*temp;
    
    [K, ~] = size(C);   % number of neurons
    C_overlap = corr(C')-eye(K);

    flag_merge = (A_overlap > A_thr) & (C_overlap > C_thr);
    
elseif strcmp(merge_type,'spatial')
    
    temp = bsxfun(@times, A>0, 1./sqrt(sum(A>0))); % binarize spatial components and normalize by Frobenius energy
    A_overlap = temp'*temp;
    
    flag_merge = (A_overlap > A_thr);
    flag_merge = flag_merge-eye(size(flag_merge));
    
elseif strcmp(merge_type,'temporal')
    
    [K,~] = size(C);
    C_overlap = corr(C') - eye(K);
    
    flag_merge = (C_overlap > C_thr);
    flag_merge = flag_merge;
    
end

%% use Bron-Kerbosch algorithm to find cliques in adjacency matrix 'flag_merge'
% 
% if ~issparse(flag_merge)
%     flag_merge = sparse(flag_merge);
% end


%% maybe use connected components instead
[l,c] = graph_connected_comp(flag_merge); 
MC = bsxfun(@eq, reshape(l, [],1), 1:c);

nonMerge = find(sum(MC,1) == 1);
toMerge = MC(:,sum(MC,1)>1);
n2merge = size(toMerge,2);

if display_flag
    d1 = sz(1); d2 = sz(2);
    % optional visualization of found merge-candidates
    
    for i = 1:n2merge
        
        memberIDs = find(toMerge(:,i));
        spatial_members = A(:,memberIDs);
        spatial_members = bsxfun(@rdivide,spatial_members,max(spatial_members,[],1));
        temporal_members = C(memberIDs,:);
        
        clique_size = length(memberIDs);
        temp = sum(spatial_members,2);
        ind_nhood = find(temp);
        [rsub,csub] = ind2sub([d1,d2],ind_nhood);
        rbounds = max(min(rsub)-5,1):min(max(rsub)+5,d1);
        cbounds = max(min(csub)-5,1):min(max(csub)+5,d2);
        all_members = zeros(length(rbounds),length(cbounds));
        
        figure('pos',[200 1400 800 200]);
        for j = 1:clique_size
            Atemp = reshape(full(spatial_members(:,j)),d1,d2);
            all_members = all_members + Atemp(rbounds,cbounds);
            subplot(121)
            imagesc(all_members)
            title(['Clique No. ',num2str(i)]);
            subplot(122);
            plot(temporal_members(j,1:10000)); hold on;
            pause;
        end
  
        close(gcf)
    end
end


A_merge = sparse(zeros(size(A))); 
C_merge = zeros(size(C));

totalComponents = size(MC,2);

for i = 1:totalComponents
    
    IDs = find(MC(:, i));   % IDs of neurons within this cluster
    
    if sum(MC(:,i),1) == 1
        
        A_merge(:,IDs) = A(:,IDs);
        C_merge(IDs,:) = C(IDs,:);
        
    else
        
    
        % determine searching area
        active_pixel = find(sum(A(:,IDs), 2)>0);
        
        % update spatial/temporal components of the merged neuron
        recon = A(active_pixel, IDs)*C(IDs, :);
        ci = C(IDs(1), :);
        for miter=1:10
            ai = recon*ci'/(ci*ci');
            ci = ai'*recon/(ai'*ai);
        end
        
        A_merge(active_pixel,IDs(1)) = ai;
        C_merge(IDs(1),:) = ci;
        
        ind_del(IDs(2:end))=true;
    end
    
end

A_merge(:,ind_del) = [];
C_merge(ind_del,:) = [];

K_merge = size(A_merge,2); % new number of neurons after merging (for subsequent steps)

end

