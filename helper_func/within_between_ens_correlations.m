function [within_corrz, btwn_corrz] = within_between_ens_correlations(data2use,assembly_info,display_flag)
%within_between_ens_correlations 
% Generates correlation matrices for each ensemble, as well as a null model
% of 'between-ensemble' correlations by randomly selecting subsets of
% neurons and looking at their correlation distributions
%   INPUTS: data2use -- any N x T matrix of timeseries
%           assembly_info -- a structure array with two fields:
%           'membership' and 'r_sq'

if ~exist('display_flag','var') || isempty(display_flag)
    display_flag = false;
end

assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership)); 
[ensembles_sorted,srt] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
unique_ensembles = unique(ensembles_sorted);

membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
for neuron = 1:length(assembly_info)
    membership_matrix(neuron,assembly_info(neuron).membership) = 1;
end

within_corrz = cell(1,length(unique_ensembles));
btwn_corrz = cell(1,length(unique_ensembles));

max_norm = bsxfun(@rdivide,bsxfun(@minus,data2use,min(data2use,[],2)),(max(data2use,[],2) -  min(data2use,[],2)));

for i = 1:length(unique_ensembles)
    ens = unique_ensembles(i);
    ens_members = find(membership_matrix(:,i) == 1);
    within_corrz{i} = corrcoef(max_norm(ens_members,:)');
    btwn_corrz{i} = corrcoef(max_norm(randperm(size(max_norm,1),length(ens_members)),:)'); % random, matched subset
end

if display_flag
    for i = 1:length(unique_ensembles)
        temp1 = within_corrz{i};
        temp2 = btwn_corrz{i};
        min_corr = min([temp1(:);temp2(:)]);
        subplot(121);
        imagesc(within_corrz{i}); caxis([min_corr 1]);
        subplot(122)
        imagesc(btwn_corrz{i}); caxis([min_corr 1]);
        pause;
    end
end
    
    
    