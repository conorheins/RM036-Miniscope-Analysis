function all_ensembles_raster(data2use,assembly_info)
%all_ensembles_raster Plot time series of all ensembles, along with their
%associated correlation matrix (sorted as well)
%   Detailed explanation goes here

assembly_membership = {assembly_info.membership};
one_assembly_IDX = find(cellfun(@(x) isscalar(x),assembly_membership));

comms = [assembly_info(one_assembly_IDX).membership];

[X,Y,srt] = grid_communities(comms);

max_norm = bsxfun(@rdivide,bsxfun(@minus,data2use,min(data2use,[],2)),(max(data2use,[],2) -  min(data2use,[],2)));
data = max_norm(one_assembly_IDX(srt),:)';

% [ensembles_sorted,srt] = sort(cellfun(@(x) x,assembly_membership(one_assembly_IDX)));
% unique_ensembles = unique(ensembles_sorted);

% membership_matrix = zeros(length(assembly_info),length(unique_ensembles));
% for neuron = 1:length(assembly_info)
%     membership_matrix(neuron,assembly_info(neuron).membership) = 1;
% end

% colors = hsv(length(unique_ensembles));

% max_norm = bsxfun(@rdivide,bsxfun(@minus,data2use,min(data2use,[],2)),(max(data2use,[],2) -  min(data2use,[],2)));

% data = max_norm(one_assembly_IDX(srt),:)';

% interval
interval = 0.75;

% uses interval and number of dimensions to create shift_matrix to translate all traces by for plotting 
[num_samples,num_dims] = size(data);
shift_matrix = repmat([1:num_dims]*interval,num_samples,1);
viz_array = data - shift_matrix;

unique_ensembles = unique(comms);
colors = hsv(length(unique_ensembles));

% for i = 1:length(unique_ensembles)
%     ensemble_members = find(ensembles_sorted == unique_ensembles(i) ); % finds unique ensemble members, within indexing of already re-arranged data array
%     color_i = repmat(colors(i,:),length(ensemble_members),1);
%     cc = mat2cell(color_i,ones(length(ensemble_members),1),repmat(3,1,1));
%     set(all_lines(ensemble_members),{'color'},cc);
% end

sorted_communities = comms(srt);

%clean them up for visualization

corrMat = corrcoef(data);

for i = 1:length(unique_ensembles)
    ensemble_members = find(sorted_communities == unique_ensembles(i)); % finds unique ensemble members, within indexing of already ensemble-sorted data array
    if length(ensemble_members) == 1
        distances = squareform(pdist(corrMat));
        if ensemble_members == size(corrMat,1)
            [~,closest_idx] = min(distances(ensemble_members,[1:ensemble_members-1]));
        else
            [~,closest_idx] = min(distances(ensemble_members,[1:ensemble_members-1,(ensemble_members+1):end]));
        end
        sorted_communities(ensemble_members) = sorted_communities(closest_idx);
    end
end

[resorted_communities,srt] = sort(sorted_communities);
corrMat = corrMat(srt,srt);
unique_ensembles = unique(resorted_communities);
colors = hsv(length(unique_ensembles));
[X,Y] = grid_communities(resorted_communities);

data = max_norm(one_assembly_IDX(srt),:)';
viz_array = data - shift_matrix;

figure;        
subplot(121)
all_lines = plot(viz_array); axis tight;        
    
for i = 1:length(unique_ensembles)
    ensemble_members = find(resorted_communities == unique_ensembles(i) ); % finds unique ensemble members, within indexing of already ensemble-sorted data array
    color_i = repmat(colors(i,:),length(ensemble_members),1);
    cc = mat2cell(color_i,ones(length(ensemble_members),1),repmat(3,1,1));
    set(all_lines(ensemble_members),{'color'},cc);
end

title('Ensembles colored by membership','FontSize',16)

subplot(122)
imagesc(corrMat); caxis([min(corrMat(:)) 1]);
hold on; plot(X,Y,'r','linewidth',2)
title('Sorted correlation matrix','FontSize',16)


