
ratIDs = 8;

data_directory = 'NeuralData/PreMerge/';

for i = 1:length(ratIDs)
    fnam = fullfile(data_directory,sprintf('Rat%s/Rat%s_results_mrg.mat',num2str(ratIDs(i)),num2str(ratIDs(i))));
    currentRat = load(fnam);
    ratname = fieldnames(currentRat);
    ratname = ratname{1};
    currentRat = currentRat.(ratname);
    for session = 1:size(currentRat,2)
        spatialz = currentRat{2,session}.A;
        temporalz = currentRat{2,session}.C;
        if ~isempty(find(sum(spatialz,1)==0))
            empty_comps = find(sum(spatialz,1)==0);
            fprintf('The following components are empty: %s\n',num2str(empty_comps))
            spatialz(:,empty_comps) = [];
            temporalz(empty_comps,:) = [];
        else
            fprintf('All components have at least 1 non-zero pixel\n')
        end
        spatialz_norm = bsxfun(@rdivide,bsxfun(@minus,spatialz,min(spatialz,[],1)),max(spatialz,[],1) - min(spatialz,[],1));
        projection = sum(spatialz_norm,2);
        figure; imagesc(reshape(projection,400,400));
        include_flag = input('Would you like to select a region of the FOV to include (y/n)?\n','s');
        if strcmp(include_flag,'y')
            includeBox = getrect(gcf);
            [~,lin_coords] = max(spatialz_norm,[],1);
            [y_coords, x_coords] = ind2sub([400 400],lin_coords);
            ind_keep = find(x_coords > includeBox(1) & x_coords < includeBox(1) + includeBox(3) ...
                & y_coords > includeBox(2) & y_coords < includeBox(2) + includeBox(4));
            new_components = spatialz(:,ind_keep);
            temporalz = temporalz(ind_keep,:);
        else
            new_components = spatialz;
        end
        currentRat{2,session}.A = new_components;
        currentRat{2,session}.C = temporalz;
        close(gcf);
    end
    eval(sprintf('%s = currentRat',ratname));
    writenam = fullfile(data_directory,sprintf('Rat%s/Rat%s_results_mrg_trimmed.mat',num2str(ratIDs(i)),num2str(ratIDs(i))));
    save(writenam,sprintf('%s',ratname));
end



%confirm cell counts across A and C

for i = 1:length(ratIDs)
    ratname = ['Rat',num2str(ratIDs(i))]; 
    eval(sprintf('currentRat = %s',ratname)); 
    for sess = 1:7 
        disp(['Rat ',num2str(ratIDs(i)),':',num2str(size(currentRat{2,sess}.A,2)),' and ',num2str(size(currentRat{2,sess}.C,1))])
    end
end
        
   
            
            

