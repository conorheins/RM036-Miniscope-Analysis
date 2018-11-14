%%

clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));
addpath(genpath('helper_func'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);

%% individual cell check

read_directory = 'NeuralData/PreMerge/';
write_directory = 'NeuralData/PostMerge/';

ratIDs = [8];

all_training_cases = [];
all_labels = [];

for i = 1:length(ratIDs)
    fnam = fullfile(read_directory,sprintf('Rat%s/Rat%s_results_mrg_trimmed.mat',num2str(ratIDs(i)),num2str(ratIDs(i))));
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
        spatialz_norm = bsxfun(@rdivide,spatialz - min(spatialz,[],1),max(spatialz,[],1) - min(spatialz,[],1));
        projection = sum(spatialz_norm,2);
        figure; imagesc(reshape(projection,400,400));
        include_flag = input('Would you like to select a region of the FOV to include (y/n)?\n','s');
        if strcmp(include_flag,'y')
            includeBox = getrect(gcf);
            [~,lin_coords] = max(spatialz_norm,[],1);
            [y_coords, x_coords] = ind2sub([400 400],lin_coords);
            ind_keep = find(x_coords > includeBox(1) & x_coords < includeBox(1) + includeBox(3) ...
                & y_coords > includeBox(2) & y_coords < includeBox(2) + includeBox(4));
            spatialz = spatialz(:,ind_keep);
            temporalz = temporalz(ind_keep,:);
        end
        
        close gcf;
        
        check_flag = input('Would you like to individually check each cell component (y/n)?\n','s');
        if strcmp(check_flag,'y')
            num_cells = size(spatialz,2);
            del_inds = false(1,num_cells);
            training_cases = zeros(31^2,num_cells);
            labels = zeros(num_cells,1);

            figure('Position',[100 1400 1300 500]);
            for neur = 1:num_cells
                [~,lin_coords] = max(spatialz(:,neur));
                [y_coords, x_coords] = ind2sub([400 400],lin_coords);
                
                subplot(121);
                img_data = reshape(spatialz(:,neur),400,400);
                
                if or(or(y_coords - 15 <= 0, y_coords+15 > 400),or(x_coords - 15 <= 0, x_coords+15 > 400))
                    del_inds(neur) = true;
                    labels(neur) = NaN;
                else
                    img_data = img_data((y_coords-15:y_coords+15),(x_coords-15:x_coords+15));
                    imagesc(img_data);
                    
                    training_cases(:,neur) = reshape(img_data,31^2,1);
                    
                    subplot(122);
                    plot(temporalz(neur,10000:20000));
                    
                    keep_flag = input('Would you like to keep this component (y/n)? \n','s');
                    
                    if strcmp(keep_flag,'n')
                        del_inds(neur) = true;
                        labels(neur) = 0;
                    else
                        labels(neur) = 1;
                    end
                end
            end   
            new_spatialz = spatialz; 
            new_spatialz(:,del_inds) = [];
            new_temporalz = temporalz;
            new_temporalz(del_inds,:) = [];
        else
            new_spatialz = spatialz;
            new_temporalz = temporalz;
        end
        
        all_training_cases = [all_training_cases,training_cases];
        all_labels = [all_labels; labels];
        
        currentRat{2,session}.A = new_spatialz;
        currentRat{2,session}.C = new_temporalz;
        
        close gcf;
        
    end
    
    
    eval(sprintf('%s = currentRat',ratname));
    
    save(fullfile(write_directory,['Rat',num2str(ratIDs(i)),'_results_final.mat']),sprintf('%s',ratname));
end

