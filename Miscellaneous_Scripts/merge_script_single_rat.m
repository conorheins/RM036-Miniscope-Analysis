
tic

data_directory = 'NeuralData/PreMerge/';
load(fullfile(data_directory,'Rat8/Rat8_results.mat'))
numsessions = size(Rat8,2);

thresholds = [0.8 0.7];
sz = [400, 400];
merge_type = 'both';
display_flag = false;

for session = 1:numsessions
    sess_name = Rat8{1,session};
    A = Rat8{2,session}.A;
    C = Rat8{2,session}.C;
    [A_merge, C_merge] = conor_neuron_merge(A, C, sz, merge_type, thresholds, display_flag);
    Rat8{2,session}.A = A_merge;
    Rat8{2,session}.C = C_merge;
    fprintf('Number of cells in %s : %d\n',Rat8{1,session},size(A_merge,2));
    fprintf('-----------------\n')
end

fprintf('Time taken to merge all neurons from Rat8: %.2f minutes\n',(toc/60));

save(fullfile(data_directory,'Rat8/Rat8_results_mrg.mat'),'Rat8');

%%

filenames = {'Rat9_results_mrg.mat','Rat10_results_mrg.mat','Rat11_results_mrg.mat','Rat15_results_mrg.mat','Rat21_results_mrg.mat','Rat23_results_mrg.mat'};

cellcounts = zeros(length(filenames),7);

for file = 1:length(filenames)
    currentRat = load(filenames{file});
    ratname = fieldnames(currentRat);
    ratname = ratname{1};
    currentRat = currentRat.(ratname);
    for session = 1:size(currentRat,2)
        cellcounts(file,session) = size(currentRat{2,session}.A,2);
    end
end
    
