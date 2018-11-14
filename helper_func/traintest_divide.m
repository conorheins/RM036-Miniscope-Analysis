function [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,neur_idx,train_proportion)

if ~exist('train_test_proportion','var') || isempty(train_proportion)
    train_proportion = 0.8;
end

if ~exist('data_type','var') || isempty(data_type)
    data_type = 'spikes_conv';
end

data2use = Sess_object.(data_type);
if and(exist('neur_idx','var'),~isempty(neur_idx))
    data2use = data2use(neur_idx,:);
end
    
event_ids = zeros(1,length(event_names));

for event_i = 1:length(event_names)
    event_ids(event_i) = find(strcmp(event_names{event_i},Sess_object.event_names));
end

%% get 'null' data epochs, if desired
if strcmp(use_null_label,'n')
    fullY = [];
    fullX = [];
    label_names = event_names;
elseif strcmp(use_null_label,'y')
    window = (-window_size):(window_size); % temporal window around which to exclude timepoints from being included in 'null' (reference) examples
    event_times = sum(Sess_object.event_matrix(:,event_ids),2) > 0;
    num_categories = length(event_names) + 1; %add +1 to include one 'null events' category
    null_tmsp = conv(double(full(event_times)),ones(length(window),1),'same') == 0; % available timestamps are those that aren't within 1 second (pre or post) event timestamps
    null_tmsp = find(null_tmsp); % get the indices of the available timestamps for 'null' labels
    null_idx = null_tmsp(randperm(length(null_tmsp),1300)); % randomly draw 1300 indices from the list
    [sorted_null_idx, sort_null_idx_srt] = sort(null_idx,'ascend');
    null_idx(sort_null_idx_srt(diff(sorted_null_idx) <= window_size)) = []; % should end up with somewhere between 1000 and 1100 timestamps
    null_idx(null_idx < (window_size+1)) = [];
    nullY = num_categories*ones(length(null_idx),1); 
    nullX = zeros(length(null_idx),size(data2use,1));
    for null_i = 1:length(null_idx)
%         nullX(null_i,:) = mean(data2use(:,(null_idx(null_i)-window_size):null_idx(null_i)),2)';
        nullX(null_i,:) = mean(data2use(:,null_idx(null_i):(null_idx(null_i)+window_size)),2)';
    end
    fullY = nullY;
    fullX = nullX;
    label_names = [event_names,'NULL'];
end

%% assemble labels and corresponding data
cat_label = 1;
for event_i = 1:length(event_names)
    cat_idx = find(Sess_object.event_matrix(:,event_ids(event_i)));
    cat_idx(find(diff(cat_idx) < 10) + 1) = [];
    catY = cat_label*ones(length(cat_idx),1);
    catX = zeros(length(cat_idx),size(data2use,1));
    for y_i = 1:length(cat_idx)
%         catX(y_i,:) = mean(data2use(:,(cat_idx(y_i)-10):cat_idx(y_i)),2)';
        catX(y_i,:) = mean(data2use(:,cat_idx(y_i):(cat_idx(y_i)+window_size)),2)';
    end
    fullY = [fullY;catY];
    fullX = [fullX;catX];
    cat_label = cat_label + 1;
end


%% find label with minimum number of examples, and trim all other labels/examples to fit the smallest one 
min_examples = length(find(fullY == fullY(1)));
label_num_list = unique(fullY);
for lab_i = 1:length(label_num_list)
    label_num = label_num_list(lab_i);
    num_examples = length(find(fullY == label_num));
    if num_examples < min_examples
        min_examples = num_examples;
    end
end

trimY = [];
trimX = [];

for lab_i = 1:length(label_num_list)
    label_num = label_num_list(lab_i);
    example_idx = find(fullY == label_num);
    trimY = [trimY;fullY(example_idx(1:min_examples))];
    trimX = [trimX;fullX(example_idx(1:min_examples),:)];
end

random_idx = randperm(size(trimY,1));
trimY = trimY(random_idx);
trimX = trimX(random_idx,:);

%% separate into train and test sets

trainX = [];
trainY = [];

testX = [];
testY = [];

for lab_i = 1:length(label_num_list)
    
    example_idx = find(trimY == label_num_list(lab_i));
    train_idx = 1:floor(train_proportion*length(example_idx));
    trainX = [trainX;trimX(example_idx(train_idx),:)];
    trainY = [trainY;trimY(example_idx(train_idx))];
    
    test_idx = (ceil(train_proportion*length(example_idx))+1):length(example_idx);
    testX = [testX;trimX(example_idx(test_idx),:)];
    testY = [testY;trimY(example_idx(test_idx))];
    
end

random_idx = randperm(size(trainX,1));
trainX = trainX(random_idx,:);
trainY = trainY(random_idx);

random_idx = randperm(size(testX,1));
testX = testX(random_idx,:);
testY = testY(random_idx);


end