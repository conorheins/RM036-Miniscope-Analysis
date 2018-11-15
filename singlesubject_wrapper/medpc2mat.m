clear all; clc
%% import med pc data
[filename,filedirectory]=uigetfile('*');
cd(filedirectory)

delimiter = ' ';
startRow = 1;

formatSpec = '%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');

textscan(fileID, '%[^\n\r]', startRow, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);

fclose(fileID);

dataArray(7:end) = [];

%%  isolate specific time stamps from different numerical arrays (e.g. K, X, Y, etc.)

firstCol = dataArray{1};

numerical_data = {};
for col = 2:6
    numerical_data = [numerical_data,dataArray{col}];
end

nameindex = find(cellfun(@(x) strcmp(x,'Subject:'),firstCol));
Cindex = find(cellfun(@(x) strcmp(x,'C:'),firstCol));
Dindex = find(cellfun(@(x) strcmp(x,'D:'),firstCol));
Kindex = find(cellfun(@(x) strcmp(x,'K:'),firstCol));
Lindex = find(cellfun(@(x) strcmp(x,'L:'),firstCol));
Mindex = find(cellfun(@(x) strcmp(x,'M:'),firstCol));
Xindex = find(cellfun(@(x) strcmp(x,'X:'),firstCol));
Yindex = find(cellfun(@(x) strcmp(x,'Y:'),firstCol));
Zindex = find(cellfun(@(x) strcmp(x,'Z:'),firstCol));

ITI_timestamps = str2double(numerical_data((Kindex+1):(Lindex-1),:))';
ITI_timestamps = ITI_timestamps(:); ITI_timestamps(isnan(ITI_timestamps))=[];

TTL_starts = str2double(numerical_data((Xindex+1):(Yindex-1),:))';
TTL_starts = TTL_starts(:); TTL_starts(isnan(TTL_starts))=[];

trial_timestamps = str2double(numerical_data((Yindex+1):(Zindex-1),:))';
trial_timestamps = trial_timestamps(:); trial_timestamps(isnan(trial_timestamps)) = [];

press_timestamps = str2double(numerical_data((Lindex+1):(Mindex-1),:))';
press_timestamps = press_timestamps(:); press_timestamps(isnan(press_timestamps)) = [];


%% create struct for the total session, including all trials as sub-structs within it

session_struct = struct;
num_trials = str2double(numerical_data{Cindex+1,1});

trial_starts = TTL_starts(1:2:end);

for i = 1:num_trials
    temp_trial_struct = struct('frame_tmsp',[],'Press_tmsp',[]);
    frame_grabs = trial_starts(i);
    neural_mapper_name = ['image_timestamp_exp_',num2str(i),'.txt'];
    temp = read_in_timestamps(neural_mapper_name);
    frame_grabs = [frame_grabs; trial_starts(i) + temp];
    temp_trial_struct.frame_tmsp = frame_grabs;
    start_t = min(frame_grabs); end_t = max(frame_grabs);
    temp_trial_struct.Press_tmsp = press_timestamps(press_timestamps > start_t & press_timestamps < end_t);
    session_struct.(['Trial_',num2str(i)]) = temp_trial_struct;
end
    



    
    
