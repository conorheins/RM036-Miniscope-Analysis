function [ trial_idx ] = obtain_trial_idx( tmsp_folder )
%obtain_trial_idx Summary of this function goes here
%   Detailed explanation goes here


tmsp_files = dir(fullfile(tmsp_folder,'*.txt'));
tmsp_files = {tmsp_files.name};

trial_idx = zeros(1,length(tmsp_files));


for i = 1:length(tmsp_files)
    
    first_split = strsplit(tmsp_files{i},'_');
    second_split = strsplit(first_split{end},'.');
    trial_idx(i) = str2double(second_split{1});
    
end

trial_idx = sort(trial_idx,'ascend');


end

