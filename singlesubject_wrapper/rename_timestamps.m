function trial_idx = rename_timestamps(tmsp_folder)
%rename_timestamps Given a timestamp folder (with potentially multiple
%sub-folders of timestamps), get rid of the folders and populate main
%tmsp_folder with all the appropriately named trial timestamp files
%   Detailed explanation goes here

 trial_idx = [];
 
 subfolders = dir(fullfile(tmsp_folder));
 subfolders = {subfolders.name};
 contains_numeric = @(x) ~isnan(str2double(x));
 keep_folder_idx = false(1,length(subfolders));
 for i = 1:length(subfolders)
     keep_folder_idx(i) = contains_numeric(subfolders{i}(1));
 end
 subfolders = subfolders(keep_folder_idx);
 
 % sort the subfolders by earliest trials
 
 start_trials = zeros(1,length(subfolders));
 for i = 1:length(subfolders)
     if ~contains_numeric(subfolders{i}(1:2))
         start_trials(i) = str2double(subfolders{i}(1));
     else
         start_trials(i) = str2double(subfolders{i}(1:2));
     end
 end
 [start_trials,srt] = sort(start_trials,'ascend');
 subfolders = subfolders(srt);
 
for i = 1:length(subfolders)
    
    fullpath = fullfile(tmsp_folder,subfolders{i});
    tmsp_files = dir(fullfile(fullpath,'*.txt'));
    first_trial = start_trials(i);
    last_trial = start_trials(i) + length(tmsp_files) - 1;
    absolute_trial_idx = first_trial:last_trial;
    
    for local_index = 1:length(tmsp_files)
  
        local_Nam = fullfile(fullpath,['image_timestamp_exp_',num2str(local_index),'.txt']);
        movefile(local_Nam, fullfile(tmsp_folder,...
            sprintf('image_timestamp_exp_%s.txt',num2str(absolute_trial_idx(local_index)))));
        
    end
    
    rmdir(fullpath)

    trial_idx = [trial_idx,absolute_trial_idx];
 
 
end

