function [ cell_to_index_map ] = load_cell_map(registration_folder,RatID)
%load_cell_map Load the registration data (output of CellReg) for a
%particular rat
%   Detailed explanation goes here

subfolders = dir(registration_folder);
subfolders = {subfolders(:).name};

rat_folder = subfolders{cellfun(@(x) ~isempty(x),strfind(subfolders,num2str(RatID)))};
reg_data = dir([registration_folder,filesep,rat_folder,filesep,'Results']);
reg_data = {reg_data(:).name};

reg_file = reg_data{cellfun(@(x) ~isempty(x),strfind(reg_data,'cellRegistered'))};

reg_data = load([registration_folder,filesep,rat_folder,filesep,'Results',filesep,reg_file]);

cell_to_index_map = reg_data.cell_registered_struct.cell_to_index_map;

end

