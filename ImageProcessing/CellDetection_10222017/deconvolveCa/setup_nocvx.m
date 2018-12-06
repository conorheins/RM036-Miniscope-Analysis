oasis_folder = fileparts(mfilename('fullpath')); 
addpath(sprintf('%s', oasis_folder)); 
addpath(sprintf('%s%sfunctions', oasis_folder, filesep)); 
addpath(sprintf('%s%soasis', oasis_folder, filesep)); 
addpath(sprintf('%s%soasis_kernel', oasis_folder, filesep)); 
addpath(sprintf('%s%sMCMC', oasis_folder, filesep)); 
addpath(sprintf('%s%sMCMC%sutilities', oasis_folder, filesep, filesep)); 

%% install convex optimization solvers
optimization_folder = sprintf('%s%soptimization', oasis_folder, filesep); 
if ~exist(optimization_folder, 'dir')
    mkdir(optimization_folder);
end
