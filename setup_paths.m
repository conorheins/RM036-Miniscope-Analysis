%% add path

main_dir = fileparts(which('setup_paths.m'));
addpath(fullfile(main_dir, 'ca_source_extraction'));
addpath(genpath(fullfile(main_dir, 'ca_source_extraction', 'utilities')));
addpath(genpath(fullfile(main_dir, 'ca_source_extraction', 'endoscope')));
addpath(genpath(fullfile(pwd,'NoRMCorre-master')))
addpath(genpath(fullfile(pwd,'usefulfunctions')));
addpath(genpath(fullfile(pwd,'chronux_2_12')));


%% deconvolution 
run(fullfile(main_dir, 'deconvolveCa', 'setup_nocvx.m'));
