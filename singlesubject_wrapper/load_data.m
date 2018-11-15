function [ratname, C, num_trials, sortedIndices, event_names, event_matrix,empiricalFs] = load_data(behav_fnam, session_name, neural_fnam, tmsp_fdir)
% Load fluorescence data and timestamps, given behavioral data file, neural
% data file, and neuralmapper timestamp files
%   INPUTS: 
%   -behav_fnam: absolute path to .mpc (MED-PC) file with relevant
%       behavioral information for a single session. If there are mutiple
%       .mpc files per session, will read each separately and combine
%       information later
%   -neural_fnam: absolute path to .mat file containing neural data
%   (activity/fluorescence, stored in array 'C')
%   -tmsp_fdir: absolute path to directory containing NeuralMapper
%   timestamps (a series of .txt files)
%   RETURNS:
%   -ratname: string, name of the subject (e.g. 'Rat21')
%   -C: 

[ratname, C] = load_session(neural_fnam,session_name,'temporal');

if iscell(behav_fnam)
    
    [num_trials,trial_starts,cue_starts,cue_ends,HL_on,HL_off,lever_out,lever_in,pellet_tmsp,beambreak_tmsp,...
        ~,~,~,press_timestamps] = combine_multi_medpc(behav_fnam);
        
else     
    
    dataArray = medpc_reader(behav_fnam);
    
    [num_trials,trial_starts,cue_starts,cue_ends,HL_on,HL_off,lever_out,lever_in,pellet_tmsp,beambreak_tmsp,...
        ~,~,~,press_timestamps] = retrieve_tmsp_array(dataArray);
    
end

sortedIndices = [reshape(repmat(1:num_trials,799,1),num_trials*799,1), repmat([1:799]',num_trials,1)];

event_tmsp_array = {lever_out,lever_in,press_timestamps,cue_starts,cue_ends,HL_on,HL_off,beambreak_tmsp,pellet_tmsp};
event_names = {'leverOUT','leverIN','press','cueON','cueOFF','HLON','HLOFF','bbk','pellets'};

empty_event_idx = find(cellfun(@(x) isempty(x), event_tmsp_array));
event_tmsp_array(empty_event_idx) = []; event_names(empty_event_idx) = [];

[event_matrix,empiricalFs,sortedIndices] = create_filter_matrix(sortedIndices,tmsp_fdir,num_trials,trial_starts,event_tmsp_array);
num_trials = length(unique(sortedIndices(:,1)));

end

