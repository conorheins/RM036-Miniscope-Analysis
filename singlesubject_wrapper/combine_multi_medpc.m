function [num_trials_all,trial_starts_all,cue_starts_all,cue_ends_all,HL_on_all,HL_off_all,lever_out_all,lever_in_all,pellet_tmsp_all,beambreak_tmsp_all,...
    ITI_timestamps_all,TTL_starts_all,trial_timestamps_all,press_timestamps_all] = combine_multi_medpc( behav_fnam )
%combine_multi_medpc Takes combination of med-pc files and spits out
%accumulated timestamps
%   Detailed explanation goes here

num_trials_all = 0;
trial_starts_all = [];
cue_starts_all = [];
cue_ends_all = [];
HL_on_all = [];
HL_off_all = [];
lever_out_all = [];
lever_in_all = [];
pellet_tmsp_all = [];
beambreak_tmsp_all = [];
ITI_timestamps_all = [];
TTL_starts_all = [];
trial_timestamps_all = [];
press_timestamps_all =  [];

for file_id = 1:length(behav_fnam)
    
    dataArray = medpc_reader(behav_fnam{file_id});
    
    [num_trials,trial_starts,cue_starts,cue_ends,HL_on,HL_off,lever_out,lever_in,pellet_tmsp,beambreak_tmsp,...
        ITI_timestamps,TTL_starts,trial_timestamps,press_timestamps] = retrieve_tmsp_array(dataArray);
    
    if file_id == 1
        offset = 0;
    else
        offset = trial_starts_all(end) + mean(diff(trial_starts_all));
    end
    
    num_trials_all = num_trials_all + num_trials;
    trial_starts_all = [trial_starts_all; (trial_starts + offset)];
    cue_starts_all = [cue_starts_all; (cue_starts + offset)];
    cue_ends_all = [cue_ends_all; (cue_ends + offset)];
    HL_on_all = [HL_on_all; (HL_on + offset)];
    HL_off_all = [HL_off_all;(HL_off + offset)];
    lever_out_all = [lever_out_all;(lever_out + offset)];
    lever_in_all = [lever_in_all;(lever_in + offset)];
    pellet_tmsp_all = [pellet_tmsp_all;(pellet_tmsp + offset)];
    beambreak_tmsp_all = [beambreak_tmsp_all;(beambreak_tmsp + offset)];
    ITI_timestamps_all = [ITI_timestamps_all;(ITI_timestamps + offset)];
    TTL_starts_all = [TTL_starts_all;(TTL_starts + offset)];
    trial_timestamps_all = [trial_timestamps_all;(trial_timestamps + offset)];
    press_timestamps_all = [press_timestamps_all;(press_timestamps + offset)];
     
    
end

end

