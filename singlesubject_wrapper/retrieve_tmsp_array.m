function [num_trials,trial_starts,cue_starts,cue_ends,HL_on,HL_off,lever_out,lever_in,pellet_tmsp,beambreak_tmsp,...
    ITI_timestamps,TTL_starts,trial_timestamps,press_timestamps] = retrieve_tmsp_array( dataArray )
% retrieve_tmsp_array This function takes the cell array of med-pc
% timestamps returned by medpc_reader and returns the behaviorally relevant
% timestamps (e.g. cue on, cue off, HL on, HL off, etc.)

firstCol = dataArray{1};

numerical_data = {};
for col = 2:6
    numerical_data = [numerical_data,dataArray{col}];
end

Cindex = find(cellfun(@(x) strcmp(x,'C:'),firstCol));
Dindex = find(cellfun(@(x) strcmp(x,'D:'),firstCol));
Findex = find(cellfun(@(x) strcmp(x,'F:'),firstCol));
Gindex = find(cellfun(@(x) strcmp(x,'G:'),firstCol));
Kindex = find(cellfun(@(x) strcmp(x,'K:'),firstCol));
Lindex = find(cellfun(@(x) strcmp(x,'L:'),firstCol));
Mindex = find(cellfun(@(x) strcmp(x,'M:'),firstCol));
Pindex = find(cellfun(@(x) strcmp(x,'P:'),firstCol));
Rindex = find(cellfun(@(x) strcmp(x,'R:'),firstCol));
Xindex = find(cellfun(@(x) strcmp(x,'X:'),firstCol));
Yindex = find(cellfun(@(x) strcmp(x,'Y:'),firstCol));
Zindex = find(cellfun(@(x) strcmp(x,'Z:'),firstCol));

if length(Cindex) > 1
    
    num_trials = 0;
    cue_starts = [];
    cue_ends = [];
    HL_on = [];
    HL_off = [];
    lever_out = [];
    lever_in = [];
    beambreak_tmsp = [];
    ITI_timestamps = [];
    TTL_starts = [];
    trial_timestamps = [];
    press_timestamps = [];
    pellet_tmsp = [];
    trial_starts = [];
    
    % temporarily get approximate inter-trial-intervals to add 
    temp = str2double(numerical_data((Xindex(1)+1):(Yindex(1)-1),:))';
    temp = temp(:); temp(isnan(temp))=[];
   
    for i = 1:length(Cindex)
        
        if i == 1
            offset = 0;
        else
            offset = mean(diff(temp(1:2:end))) + trial_starts(end);
        end
        
        D_array = str2double(numerical_data((Dindex(i)+4):(Findex(i)-1),:));
        
        cue_timestamps = D_array(2:3:end,4:5);
        cue_starts_i = cue_timestamps(:,1); cue_starts_i(cue_starts_i==0) = [];
        cue_ends_i = cue_timestamps(:,2); cue_ends_i(cue_ends_i==0) = [];
        cue_starts = [cue_starts;(cue_starts_i + offset)]; cue_ends = [cue_ends;(cue_ends_i + offset)];
        
        HL_timestamps = D_array(3:3:end,2:3);
        HL_on = [HL_on;(HL_timestamps(:,1) + offset)]; HL_off = [HL_off;(HL_timestamps(:,2)+offset)];
        
        lever_timestamps = D_array(1:3:end,4:5);
        lever_out = [lever_out;(lever_timestamps(:,1)+offset)]; lever_in = [lever_in;(lever_timestamps(:,2)+offset)];
        
        beambreak_tmsp_i = str2double(numerical_data((Findex(i)+1):(Gindex(i)-1),:))';
        beambreak_tmsp_i = beambreak_tmsp_i(:); beambreak_tmsp_i(isnan(beambreak_tmsp_i)) = [];
        beambreak_tmsp = [beambreak_tmsp;(beambreak_tmsp_i+offset)]; 
        
        ITI_timestamps_i = str2double(numerical_data((Kindex(i)+1):(Lindex(i)-1),:))';
        ITI_timestamps_i = ITI_timestamps_i(:); ITI_timestamps_i(isnan(ITI_timestamps_i)) = [];
        ITI_timestamps = [ITI_timestamps;(ITI_timestamps_i + offset)];
        
        TTL_starts_i = str2double(numerical_data((Xindex(i)+1):(Yindex(i)-1),:))';
        TTL_starts_i  = TTL_starts_i(:); TTL_starts_i(isnan(TTL_starts_i)) = [];
        TTL_starts = [TTL_starts;(TTL_starts_i+offset)]; 
        
        trial_timestamps_i = str2double(numerical_data((Yindex(i)+1):(Zindex(i)-1),:))';
        trial_timestamps_i = trial_timestamps_i(:); trial_timestamps_i(isnan(trial_timestamps_i)) = [];
        trial_timestamps = [trial_timestamps;(trial_timestamps_i + offset)];
        
        press_timestamps_i = str2double(numerical_data((Lindex(i)+1):(Mindex(i)-1),:))';
        press_timestamps_i = press_timestamps_i(:); press_timestamps_i(isnan(press_timestamps_i)) = [];
        press_timestamps = [press_timestamps;(press_timestamps_i + offset)];
        
        pellet_tmsp_i = str2double(numerical_data((Pindex(i)+1):(Rindex(i)-1),:))';
        pellet_tmsp_i = pellet_tmsp_i(:); pellet_tmsp_i(isnan(pellet_tmsp_i)) = [];
        pellet_tmsp = [pellet_tmsp;(pellet_tmsp_i + offset)];
        
        num_trials = num_trials + str2double(numerical_data{Cindex(i)+1,1});
        trial_starts = [trial_starts;(TTL_starts(1:2:end) + offset)];
        
    end
       
else
    
    
    D_array = str2double(numerical_data((Dindex+4):(Findex-1),:));
    
    cue_timestamps = D_array(2:3:end,4:5);
    cue_starts = cue_timestamps(:,1); cue_ends = cue_timestamps(:,2);
    
    HL_timestamps = D_array(3:3:end,2:3);
    HL_on = HL_timestamps(:,1); HL_off = HL_timestamps(:,2);
    
    lever_timestamps = D_array(1:3:end,4:5);
    lever_out = lever_timestamps(:,1); lever_in = lever_timestamps(:,2);
    
    beambreak_tmsp = str2double(numerical_data((Findex+1):(Gindex-1),:))';
    beambreak_tmsp = beambreak_tmsp(:); beambreak_tmsp(isnan(beambreak_tmsp))=[];
    
    ITI_timestamps = str2double(numerical_data((Kindex+1):(Lindex-1),:))';
    ITI_timestamps = ITI_timestamps(:); ITI_timestamps(isnan(ITI_timestamps))=[];
    
    TTL_starts = str2double(numerical_data((Xindex+1):(Yindex-1),:))';
    TTL_starts = TTL_starts(:); TTL_starts(isnan(TTL_starts))=[];
    
    trial_timestamps = str2double(numerical_data((Yindex+1):(Zindex-1),:))';
    trial_timestamps = trial_timestamps(:); trial_timestamps(isnan(trial_timestamps)) = [];
    
    press_timestamps = str2double(numerical_data((Lindex+1):(Mindex-1),:))';
    press_timestamps = press_timestamps(:); press_timestamps(isnan(press_timestamps)) = [];
    
    pellet_tmsp = str2double(numerical_data((Pindex+1):(Rindex-1),:))';
    pellet_tmsp = pellet_tmsp(:); pellet_tmsp(isnan(pellet_tmsp)) = [];
    
    num_trials = str2double(numerical_data{Cindex(1)+1,1});
    trial_starts = TTL_starts(1:2:end);
    
end




end

