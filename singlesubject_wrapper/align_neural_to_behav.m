function [ event_locked,relative_event_tmsp ] = align_neural_to_behav(event_names,before_align,after_align,event_matrix,data2use)
%align_neural_to_behav Aligns neural data to behavioral events and stores
%in a cell array of event-locked matrices. Uses event_matrix to extract
%event timestamps

num_events = length(event_names);
event_locked = cell(num_events,1);

minimum_tmsp_ISI = max(before_align,after_align);
relative_event_tmsp = before_align + 1;

for event_i = 1:num_events
    fprintf('Now extracting calcium aligned to event: %s\n',event_names{event_i});
    
    event_times = find(event_matrix(:,event_i) == 1);
    
    event_locked_calcium = [];
    
    if length(event_times) > 1
        
        curr_tmsp = event_times(1);
        tmsp = 2;
        
        throw_out_idx = false(length(event_times),1);
        
        for tmsp = tmsp:length(event_times)
            if (event_times(tmsp) - curr_tmsp) <= minimum_tmsp_ISI
                throw_out_idx(tmsp) = true;
            elseif (event_times(tmsp) - curr_tmsp) > minimum_tmsp_ISI
                curr_tmsp = event_times(tmsp);
                tmsp = tmsp + 1;
            end
        end
        
        fprintf('Threw out a total of %d timestamps from event %s\n',sum(throw_out_idx),event_names{event_i});
        event_times = event_times(~throw_out_idx);
        
    end
    
    if strcmp(event_names{event_i},'HLOFF') % in case of HLOFF, where last time-stamp of trial comes within within ~11.7 seconds of trial_end, use shorter 'after_align' 
        after_align_HLOFF = 100;
        for event_time_i = 1:length(event_times)
            event_tmsp = event_times(event_time_i);
            if ~( (event_tmsp + after_align_HLOFF) > size(data2use,2))
                event_locked_calcium = cat(3,event_locked_calcium,data2use(:, (event_tmsp - before_align):(event_tmsp + after_align_HLOFF)));
            end
        end
    else
        for event_time_i = 1:length(event_times)
            event_tmsp = event_times(event_time_i);
            event_locked_calcium = cat(3,event_locked_calcium,data2use(:, (event_tmsp - before_align):(event_tmsp + after_align)));
        end
    end
    event_locked{event_i} = event_locked_calcium;
    fprintf('Done extracting calcium aligned to event: %s\n',event_names{event_i});

end

end

