function [ ratname, data ] = load_session( filename, which_session, which_components )
% load_session Extracts spatial/temporal components or both spatial and temporal components from
% specific session in a neural data file

temp_struct = load(filename);
ratname = fieldnames(temp_struct);
ratname = ratname{1};

all_sessions = temp_struct.(ratname);

% pull out SA1, temporal components
col_index = find(strcmp(all_sessions(1,:),which_session));

switch lower(which_components)
    case 'spatial'
        data = all_sessions{2,col_index}.A;
    case 'temporal'
        data = all_sessions{2,col_index}.C;
    case 'both'
        data = all_sessions{2,col_index};
end

end

