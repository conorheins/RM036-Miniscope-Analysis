function [sessions_array] = make_sessions_array(subjectArray,data_type)
% takes cell array of individual SubjectObj class objects (one for a given
% data-collection session, e.g. SA1) and puts their trace/behavioral data into 
% a single cell array

numSess = length(subjectArray);
sessions_array = cell(2,numSess);

if ~exist('data_type','var') || isempty(data_type)
    data_type = 'spikes_denoised';
end

for i = 1:numSess
    sessions_array{1,i} = subjectArray{i}.(data_type);
    sessions_array{2,i} = {subjectArray{i}.event_names, subjectArray{i}.event_matrix};
    
end

end