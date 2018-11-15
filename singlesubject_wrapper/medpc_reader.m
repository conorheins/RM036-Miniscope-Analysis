function [ dataArray ] = medpc_reader( filename )

%% import med pc data

delimiter = ' ';
startRow = 1;

formatSpec = '%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');

textscan(fileID, '%[^\n\r]', startRow, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);

fclose(fileID);

dataArray(7:end) = [];


end

