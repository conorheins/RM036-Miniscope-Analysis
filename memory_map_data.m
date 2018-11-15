%% script to memory map and (optionally) filter data

% load data

if strcmp(already_loaded,'y')
    fprintf('Choose an existing .mat file to analyze or add data to...\n');
    [fnam, fdir] = uigetfile('*.mat');
    dataFile = fullfile(fdir,fnam);
    data = matfile(dataFile,'Writable',true);
%     add_data = input(sprintf('Would you like to add more data to existing file %s ? (y/n, default n)\n',fnam),'s');
    add_data = 'n';
    if strcmp(add_data,'y')
        add_data = 1;
        
        tic
        [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,add_data,trial_based,trial_length,delFrames);
        fprintf('Time taken to add data: %.2f minutes\n',(toc/60))
        
        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origIndices = data.sortedIndices;
        sortedIndices_new = [origIndices;frameIndices];
        data.sortedIndices = sortedIndices_new; clear origIndices sortedIndices_new;
    else
        add_data = 0;
        if strcmp(filter_flag,'y')
            
            chunkSize = 10000;
            
            tic
            data = filter_memmap(data,gaussFilt,chunkSize);
            fprintf('Time taken to filter data: %.2f minutes\n',(toc/60))
            
        end       
    end
elseif strcmp(already_loaded,'n')
    
    fprintf('Choose a folder to determine name of data file...\n');
    dataFile = [uigetdir(),'.mat'];
    [fdir, fnam, temp] = fileparts(dataFile);
    fnam = strcat(fnam,temp); clear temp;
    
%     add_data = input(sprintf('Would you like to add more data after loading new file %s ? (y/n, default n)\n',fnam),'s');
    add_data = 'n';
    
    Y = [];
    if strcmp(filter_flag,'y')
        Yfilt = []; save(dataFile,'Y','Yfilt','-v7.3')
    else
        save(dataFile,'Y','-v7.3');
    end
    
    data = matfile(dataFile,'Writable',true);
    
    tic
    [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,0,trial_based,trial_length,delFrames);
    fprintf('Time taken to load new data: %.2f minutes\n',(toc/60))

    data.Ysiz = Ysiz; data.sortedIndices = frameIndices;
    
    
    if strcmp(add_data,'y')
        add_data = 1;
        
        tic
        [data,Ysiz,frameIndices] = sequence2mat_mod(data,filter_flag,gaussFilt,add_data,trial_based,trial_length,delFrames);
        fprintf('Time taken to add data: %.2f minutes\n',(toc/60))

        origYsiz = data.Ysiz;
        Ysiz_new = [origYsiz(1:2);origYsiz(3)+Ysiz(3)];
        data.Ysiz = Ysiz_new; clear origYsiz Ysiz_new;
        origIndices = data.sortedIndices;
        sortedIndices_new = [origIndices;frameIndices];
        data.sortedIndices = sortedIndices_new; clear origIndices sortedIndices_new;
    else
        add_data = 0; 
    end
end

%get information about the data
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrames = Ysiz(3);    %total number of frames
frameIndices = data.sortedIndices; % trials and frames array 

if strcmp(trial_based,'y')
    fprintf('\nThe data has been mapped to a hard disk. It has %d X %d pixels X %d frames, and is comprised of %d total trials. \nLoading all data requires %.2f GB RAM\n\n', ...
        d1, d2, numFrames,length(unique(frameIndices(:,1))),prod(Ysiz)*8/(2^30));
else
     fprintf('\nThe data has been mapped to a hard disk. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', ...
        d1, d2, numFrames,prod(Ysiz)*8/(2^30));
end