function [cellreg_fnam] = convert_to_cellreg(filenames,session_name,sz,savedir)
% convert_to_cellreg: helper function that converts spatial footprints
% stored in results files (cell arrays of 2 x numSessions, where each
% column is a session, first row is session name, second row is a struct
% with the A (spatial) and C (temporal) matrices

% Input args: 1. filename: string or cell-array of strings with filenames to
%             be converted
%            2. session_name: string or cell_array of strings with session
%            names to convert

% Outputs: 1. cellreg_fnam: string or cell-array of strings with filenames
%               to output of conversion

d1 = sz(1); d2 = sz(2);

if ~exist('savedir','var') || isempty(savedir)
    savedir = pwd;
end

if iscell(filenames)
    
    cellreg_fnam = cell(length(filenames),1);

    for fid = 1:length(filenames)
        currNam = filenames{fid};
        dat = load(currNam);
        ratnam = fieldnames(dat);
        ratnam = ratnam{1};
        dat = dat.(ratnam);
        
        write_dir = fullfile(savedir,[ratnam,'_data']);
        cellreg_fnam{fid} = write_dir;
        
        if ~isdir(write_dir)
            mkdir(write_dir)
        end
        
        if ~exist('session_name','var') || isempty(session_name)
            for sess = 1:size(dat,2)
                spatialz = full(dat{2,sess}.A)';
                K = size(spatialz,1);
                reshaped = reshape(spatialz,K,d1,d2);
                savename = fullfile(write_dir,['spatial_footprints_0',num2str(sess),'.mat']);
                save(savename,'reshaped')
                cellreg_fnam{fid}(1,sess) = ['spatial_footprints_0',num2str(sess),'.mat'];
            end
        else
            for sess = 1:length(session_name)
                idx = find(cellfun(@(x) ~isempty(x), strfind(dat(1,:),session_name{sess})));
                spatialz = full(dat{2,idx}.A)';
                K = size(spatialz,1);
                reshaped = reshape(spatialz,K,d1,d2);
                savename = fullfile(write_dir,['spatial_footprints_0',num2str(idx),'.mat']);
                save(savename,'reshaped');
                cellreg_fnam{fid}(1,sess) = ['spatial_footprints_0',num2str(idx),'.mat'];
            end
        end
    end
else
    
    dat = load(filenames);
    ratnam = fieldnames(dat);
    ratnam = ratnam{1};
    dat = dat.(ratnam);
    
    write_dir = fullfile(savedir,[ratnam,'_data']);
    
    if ~isdir(write_dir)
        mkdir(write_dir)
    end
    
    if ~exist('session_name','var') || isempty(session_name)
        cellreg_fnam = cell(1,size(dat,2));
        for sess = 1:size(dat,2)
            spatialz = full(dat{2,sess}.A)';
            K = size(spatialz,1);
            reshaped = reshape(spatialz,K,d1,d2);
            savename = fullfile(write_dir,['spatial_footprints_0',num2str(sess),'.mat']);
            save(savename,'reshaped')
            cellreg_fnam{sess} = savename;
        end
    else
        cellreg_fnam = cell(1,length(session_name));
        for sess = 1:length(session_name)
            idx = find(cellfun(@(x) ~isempty(x), strfind(dat(1,:),session_name{sess})));
            spatialz = full(dat{2,idx}.A)';
            K = size(spatialz,1);
            reshaped = reshape(spatialz,K,d1,d2);
            savename = fullfile(write_dir,['spatial_footprints_0',num2str(idx),'.mat']);
            save(savename,'reshaped');
            cellreg_fnam{sess} = savename;
        end
    end
end
                
                
                
            
                
            
        
        
