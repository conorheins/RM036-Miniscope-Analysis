%% set paths

clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));
addpath(genpath('helper_func'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);


%%

leverOUT_idx = 250:300; % indices within a single trial to investigate

RatIDs = [9,10,11,21,23];
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

rat_folder = 'rat_processed';

registration_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/CellReg-1.3.4/data';
sess_range = 3:7;

all_rats_sessVects = cell(1,length(RatIDs));

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    cell_map = cell_map(:,sess_range);
   
    all_sessions_badneur = cell(1,length(sess_range));
    
    for jj = 1:length(sess_range)
        
        sess_id = session_names{sess_range(jj)};
                
        rat_dir = fullfile('rats_processed',sprintf('Rat%d',rat_id));
        fnam = fullfile(rat_dir,sprintf('Rat%d_%s.mat',rat_id,sess_id));
        load(fnam);
        
        all_sessions_badneur{jj} = Sess_object.bad_neurons;
        
    end
    
    [cell_map,~] = calculate_global_excludeIDX(cell_map,all_sessions_badneur);
    common_neurons = find(sum(cell_map > 0,2) == length(sess_range));
    
    num_common = length(common_neurons);
    
    temp_rat_data = cell(length(sess_range),1);
    
    session_vectors = zeros(num_common,length(sess_range),2);
    
    for jj = 1:length(sess_range)
        
        sess_id = session_names{sess_range(jj)};
        
        rat_dir = fullfile('rats_processed',sprintf('Rat%d',rat_id));
        fnam = fullfile(rat_dir,sprintf('Rat%d_%s.mat',rat_id,sess_id));
        load(fnam);
        
        common_neuron_tempidx = cell_map(common_neurons,jj);

        reshaped_events = permute(reshape(full(Sess_object.event_matrix),[],Sess_object.num_trials,length(Sess_object.event_names)),[1,3,2]);
        press_idx = find(cellfun(@(x) ~isempty(x),strfind(Sess_object.event_names,'press')));
        noresponse_trial_idx = sum(squeeze(reshaped_events(:,press_idx,:)),1) == 0;
        
        data = Sess_object.spikes_conv(common_neuron_tempidx,:);
        
        data = bsxfun(@rdivide,data - min(data,[],2),max(data,[],2) - min(data,[],2));
        reshaped_data = reshape(data,num_common,799,Sess_object.num_trials);
        
        if length(find(noresponse_trial_idx)) < 20
            fprintf('Not enough no-response trials to compute no-response vectors separately, only doing response trials\n')
            session_vectors(:,jj,1) = mean(mean(reshaped_data(:,leverOUT_idx,~noresponse_trial_idx),3),2);
        else
            session_vectors(:,jj,1) = mean(mean(reshaped_data(:,leverOUT_idx,~noresponse_trial_idx),3),2);
            session_vectors(:,jj,2) = mean(mean(reshaped_data(:,leverOUT_idx,noresponse_trial_idx),3),2);
        end
        
    end
    
    all_rats_sessVects{ii} = session_vectors;
    
end

%% visualize results

for ii = 1:length(RatIDs)
    
    figure('Position',[400 400 1000 1000])    
    subplot(221)
    imagesc(all_rats_sessVects{ii}(:,:,1))
    title('Population vectors aligned to Lever-Entry, R Trials')
    xticklabels(session_names(sess_range));
    ylabel('Neuron #')
     
    subplot(222)
    imagesc(all_rats_sessVects{ii}(:,:,2))
    title('Population vectors aligned to Lever-Entry, NR Trials')
    xticklabels(session_names(sess_range));
    ylabel('Neuron #')
    
    subplot(223)
    imagesc(1-pdist2(all_rats_sessVects{ii}(:,:,1)',all_rats_sessVects{ii}(:,:,1)','correlation'))
    colorbar;
    title('Pattern similarity across sessions: R Trials')
    xticklabels(session_names(sess_range));
    yticklabels(session_names(sess_range));
    
    subplot(224)
    imagesc(1-pdist2(all_rats_sessVects{ii}(:,:,2)',all_rats_sessVects{ii}(:,:,2)','correlation'))
    colorbar;
    title('Pattern similarity across sessions: NR Trials')
    xticklabels(session_names(sess_range));
    yticklabels(session_names(sess_range));
    
    pause; close gcf;
end
    
    
        
        
        
        
    
    




