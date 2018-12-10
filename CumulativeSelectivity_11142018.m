%% Script 11.14.2018
% edited 11.15.2018
% edited 12.01.2018


%% paths
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

registration_folder = fullfile(base_directory,'CellReg-1.3.4/data');


%% Read in data (already saved as 'subject-objects' for each session), rat by rat

%% parameters for plotting/rats to use/etc.
RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

event_name = {'leverOUT','leverIN','HLON','HLOFF'}; 
data_type = 'spikes_conv'; 

%%  main-loop through rats/neurons/sessions

daily_proportion_sig_neurons = zeros(length(RatIDs),length(session_names),length(event_name));

all_sig_arrays = cell(2,length(RatIDs),length(event_name));

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));
    exist_sessions = {};
    all_data_files = dir([rat_folder,filesep,'*.mat']);
    all_data_files = {all_data_files(:).name};
    for dat_i = 1:length(all_data_files)
        num_remove = length(sprintf('Rat%d_',rat_id))+1;
        exist_sessions{dat_i} = all_data_files{dat_i}(num_remove:end-4);
    end
    
    sess_ids = [];
    for kk = 1:length(session_names)
        if sum(strcmp(exist_sessions,session_names{kk})) > 0
            sess_ids = [sess_ids,kk];
        end
    end
    
    
    sig_map = zeros([size(cell_map),length(event_name)]);
    
    for kk = 1:length(sess_ids)
        
        sess_name = session_names{sess_ids(kk)};
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        Sess_object.lock2events(data_type,52,52);
        Sess_object.compute_selectivity('permutation',[52 52],1000,0.05);
        
        for event_i = 1:length(event_name)
        
            sig_matrix = Sess_object.stats_results.SigMatrix(:,strcmp(Sess_object.event_names,event_name{event_i}));
            pos_mod_matrix = Sess_object.stats_results.Modulation_direction(:,strcmp(Sess_object.event_names,event_name{event_i}));
            
            combined = sig_matrix == 1 & and(~isnan(pos_mod_matrix),pos_mod_matrix == 1);
            
            sig_map(find(cell_map(:,kk)),kk,event_i) = combined(cell_map(find(cell_map(:,kk)),kk));
            
        end
        
    end
   
    for event_i = 1:length(event_name)
        daily_proportion_sig_neurons(ii,sess_ids,event_i) = sum(squeeze(sig_map(:,:,event_i)),1)./sum(cell_map>0,1);
        all_sig_arrays{1,ii,event_i} = squeeze(sig_map(:,:,event_i));
        all_sig_arrays{2,ii,event_i} = cell_map;
    end
    
end


%% some graphics

% 1. bar graph of proportion of neurons (relative to that day's count) that
% are selective for the event-of-interest on that session

% figure(1)
% bar(sum(sig_map,1)./sum(cell_map > 0,1));
% 
% 
% % 2. Histogram of how many days neurons were selective for leverOUT
% 
% figure(2)
% hist(sum(sig_map,2),15);
% 
% % 3. Jaccard (dis)similarity matrix between selectivity vectors on
% % all sessions
% 
% figure(3);
% imagesc(1-pdist2(sig_map',sig_map','jaccard'))

% 4. visualize significant neurons (modulated by event_name) over time

figure();

lineProps.col = {[0 0 1];[0 0.25 0.6];[0.6 0.15 0];[0.9 0.05 0]};
mseb(1:7,squeeze(mean(daily_proportion_sig_neurons,1))',squeeze(std(daily_proportion_sig_neurons,0,1)./sqrt(size(daily_proportion_sig_neurons,1)))',lineProps)
xticklabels(session_names)
set(gca,'FontSize',16);

xlabel('Session','FontSize',20)
ylabel('Proportion of neurons detected per session','FontSize',20)
title('Session-relative proportion of neurons selective for different behavioral events','FontSize',24)

legend(event_name,'FontSize',18)

%% Plot average traces of that were lever-selective neurons on Extinction 4, on different sessions

event_nam_i = 'leverOUT';
initial_sess = 'Ext4';

for ii = 1:length(RatIDs)

    rat_id = RatIDs(ii);
    rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));
          
    exist_sessions = {};
    all_data_files = dir([rat_folder,filesep,'*.mat']);
    all_data_files = {all_data_files(:).name};
    for dat_i = 1:length(all_data_files)
        num_remove = length(sprintf('Rat%d_',rat_id))+1;
        exist_sessions{dat_i} = all_data_files{dat_i}(num_remove:end-4);
    end
    
    sess_ids = 1:7; % but figure out this variable 'sess_ids' using available data from given rat
    for kk = 1:length(session_names)
        if any(strcmp(exist_sessions,session_names{kk}))
            sess_ids = [sess_ids,kk];
        end
    end
    
    window = 230:314;
    x_ax = window./10.49;
    
    neurons2track = find(all_sig_arrays{2,strcmp(event_name,event_nam_i)}(:,strcmp(initial_sess,session_names(sess_ids))));
    
    figure;
    
    for kk = 3:6
        
        sess_specific_idx = all_sig_arrays{2,ii,strcmp(event_name,event_nam_i)}(neurons2track,kk);
        sess_specific_idx = sess_specific_idx(sess_specific_idx ~= 0);
        sess_name = session_names{sess_ids(kk)};
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        R_trial_bool = logical(sum(full(reshape(Sess_object.event_matrix(:,find(strcmp('cueON',Sess_object.event_names))),799,Sess_object.num_trials)),1));
        NR_trial_bool = ~R_trial_bool;
        
        reshaped = reshape(Sess_object.spikes_conv(sess_specific_idx,:),[length(sess_specific_idx), 799,Sess_object.num_trials]);
        
        R_trial_avg = mean(mean(reshaped(:,:,R_trial_bool),3),1);
        R_trial_std = std(mean(reshaped(:,:,R_trial_bool),3),0,1)./sqrt(length(sess_specific_idx));
        
        NR_trial_avg = mean(mean(reshaped(:,:,NR_trial_bool),3),1);
        NR_trial_std = std(mean(reshaped(:,:,NR_trial_bool),3),0,1)./sqrt(length(sess_specific_idx));
        
        subplot(4,1,kk-2)
        mseb(x_ax,[R_trial_avg(window);NR_trial_avg(window)],[R_trial_std(window);NR_trial_std(window)]); axis tight;
        set(gca,'FontSize',14);
        title(sprintf('%s - %d neurons in common with %s',session_names{sess_ids(kk)},length(sess_specific_idx),initial_sess),'FontSize',16)
        legend({sprintf('R Average (%d trials)',sum(R_trial_bool)),sprintf('NR Average (%d trials)',sum(NR_trial_bool))})
        legend('show')
        
        
    end
    
    xlabel('Time (seconds)','FontSize',16)
    
    saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_11162018/',...
        sprintf('Relative_to_%s',initial_sess),sprintf('Rat%d_all%sNeurons',rat_id,initial_sess)),'png');

end

%% Persistent vs. unique neurons: activity comparison:
% Start with a given session (Session X):
% -for every session, plot trial-averaged trace (event-locked or not) of: 
% [neurons that were present on Session X] vs. [neurons that were not present on Session X]

RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

% left off with 'Ext2' -- start from 'Ext2' next time. 12.02.2018 21:17 
% left off with 'Ext4' -- start from 'Ext4', ii = 3 (Rat10) next time. 12.03.2018 13:01
% Finito
initial_sess = 'Reinstatement';

%%

for ii = 1:length(RatIDs)

    %% set up params
    rat_id = RatIDs(ii);
    rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    exist_sessions = {};
    all_data_files = dir([rat_folder,filesep,'*.mat']);
    all_data_files = {all_data_files(:).name};
    for dat_i = 1:length(all_data_files)
        num_remove = length(sprintf('Rat%d_',rat_id))+1;
        exist_sessions{dat_i} = all_data_files{dat_i}(num_remove:end-4);
    end
    
    sess_ids = [];
    for kk = 1:length(session_names)
        if any(strcmp(exist_sessions,session_names{kk}))
            sess_ids = [sess_ids,kk];
        end
    end
    
    initial_sess_id = strcmp(initial_sess,session_names(sess_ids));
    neurons2track = logical(cell_map(:,initial_sess_id)>0);

    window = 10:780; % restrict averaging to ~1 second after trial start, and ~2 seconds from trial end
    x_ax = window./10.49;
    
    figure(ii);
    
%     sess2inspect = find(~ismember(initial_sess_id,sess_ids));
    
    %% go through all sessions
    
    for kk = 1:length(sess_ids)
%     for kk = 1:length(sess2inspect)
        
%         sess_id = sess2inspect(kk);
        sess_id = sess_ids(kk);
        
        persist_idx = cell_map(neurons2track,sess_id);
        persist_idx = persist_idx(persist_idx ~= 0);
        
        new_idx = cell_map(~neurons2track,sess_id);
        new_idx = new_idx(new_idx ~= 0);
             
        sess_name = session_names{sess_id};
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        R_trial_bool = logical(sum(full(reshape(Sess_object.event_matrix(:,find(strcmp('cueON',Sess_object.event_names))),799,Sess_object.num_trials)),1));
        NR_trial_bool = ~R_trial_bool;
        
        %% persistent neurons, trial average and CIs
        reshaped = reshape(Sess_object.spikes_conv(persist_idx,1:(799*Sess_object.num_trials)),[length(persist_idx), 799,Sess_object.num_trials]);
        
        R_trial_avg_persist = mean(mean(reshaped(:,:,R_trial_bool),3),1);
        R_trial_sem_persist = bootci(1000,{@mean,mean(reshaped(:,:,R_trial_bool),3)},'type','per','alpha',0.32);
        
        if sum(NR_trial_bool) >= 10
            NR_trial_avg_persist = mean(mean(reshaped(:,:,NR_trial_bool),3),1);
            NR_trial_sem_persist = bootci(1000,{@mean,mean(reshaped(:,:,NR_trial_bool),3)},'type','per','alpha',0.32);
        end
        
        
        if or(strcmp(sess_name,initial_sess),isempty(new_idx))
            
            subplot(length(sess_ids),2,2*(kk-1)+1)
            mseb(x_ax,R_trial_avg_persist(window),R_trial_sem_persist(1,window));
            axis tight; xlabel('Time (seconds)');
            title(sprintf('Response average (%d trials) on %s',sum(R_trial_bool),sess_name))
            legend(sprintf('%d Neurons Total',length(persist_idx)));
            legend('show')
            legend('boxoff')
            
            subplot(length(sess_ids),2,2*kk)
            if sum(NR_trial_bool) >= 10
                mseb(x_ax,NR_trial_avg_persist(window),NR_trial_sem_persist(1,window));
                axis tight; xlabel('Time (seconds)');
                title(sprintf('Non-response average (%d trials) on %s',sum(NR_trial_bool),sess_name))
                legend(sprintf('%d Neurons Total',length(persist_idx)));
                legend('show')
                legend('boxoff')
            end
            
        else
        
            %% session-specific, non-persistent neurons, trial_average and CIs
            reshaped = reshape(Sess_object.spikes_conv(new_idx,1:(799*Sess_object.num_trials)),[length(new_idx), 799,Sess_object.num_trials]);
            
            R_trial_avg_new = mean(mean(reshaped(:,:,R_trial_bool),3),1);
            R_trial_sem_new = bootci(1000,{@mean,mean(reshaped(:,:,R_trial_bool),3)},'type','per','alpha',0.32);
            
            if sum(NR_trial_bool) >= 10
                NR_trial_avg_new = mean(mean(reshaped(:,:,NR_trial_bool),3),1);
                NR_trial_sem_new = bootci(1000,{@mean,mean(reshaped(:,:,NR_trial_bool),3)},'type','per','alpha',0.32);
            end
            
            %         subplot(length(sess2inspect),2,2*(kk-1)+1)
            subplot(length(sess_ids),2,2*(kk-1)+1)
            mseb(x_ax,[R_trial_avg_persist(window);R_trial_avg_new(window)],[R_trial_sem_persist(1,window);R_trial_sem_new(1,window)]);
            axis tight; xlabel('Time (seconds)');
            title(sprintf('Response average (%d trials) on %s',sum(R_trial_bool),sess_name))
            legend({sprintf('%d Shared Neurons',length(persist_idx)),...
                sprintf('%d Unshared neurons',length(new_idx))});
            legend('show')
            legend('boxoff')
            
            %         subplot(length(sess_ids),2,2*kk)
            subplot(length(sess_ids),2,2*kk)
            if sum(NR_trial_bool) >= 10
                mseb(x_ax,[NR_trial_avg_persist(window);NR_trial_avg_new(window)],[NR_trial_sem_persist(1,window);NR_trial_sem_new(1,window)]);
                axis tight; xlabel('Time (seconds)');
                title(sprintf('Non-response average (%d trials) on %s',sum(NR_trial_bool),sess_name))
                legend({sprintf('%d Shared Neurons',length(persist_idx)),...
                    sprintf('%d Unshared neurons',length(new_idx))});
                legend('show')
                legend('boxoff')
            end
        end
            
    end
    %% save figure
    
    saveas(gcf,fullfile('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_12012018/',...
        sprintf('Relative_to_%s',initial_sess),sprintf('Rat%d',rat_id)),'png');

end


%% accumulate average activity around events

event_name = {'leverOUT','leverIN','HLON'};
data_type = 'spikes_conv';
stat_window = [-10 52];

ii = 1;
sess2check = 3:6; % look at extinction sessions

rat_id = RatIDs(ii);
rat_folder = fullfile('rats_processed',sprintf('Rat%d',rat_id));

cell_map = load_cell_map(registration_folder,rat_id);

neurons2track = find(sum(cell_map(:,sess2check) > 0,2)==length(sess2check)); % this line extracts all neurons detected on that day

all_events = [];

for kk = 1:length(sess2check)
    
    sess_specific_idx = cell_map(neurons2track,sess_ids(sess2check(kk)));
    sess_specific_idx = sess_specific_idx(sess_specific_idx ~= 0);
    sess_name = session_names{sess_ids(sess2check(kk))};
    load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
    
%     R_trial_bool = logical(sum(full(reshape(Sess_object.event_matrix(:,find(strcmp('cueON',Sess_object.event_names))),799,Sess_object.num_trials)),1));
%     NR_trial_bool = ~R_trial_bool;
    
    for event_i = 1:length(event_name)
        
        event_times = find(Sess_object.event_matrix(:,strcmp(Sess_object.event_names,event_name{event_i})));
        
        for t = 1:length(event_times)
            all_events = [all_events; [sess_ids(sess2check(kk)),event_i,t,...
                mean(Sess_object.(data_type)(sess_specific_idx,(event_times(t)+stat_window(1)):(event_times(t)+stat_window(2))),2)'] ];
        end
        
    end
        
    
end

%% do some tSNE visualizations


all_data = calcTFIDF(all_events(:,4:end)')';

mappedX = tsne(all_data,[],2,30,25);


figure(1);
% color by event type
scatter(mappedX(all_events(:,2)==1,1),mappedX(all_events(:,2)==1,2),30,'r','filled');
hold on;
scatter(mappedX(all_events(:,2)==2,1),mappedX(all_events(:,2)==2,2),30,'b','filled');
hold on;
scatter(mappedX(all_events(:,2)==3,1),mappedX(all_events(:,2)==3,2),30,'g','filled');

% color by session


figure(2);
Ext1_leverOUT = and(all_events(:,1)==3,all_events(:,2)==1);
scatter(mappedX(Ext1_leverOUT,1),mappedX(Ext1_leverOUT,2),30,'r','filled');

Ext2_leverOUT = and(all_events(:,1)==4,all_events(:,2)==1);
hold on;
scatter(mappedX(Ext2_leverOUT,1),mappedX(Ext2_leverOUT,2),30,'b','filled');

Ext3_leverOUT = and(all_events(:,1)==5,all_events(:,2)==1);
hold on;
scatter(mappedX(Ext3_leverOUT,1),mappedX(Ext3_leverOUT,2),30,'g','filled');

Ext4_leverOUT = and(all_events(:,1)==6,all_events(:,2)==1);
hold on;
scatter(mappedX(Ext4_leverOUT,1),mappedX(Ext4_leverOUT,2),30,'k','filled');


%%

cell_map_binary = logical(cell_map > 0);

overlap_graph = zeros(size(cell_map_binary,2));

for sess_i = 1:7
    
    sess_i_idx = cell_map_binary(:,sess_i);
    
    for sess_j = 1:7
        
        sess_j_idx = cell_map_binary(:,sess_j);
        
        overlap_graph(sess_i,sess_j) = 0;
    end
end

    
    
    


        
    
    
    
    
    
    
    
    






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







    