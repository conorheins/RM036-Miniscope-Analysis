%% Line up calcium traces to behaviorally-relevant events

clear all; clc

%% read in med-pc data

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('alignment_functions'))
addpath(genpath('OASIS_matlab-master'))
addpath(genpath('usefulfunctions'))
addpath(genpath('toolboxes'));

cd([base_directory,filesep,'OASIS_matlab-master'])
setup;
cd(base_directory);

fprintf('Choose behavioral data file (MED-PC output file)\n')
[behav_fnam,behavior_data_folder] = uigetfile('*');

fprintf('Choose neural data file (.mat file)\n')
[results_fnam, results_folder] = uigetfile('*.mat');

fprintf('Choose neuralmapper timestamp file (.txt file)\n')
[tmsp_folder] = uigetdir();

dataArray = medpc_reader(fullfile(behavior_data_folder,behav_fnam));

%%  isolate specific time stamps from different numerical arrays (e.g. K, X, Y, etc.)

tic
[num_trials,trial_starts,cue_starts,cue_ends,HL_on,HL_off,lever_out,lever_in,pellet_tmsp,beambreak_tmsp,...
    ITI_timestamps,TTL_starts,trial_timestamps,press_timestamps] = retrieve_tmsp_array(dataArray);
clear dataArray;
fprintf('Time taken to assemble timestamps: %.2f seconds\n',toc)

%% load neural data, and create struct for the total session, including all trials as sub-structs within it

tic
[ratname, C] = load_session(fullfile(results_folder,results_fnam),'SA1','temporal');
fprintf('Time taken to load data: %.2f seconds\n',toc);

%% smooth data in time (200 ms moving window, 100 ms overlap)

tic

sortedIndices = [reshape(repmat(1:num_trials,799,1),num_trials*799,1), repmat([1:799]',num_trials,1)];
Fs = 10;
window_length = 2;
win_overlap = 1;

C_smooth = smooth_dFF(C,sortedIndices,num_trials,Fs,window_length,win_overlap);

fprintf('Time taken to smooth data: %.2f seconds\n',toc);
%% LATEST PIPELINE FOR DECONVOLUTION AS OF 01.21.2018 

% OVERVIEW:
% 1. Estimate time-constants, baseline fluorescence, and noise STD from
%   entire session's trace data using OASIS
% 2. use learned parameters to deconvolve trace data of individual trials.
%   You do NOT need re-estimate the baseline fluorescence for each trial; 
%   the timescale of fluorescence changes is sometimes so slow that it introduces wild jumps
%   into the baseline estimate if you do it trial-by-trial-- you can use time-constants, noise, &
%   baseline values for entire session to deconvolve individual traces, since
%   baseline across the session is relatively stationary (due to time-lapsed delta-F/F
%   calculation during original fluorescence extraction)

% Time cost: ~1.9 minutes for 118 neurons (79,900 time-steps each), whereas
% foopsi took 7.6 hours with this approach
% Moral of the story: <33333 OASIS 

tic

[C_decon,spikes,noise_vals] = trial_wise_deconvolve(C,sortedIndices,num_trials);

num_neurons = size(C_decon,1);
fprintf('Total time taken to deconvolve %d neurons: %.2f minutes\n',num_neurons,toc/60);

%% create frame-by-event matrix

tic
event_tmsp_array = {lever_out,lever_in,press_timestamps,cue_starts,cue_ends,HL_on,HL_off,beambreak_tmsp,pellet_tmsp};
event_names = {'leverOUT','leverIN','press','cueON','cueOFF','HLON','HLOFF','bbk','pellets'};
event_matrix = create_filter_matrix(sortedIndices,tmsp_folder,num_trials,trial_starts,event_tmsp_array);
fprintf('Time taken to assemble event-filter matrix: %.2f seconds\n',toc)

%% event-detection (i.e. clean up spikes)

tic
fluorescence_array = C_decon;
point_process_array = spikes;
dFF_thr = 0.3;
look_back_f = 10;
look_ahead_f = 30;
spikes_denoised = denoise_spikes(sortedIndices,fluorescence_array,point_process_array,...
    dFF_thr,look_back_f,look_ahead_f);
fprintf('Time taken to denoise spikes: %.2f seconds\n',toc);
    
% figure;
% row = 0.5;
% for neur = 1:num_neurons
%     spks_temp = find(spikes_denoised(neur,:));
%     plot(spks_temp,repmat(row,length(spks_temp),1),'k.','MarkerSize',2);
%     hold on; row = row + 0.5;
% end

kernel = normpdf(-3:3,0,1.5);
spikes_convolved = convolve_timeseries(spikes_denoised,0,kernel);

%% generate array of event-locked calcium traces

tic
data2use = spikes_denoised;

before_align = 50;
after_align = 50;

[event_locked] = align_neural_to_behav(event_names,before_align,after_align,event_matrix,data2use);
fprintf('Time taken to extract neural data aligned to timestamps: %.2f seconds\n',toc);


%% visualizations for Tarun's talk on 09.02.2018

cue_locked = event_locked{4};
cue_locked_mean = mean(cue_locked,3);
cue_locked_mean = bsxfun(@rdivide,cue_locked_mean,max(cue_locked_mean,[],2));

[~,all_maxIDX] = max(cue_locked_mean,[],2);
[~,all_srt] = sort(all_maxIDX,'ascend');

all_sorted = cue_locked_mean(all_srt,:);

imagesc(cue_locked_mean(all_srt,:))

trial_locked = reshape(spikes_convolved,[118 799 100]);
trial_locked_mean = mean(trial_locked,3);
trial_locked_mean = bsxfun(@rdivide,trial_locked_mean,max(trial_locked_mean,[],2));
trial_locked_sorted = trial_locked_mean(all_srt,:);
imagesc(trial_locked_sorted);

full_array = [all_srt,trial_locked_sorted];
full_array = [0:799;full_array];

significant2cue = find(stats_results.SigMatrix(:,4) == 1 & stats_results.Modulation_direction(:,4) == 1);

all_others = true(118,1);
all_others(significant2cue) = false;
all_others = find(all_others);

cue_locked_sig = cue_locked_mean(significant2cue,:);
cue_locked_others = cue_locked_mean(all_others,:);

[~,sig_maxIDX] = max(cue_locked_sig,[],2);
[~,sig_srt] = sort(sig_maxIDX,'ascend');

[~,others_maxIDX] = max(cue_locked_others,[],2);
[~,others_srt] = sort(others_maxIDX,'ascend');

sep_sorted = [cue_locked_sig(sig_srt,:);cue_locked_others(others_srt,:)];

sep_sorted = [ [significant2cue(sig_srt);all_others(others_srt)],sep_sorted];

full_trial_sepsorted = trial_locked_mean(sep_sorted(:,1),:);

sig_order = [significant2cue(sig_srt);all_others(others_srt)];

for trial = 1:100
    temp = cue_locked(sig_order,:,trial);
    imagesc(bsxfun(@rdivide,temp,max(temp,[],2)));
    fprintf('Trial: %d\n',trial);
    pause;
end

good_trials = [4,13,17,37,59];

press_hist = full(mean(reshape(event_matrix(:,3),799,100),2));
bbk_hist = full(mean(reshape(event_matrix(:,8),799,100),2));

press_bytrial = full(reshape(event_matrix(:,3),799,100));
bbk_bytrial = full(reshape(event_matrix(:,8),799,100));

good_trials_press = press_bytrial(:,good_trials);
good_trials_bbk = bbk_bytrial(:,good_trials);

%% visualize trial-by-trial sequencing in select cells
cells2visualize = [17,29,33,37,43,60,68,82,91,106,112,113];

cue_locked = event_locked{4}(cells2visualize',:,:);

cue_locked_conv = cue_locked;
for i = 1:size(cue_locked_conv,1)
    for trial = 1:num_trials
        cue_locked_conv(i,:,trial) = conv(cue_locked(i,:,trial),kernel,'same');
    end
end

trial_avg = mean(cue_locked_conv,3);
trial_avg = bsxfun(@rdivide,trial_avg,max(trial_avg,[],2));

[~,max_inds] = max(trial_avg,[],2);
[~,srt] = sort(max_inds,'ascend');

colorz = parula(ceil(2*length(cells2visualize)));
colorz = colorz(1:length(cells2visualize),:);
% colorz = cool(length(cells2visualize));
plot_ndarray(trial_avg(srt,:),1,1,colorz);

for trial = 1:num_trials
    
    temp = cue_locked_conv(srt,:,trial);
    imagesc(temp);
    
%     [units_i,tmsps_i] = find(temp);
%     units_list = unique(units_i);
%     spikeTimes = cell(length(units_list),1);
%     for j = 1:length(units_list)
%         spikeTimes{j} = (tmsps_i(units_i == units_list(j))./10)';
%     end
%     plotSpikeRaster(spikeTimes,'PlotType','vertline','Autolabel',true,'SpikeDuration',0.1);
    
    pause;
    close gcf;
end

% study response variance across neurons

cells2visualize = [17,29,33,37,43,60,68,82,91,106,112,113];

response_latencies = zeros(size(cue_locked_conv,1),num_trials);
for trial = 1:num_trials
    [~,response_latencies(:,trial)] = max(cue_locked_conv(:,:,trial),[],2);
    response_latencies( sum(cue_locked_conv(:,:,trial),2) == 0 ,trial) = size(cue_locked_conv,2);
end

segments = 1:9:100;

chunked_variance = zeros(size(cue_locked_conv,1),length(segments)-1);
for chunk = 1:length(segments)-1
    chunked_variance(:,chunk) = std(response_latencies(:,segments(chunk): segments(chunk+1)),[],2);
%     chunked_variance(:,chunk) = std(response_latencies(:,segments(chunk): (segments(chunk) + 19)),[],2);
end
    
    
%% Calculate change-scores for every trial and accumulate Change_score_real and change_score_random vectors for each event type

surround_time = [30 30]; % window of frames before and after event (or random timestamp) to use for calculating change score
data2use = spikes_denoised;
[change_scores] = calculate_change_scores(sortedIndices,event_names,event_matrix,data2use,surround_time);

surround_time = [50 50];
event_tmsp = 51;
[stats_results] = permutation_test(event_locked,surround_time,event_tmsp,precision);
% for event_i = 1:length(event_names)
%     for neuron = 1:num_neurons
%         change_real = change_scores{event_i}{1}(neuron,:);
%         change_rand = change_scores{event_i}{2}(neuron,:);
%         
%         [temp1,bincenters1] = hist(change_real,25);
%         [temp2,bincenters2] = hist(change_rand,25);
%         plot(bincenters1,temp1,'b'); hold on;
%         plot(bincenters2,temp2,'r'); hold on;
%         title(sprintf('Change scores for event: %s',event_names{event_i}));
%         pause; hold off;
%     end  
% end

%% test for significance with 2-sided ranksum and/or KS tests, then 1-sided if 2-sided is significant

surround_time = [50 100];
event_tmsp = before_align+1;
[stats_results] = wilcoxon_test_firingrates(event_locked,surround_time,event_tmsp);

[stats_results] = permutation_test(event_locked,surround_time,event_tmsp,precision);

[stats_results] = kstest_change_scores(change_scores,event_names,num_neurons);

significant_modulation = stats_results.SigMatrix;

for event_id = 1:length(event_names)
    
    mkdir(fullfile(base_directory,ratname),event_names{event_id})
    event_directory = fullfile(base_directory,ratname,event_names{event_id});
    
    sig_responding_neurons = find(significant_modulation(:,event_id)); 
    fprintf('%d total neurons are responsive to event: %s\n',length(sig_responding_neurons),event_names{event_id});
    
    for i = 1:length(sig_responding_neurons)
%     for i = 1:length(cells2visualize)
        
        if stats_results.Modulation_direction(sig_responding_neurons(i),event_id) == -1
            mod_direction = 'decreasing';
        elseif stats_results.Modulation_direction(sig_responding_neurons(i),event_id) == 1
            mod_direction = 'increasing';
        else
            mod_direction = 'Conah ya done fucked up';
        end
        fprintf('Displaying neuron: %d, direction of firing rate change: %s\n',sig_responding_neurons(i),mod_direction);
        temp = squeeze(event_locked{event_id}(sig_responding_neurons(i),:,:));

%         fprintf('Displaying neuron: %d\n',cells2visualize(i));
%         temp = squeeze(event_locked{event_id}(cells2visualize(i),:,:));
        
        
        [tmsps_i,trials_i] = find(temp);
        trials_list = unique(trials_i);
        spikeTimes = cell(length(trials_list),1);
        for j = 1:length(trials_list)
            spikeTimes{j} = (tmsps_i(trials_i == trials_list(j))./10)';
        end
        plotSpikeRaster(spikeTimes,'PlotType','vertline','Autolabel',true,'SpikeDuration',0.1);
        figname = fullfile(event_directory,['Neuron',num2str(sig_responding_neurons(i)),'_',event_names{event_id},'.png']);
        saveas(gcf,figname);
        
        

%         temp = temp(sum(temp,1) > 1,:);
%         temp = logical(temp');
%         plotSpikeRaster(temp,'PlotType','vertline');
%         
%         plot_ndarray(temp,0,1,repmat([0 0 0.25],size(temp,2),1));

%         save_flag = input('Would you like to save this plot for display purposes? (y/n)\n','s');
%         if strcmp(save_flag,'y')
%             figname = ['Neuron',num2str(sig_responding_neurons(i)),'_',event_names{event_id},'.png'];
%             saveas(gcf,figname);
%         end

        close gcf;
    end
    
end

% find some of the nice-looking cue-locked neural responses and pull out
% some trial-averages, other shit, etc.
neurons_to_note = [17,29,33,43,106];

event_id = 4;
temp = event_locked{event_id}(neurons_to_note,:,:);

for i = 1:size(temp,1)
    for trial = 1:num_trials
        temp(i,:,trial) = conv(temp(i,:,trial),kernel,'same');
    end
end

trial_mean = mean(temp,3);
[~,max_inds] = max(trial_mean,[],2);
[~,new_order] = sort(max_inds,'ascend');    

trial_mean = trial_mean(new_order,:);
neurons_to_note = neurons_to_note(new_order);
%%        
    
preCue_data = zeros(num_neurons,num_trials); postCue_data = zeros(num_neurons,num_trials);

for trial = 1:num_trials
    trial_struct = session_struct.(['Trial_',num2str(trial)]);
    frame_times = trial_struct.frame_times;
    
    event_times = trial_struct.bbk_tmsp(1);
    [~,align_ind] = min(abs(frame_times - event_times));
    preCue_data(:,trial) = mean(trial_struct.Activity_smooth(:,(align_ind-30):(align_ind)),2);
    postCue_data(:,trial) = mean(trial_struct.Activity_smooth(:,(align_ind+1):(align_ind+30)),2);
end

for neuron = 1:num_neurons
    [precue,bincenters1] = hist(preCue_data(neuron,:),25);
    [postcue,bincenters2] = hist(postCue_data(neuron,:),25);
    
    plot(bincenters1,precue);
    hold on; plot(bincenters2,postcue);
    pause; hold off;
    
end


%% assemble population vectors aligned to every behavioral event (use mean activity/firing rate)

num_neurons = size(C,1);

population_vectors = [];

for trial = 1:num_trials
    
    trial_struct =  session_struct.(['Trial_',num2str(trial)]);
    frame_times = trial_struct.frame_times;
    
    event_types = fieldnames(trial_struct);
    event_types = event_types(~cellfun(@isempty,strfind(event_types,'tmsp')));
    
    spk_tms_temp = trial_struct.Activity_spikes_kept;
    dFF_temp = trial_struct.Activity_smooth;

    for event = 1:length(event_types)
        
        event_nam = event_types{event};
        
        if ~isempty(trial_struct.(event_nam))
            
            within_trial_events = trial_struct.(event_nam);
            pop_vect_event_temp = zeros(num_neurons,length(within_trial_events));
            
            for event_i = 1:length(within_trial_events)
                pop_vect_temp = zeros(num_neurons,1);
                event_times = within_trial_events(event_i);
                [~,align_ind] = min(abs(frame_times - event_times));
                for neuron = 1:num_neurons
                    spks = find(spk_tms_temp(neuron,(align_ind: (align_ind + 20))));
                    if ~isempty(spks)
                        for spk_i = 1:length(spks)
                            spk = spks(spk_i);
                            pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1) + mean(dFF_temp(neuron,(spk:spk+5)));
                        end
                        pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1)./length(spks);
                    else
                        pop_vect_temp(neuron,1) = 0;
                    end
                end
                pop_vect_event_temp(:,event_i) = pop_vect_temp;
            end
            
            population_vectors = [population_vectors, [repmat(event,1,length(within_trial_events)); pop_vect_event_temp]];
                
                
%             event_time = trial_struct.(event_nam)(1);
%             
%             [~,align_ind] = min(abs(frame_times - event_time)); % this line finds the frame time-stamp that is closest in time to the event time-stamp
%             
%             pop_vect_temp = zeros(num_neurons,1);
%             
%             for neuron = 1:num_neurons
%                 spks = find(spk_tms_temp(neuron,(align_ind: (align_ind + 50))));
%                 if ~isempty(spks)
%                     for spk_i = 1:length(spks)
%                         spk = spks(spk_i);
%                         pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1) + mean(dFF_temp(neuron,(spk:spk+20)));
%                     end
%                     pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1)./length(spks);
%                 else
%                     pop_vect_temp(neuron,1) = 0;
%                 end
%             end
%             
%             population_vectors = [population_vectors, [event; pop_vect_temp]];
               
        end
    end
end


vectors_and_labels = population_vectors'; 
labels = vectors_and_labels(:,1); vects = vectors_and_labels(:,2:end);

pattern_dists = 1-squareform(pdist(vects,'cosine'));
[sorted_labels,sort_inds] = sort(labels);
imagesc(pattern_dists(sort_inds,sort_inds));

sorted_vects = vects(sort_inds,:);
[~,scores] = pca(sorted_vects);

colors = hsv(1000);
colors = colors(randperm(1000,length(event_types)),:);
for lab = 1:length(event_types)
    num_labs = length(find(sorted_labels==lab));
    scatter3(scores(sorted_labels==lab,1),scores(sorted_labels==lab,2),scores(sorted_labels==lab,3),10,repmat(colors(lab,:),num_labs,1),'DisplayName',event_types{lab});
    hold on;
end
legend('show')

%tsne it 
tsne_map = tsne(sorted_vects,sorted_labels,2,5,50);

colors = hsv(500);
colors = colors(randperm(500,length(event_types)),:);
for lab = 1:length(event_types)
    num_labs = length(find(sorted_labels==lab));
    scatter(tsne_map(sorted_labels==lab,1),tsne_map(sorted_labels==lab,2),10,repmat(colors(lab,:),num_labs,1),'DisplayName',event_types{lab});
    hold on;
end
legend('show')

%% assemble population vectors aligned to every behavioral event (represent in terms of spike latencies)

tic

Fs = 10; %sampling frequency of camera
num_neurons = size(C,1);

population_vectors = [];

kernel = normpdf(-3:3,0,1.5);

time_after_sec = 2;
time_after_f = floor(time_after_sec * Fs);

for trial = 1:num_trials
    
    trial_struct =  session_struct.(['Trial_',num2str(trial)]);
    frame_times = trial_struct.frame_times;
    
    event_types = fieldnames(trial_struct);
    event_types = event_types(~cellfun(@isempty,strfind(event_types,'tmsp')));
    
    spk_tms_temp = trial_struct.Activity_spikes_kept;
    dFF_temp = trial_struct.Activity_smooth;
    spks_conv_temp = zeros(size(spk_tms_temp));
    
    for neuron = 1:num_neurons
        spks_conv_temp(neuron,:) = conv(spk_tms_temp(neuron,:),kernel,'same');
    end
    
    for event = 1:length(event_types)
        
        event_nam = event_types{event};
        
        if ~isempty(trial_struct.(event_nam))
            
            within_trial_events = trial_struct.(event_nam);
            
            
            % below for loop collects all time-stamps at intervals of every
            % (time_after_sec) to not collect overlapping spikes 
            
            if length(within_trial_events) > 1
                
                curr_tmsp = within_trial_events(1);
                tmsp = 2;
                
                throw_out_idx = false(length(within_trial_events),1);
                
                for tmsp = tmsp:length(within_trial_events)     
                    if (within_trial_events(tmsp) - curr_tmsp) <= time_after_sec
                        throw_out_idx(tmsp) = true;
                    elseif (within_trial_events(tmsp) - curr_tmsp) > time_after_sec
                        curr_tmsp = within_trial_events(tmsp);
                        tmsp = tmsp + 1;
                    end 
                end
                
                within_trial_events = within_trial_events(~throw_out_idx);
                
            end
                            
            pop_vect_event_temp = zeros(num_neurons,length(within_trial_events));
            
            
            % below determines the encoding of population events -- either
            % use some measure of the calcium concentration (i.e. dF/F0 averaged
            % after each spike) or the latency code (time of max activity
            % -- found by max in convolved spike train. Comment/uncomment
            % according block of code
            
            for event_i = 1:length(within_trial_events)
                pop_vect_temp = zeros(num_neurons,1);
                event_times = within_trial_events(event_i);
                [~,align_ind] = min(abs(frame_times - event_times));
%                 for neuron = 1:num_neurons
%                     spks = find(spk_tms_temp(neuron,(align_ind: (align_ind + 20))));
%                     if ~isempty(spks)
%                         for spk_i = 1:length(spks)
%                             spk = spks(spk_i);
%                             pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1) + mean(dFF_temp(neuron,(spk:spk+5)));
%                         end
%                         pop_vect_temp(neuron,1) = pop_vect_temp(neuron,1)./length(spks);
%                     else
%                         pop_vect_temp(neuron,1) = 0;
%                     end
%                 end
                [~,pop_vect_event_temp(:,event_i)] = max(spks_conv_temp(:,align_ind : (align_ind + time_after_f)),[],2);
            end
            
            population_vectors = [population_vectors, [repmat(event,1,length(within_trial_events)); pop_vect_event_temp]];
        end
    end
end

fprintf('Time taken to assemble population vectors: %.2f seconds\n',toc);

%% do geometry on them

vectors_and_labels = population_vectors'; 
labels = vectors_and_labels(:,1); vects = vectors_and_labels(:,2:end);

vects = vects(:,sum(vects > 1) > prctile(sum(vects > 1),25));
active_idx = sum(vects > 1,2) > prctile( sum(vects > 1,2),25);

patternSimilarity = 1-squareform(pdist(vects,'cosine'));

[sorted_labels,sort_idx] = sort(labels,'ascend');
imagesc(patternSimilarity(sort_idx,sort_idx));

sorted_vects = vects(sort_idx,:);
[~,scores] = pca(sorted_vects);

num_events_firing = sum(vects > 1); %when index (i,j) of vects is 1, indicates that neuron j did not fire in the 2 seconds following event i

rarely_active = num_events_firing <= prctile(num_events_firing,10);

active_vects = vects(:,~rarely_active);
sorted_vects = active_vects(sort_idx,:);
[~,scores] = pca(sorted_vects);


tsne_map = tsne(sorted_vects,sorted_labels,2,5,50);

colors = hsv(500);
colors = colors(randperm(500,length(event_types)),:);
for lab = 1:length(event_types)
    num_labs = length(find(sorted_labels==lab));
    scatter(tsne_map(sorted_labels==lab,1),tsne_map(sorted_labels==lab,2),10,repmat(colors(lab,:),num_labs,1),'DisplayName',event_types{lab});
    hold on;
end
legend('show')     

for lab = 1:length(event_types)
    
    if lab == 4
        scatter(tsne_map(sorted_labels == lab,1),tsne_map(sorted_labels==lab,2),30,'b');
    else
        scatter(tsne_map(sorted_labels == lab,1),tsne_map(sorted_labels == lab,2),30,'r');
    end
    
    hold on;
end

%% create struct with event-locked data (one sub-struct for each of the event types)

num_neurons = size(C,1);

event_types = {'leverOUT_tmsp','leverIN_tmsp','press_tmsp','cueON_tmsp','cueOFF_tmsp','HLON_tmsp',...
    'HLOFF_tmsp','bbk_tmsp'};

event_locked_data = struct;
for neuron = 1:num_neurons
    event_locked_data.(['Neuron',num2str(neuron)]) = struct;
    for event = 1:length(event_types)
        event_locked_data.(['Neuron',num2str(neuron)]).(event_types{event}) = [];
    end
end
    
for trial = 1:num_trials
    
    trial_struct =  session_struct.(['Trial_',num2str(trial)]);
    frame_times = trial_struct.frame_times;
    
    event_types = fieldnames(trial_struct);
    event_types = event_types(~cellfun(@isempty,strfind(event_types,'tmsp')));
    
    temp_data_all = trial_struct.Activity_smooth;

    for event = 1:length(event_types)
        
        event_nam = event_types{event};
        
        if ~isempty(trial_struct.(event_nam))
            event_times = trial_struct.(event_nam)(1);
            
            [~,align_ind] = min(abs(frame_times - event_times)); % this line finds the frame time-stamp that is closest in time to the event time-stamp
            
            for neuron = 1:num_neurons
                temp_event_dat = [temp_data_all(neuron, (align_ind - 49) : align_ind)', temp_data_all(neuron, (align_ind+1):(align_ind+50))'];
                event_locked_data.(['Neuron',num2str(neuron)]).(event_nam) = [event_locked_data.(['Neuron',num2str(neuron)]).(event_nam); temp_event_dat];
            end
        end
    end
   
end

%% assemble tensor of population activity after timestamp of interest ( Neuron x time_after_f x num_trials )

Fs = 10;
event_nam = 'cueON_tmsp';
time_after = 20; %number of seconds before event to look 
time_after_f = time_after * Fs; % convert seconds to frames by multiplying by frame rate

tensor = []; % don't initialize dimensions of tensor because we don't know how many trials event of interest happens on

for trial = 1:num_trials
    
    trial_struct =  session_struct.(['Trial_',num2str(trial)]);
    frame_times = trial_struct.frame_times;
    
    temp_data_raw = trial_struct.Activity;
    temp_data_smooth = trial_struct.Activity_smooth;
    temp_data_decon = trial_struct.Activity_decon;
    temp_data_spikes = trial_struct.Activity_spikes;
    temp_data_spikes_kept = trial_struct.Activity_spikes_kept;
     
    data2use = temp_data_spikes_kept;

    if ~isempty(trial_struct.(event_nam))
        event_times = trial_struct.(event_nam)(1);
        
        if length(event_times) > 1
            for j = 1:length(event_times)
                [~,align_ind] = min(abs(frame_times - event_times(j))); % this line finds the frame time-stamp that is closest in time to the event time-stamp
                if align_ind > time_before_f
                    break
                end
            end
        else
            [~,align_ind] = min(abs(frame_times - event_times)); % this line finds the frame time-stamp that is closest in time to the event time-stamp
        end
        
        tensor = cat(3,tensor,data2use(:, align_ind : (align_ind + time_after_f - 1) ));
    end
    
end

%% visualize and do stuff with the trial-by-trial tensor of event-locked data
%%%%%%%%%%%%%%%%%%%%
% center of mass calculation, sort neurons by their center of mass
%%%%%%%%%%%%%%%%%%%%

tensor_convolved = tensor;
kernel = normpdf(-3:3,0,1.5);

for trial = 1:num_trials
    for neuron = 1:num_neurons
        temp = squeeze(tensor(neuron,:,trial));
        tensor_convolved(neuron,:,trial) = conv(temp,kernel,'same');
    end
end

tensor_convolved(isnan(tensor_convolved))=0;
mean_sequence = mean(tensor_convolved,3);

% com_idx=repmat(1:size(mean_sequence,2),num_neurons,1);
% mass=sum(mean_sequence,2);
% com=sum((mean_sequence.*com_idx),2)./mass;
% com = round(com);
% 
% [~,sort_inds] = sort(com,'ascend');

[~,max_inds] = max(mean_sequence,[],2);
[~,sort_inds] = sort(max_inds,'ascend');

figure(1);
imagesc(mean_sequence(sort_inds,:));
xlabel('Seconds')
ylabel('Neuron ID')

figure(2);
for i = 1:num_trials
    imagesc(squeeze(tensor_convolved(sort_inds,:,i)));
    xlabel('Seconds since CUE ON');
    ylabel('Neuron Number'); 
    pause;
end

% sorted_mean_seq = mean_sequence(sort_inds,:);
% 
% highly_active = find(max(sorted_mean_seq,[],2) >= 0.5);
% sorted_mean_seq = sorted_mean_seq(highly_active,:);
% 
% sorted_tensor = tensor(sort_inds,:,:);
% sorted_tensor = sorted_tensor(highly_active,:,:);
% 
% figure(1);
% imagesc(sorted_mean_seq);
% xlabel('Seconds since CUE ON');
% ylabel('Neuron Number');
% 
% figure(2);
% for i = 1:num_trials
%     imagesc(squeeze(sorted_tensor(:,:,i)));
%     xlabel('Seconds since CUE ON');
%     ylabel('Neuron Number'); 
%     pause;
% end

%% look at neurons that spike at predictable times in the post-event epoch

active_or_not = false(num_neurons,num_trials);
vectorized = zeros(num_neurons,num_trials);
for trial = 1:num_trials
    active_or_not(:,trial) = sum(tensor_convolved(:,:,trial),2) > 0;
    [~,vectorized(:,trial)] = max(tensor_convolved(:,:,trial),[],2);
end

% neurons that were active at least 50% of the time
active_neurons = find(sum(active_or_not,2) >= num_trials/2);
active_neurons_firing_tms = vectorized(active_neurons,:);
active_neurons_timing_std = std(active_neurons_firing_tms,0,2);

predictable_units = active_neurons((active_neurons_timing_std <= 50));

reduced_tensor_conv = tensor_convolved(predictable_units,:,:);
reduced_tensor = tensor(predictable_units,:,:);
mean_seq = mean(reduced_tensor,3);
[~,sort_inds] = max(mean_seq,[],2);
[~,sort_inds] = sort(sort_inds,'ascend');
imagesc(mean_seq(sort_inds,:));

for i = 1:num_trials
    imagesc(squeeze(reduced_tensor(sort_inds,:,i)));
    pause;
end

%% vectorize each population response in terms of spike-timing, then compare distances between vectors


vectorized = zeros(num_neurons,num_trials);

for trial = 1:num_trials
    [~,vectorized(:,trial)] = max(tensor_convolved(:,:,trial),[],2);
end

active_neur = find(vectorized(median(vectorized,2) > 1));

vectorized = vectorized(active_neur,:);
active_neur_tensor = tensor(active_neur,:,:);
active_neur_tensor_conv = tensor_convolved(active_neur,:,:);

timing_var = std(vectorized,0,2);
tight_firing = find(timing_var <= prctile(timing_var,25));

for i = 1:num_trials
    imagesc(squeeze(active_neur_tensor(tight_firing,:,i)));
    xlabel('Seconds since CUE ON');
    ylabel('Neuron Number'); 
    pause;
end

mean_sequence = mean(active_neur_tensor_conv,3);
[~,max_inds] = max(mean_sequence,[],2);
[~,sort_inds] = sort(max_inds,'ascend');

figure(1);
imagesc(mean_sequence(sort_inds,:));
xlabel('Seconds')
ylabel('Neuron ID')

figure(2);
for i = 1:num_trials
    imagesc(squeeze(active_neur_tensor_conv(sort_inds,:,i)));
    xlabel('Seconds since CUE ON');
    ylabel('Neuron Number'); 
    pause;
end


spike_time_dists = squareform(pdist(vectorized,'euclidean'));
[V,D] = eig(spike_time_dists); % decompose distance of spike-timing between neurons to find patterns of correlated spike patterns (e.g. sequences)
projection = vectorized'*V(:,1:2); % project data onto first eigenvector of this decomposition

pattern_dists = squareform(pdist(vectorized','cosine')); %look at cosine distance between multineuronal patterns

scatter(projection,sum(pattern_dists));
    
%% permutation testing of single neuron responses

tic
selectivities = zeros(num_neurons,length(event_types));
for neuron = 1:num_neurons
    temp_neuron_dat = event_locked_data.(['Neuron',num2str(neuron)]);
    for event = 1:length(event_types)
        event_nam = event_types{event};
        temp_event_dat = temp_neuron_dat.(event_nam);
%         test_stat = mean(temp_event_dat(:,2)) - mean(temp_event_dat(:,1));
        test_stat = length(find(temp_event_dat(:,2) > 0))/length(find(temp_event_dat(:,1) > 0));
        shuffled = temp_event_dat(:);
        null_dist = zeros(10000,1);
        for null = 1:length(null_dist)
            shuffled = shuffled(randperm(length(shuffled)));
%             mean1 = mean(shuffled(1:ceil(length(shuffled)/2))); 
%             mean2 = mean(shuffled((length(shuffled)/2+1):end));
%             null_dist(null) = mean1-mean2;
            count1 = length(find(shuffled(1:ceil(length(shuffled)/2)) > 0));
            count2 = length(find(shuffled((length(shuffled)/2+1):end) > 0));
            null_dist(null) = count1/count2;
        end
        if test_stat < mean(null_dist)
            selectivities(neuron,event) = length(find(null_dist < test_stat))/10000;
        elseif test_stat > mean(null_dist)
            selectivities(neuron,event) = length(find(null_dist > test_stat))/10000;
        end
    end
end
fprintf('Time taken to do permutation testing: %.2f minutes',(toc/60) );           
            
%% trial-wise functional connectivity

num_neurons = size(C,1);

seg_length = 12; % number of frames over which to calculate functional correlations

seg_edges = 1:12:792; % segments over which to calculate FC
num_segs = length(seg_edges)-1;

FC_tensor = zeros(num_neurons,num_neurons,num_trials,num_segs);
for trial = 1:num_trials
    chunk = session_struct.(['Trial_',num2str(trial)]).Activity_smooth(:,4:end-4)';
    for seg = 1:num_segs
        FC_tensor(:,:,trial,seg) = corrcoef(chunk(seg_edges(seg):seg_edges(seg+1),:));
    end    
end

meanFC = squeeze(mean(FC_tensor,3));

% video of average FC over time
color_range = [min(meanFC(:)) max(meanFC(:))];
for seg = 1:num_segs
    imagesc(meanFC(:,:,seg));
    caxis(color_range);
    pause(1.2);
end

corrMat = corrcoef(C_smooth');
corrMat = tril(corrMat,-1);
correlations = corrMat(corrMat ~= 0);

high_correlations = corrMat > 0.4;
[rowID,colID] = find(high_correlations); 
pairs = [rowID,colID]; clear rowID colID

for pair = 1:length(pairs)
    fprintf('Neuron %d and Neuron %d\n',pairs(pair,1),pairs(pair,2))
    for seg = 1:num_segs
        hist(squeeze(FC_tensor(pairs(pair,1),pairs(pair,2),:,seg)),25);
        pause;
    end
end
    


%% choose neurons and cycle through traces
   
neuron = 2;

for trial = 1:num_trials
    
    temp_trial_struct = session_struct.(['Trial_',num2str(trial)]);
    activity_decon = temp_trial_struct.Activity_decon(neuron,:);
    activity = temp_trial_struct.Activity(neuron,:);
    frame_tmsp = temp_trial_struct.frame_tmsp;
    leverOUT_tmsp = temp_trial_struct.leverOUT_tmsp;
    leverIN_tmsp = temp_trial_struct.leverIN_tmsp;
    press_tmsp = temp_trial_struct.press_tmsp;
    cueON_tmsp = temp_trial_struct.cueON_tmsp; 
    cueOFF_tmsp = temp_trial_struct.cueOFF_tmsp; 
    HLON_tmsp = temp_trial_struct.HLON_tmsp; 
    HLOFF_tmsp = temp_trial_struct.HLOFF_tmsp; 
    bbk_tmsp = temp_trial_struct.bbk_tmsp; 
    
    figure; 
    title(['Trial ' num2str(trial)])
%     F0 = median(activity_decon(1:round(HLON_tmsp-frame_tmsp(1))));
%     activity_decon = (activity_decon - F0)./F0;
    activity_decon = activity_decon./max(activity_decon);
    activity = activity./max(activity);
    plot(frame_tmsp,activity_decon,'DisplayName','Calcium trace'); ylim([0 1.1]);axis tight
    hold on; plot(frame_tmsp,activity,'Color',[0 0 1],'DisplayName','Raw fluorescence'); 
    hold on; plot([leverOUT_tmsp leverOUT_tmsp], [0 0.5],'r-','LineWidth',2.5,'DisplayName','LeverIn')
    hold on; plot([leverIN_tmsp leverIN_tmsp], [0 0.5],'r-','LineWidth',2.5,'DisplayName','LeverOut')
    
    if ~isempty(press_tmsp)
        hold on; stem(press_tmsp,repmat(0.5,1,length(press_tmsp)),'k','LineWidth',1.5,'DisplayName','LeverPress')
        hold on; plot([cueON_tmsp cueON_tmsp], [0 0.5],'LineWidth',2.5,'Color',[0 0.5 0.5],'DisplayName','Cue ON')
        hold on; plot([cueOFF_tmsp cueOFF_tmsp], [0 0.5],'LineWidth',2.5,'Color',[0 0.5 0.5],'DisplayName','Cue OFF')
    end
    
    hold on; plot([HLON_tmsp HLON_tmsp], [0 0.5],'LineWidth',2.5,'Color',[1 0.5 0],'DisplayName','HL ON')
    hold on; plot([HLOFF_tmsp HLOFF_tmsp], [0 0.5],'LineWidth',2.5,'Color',[1 0.5 0],'DisplayName','HL OFF')
    
    if ~isempty(bbk_tmsp)
        hold on; stem(bbk_tmsp,repmat(0.5,1,length(bbk_tmsp)),'LineWidth',1.5,'Color',[0 0.5 1],'DisplayName','Beambreaks')
    end

    legend('show')
    
    pause;
    close gcf
end

%% align data to specific time-stamps

event_aligned = [];
temporal_data = C'; 
temporal_data = bsxfun(@rdivide,bsxfun(@minus,temporal_data,min(temporal_data)),(max(temporal_data)-min(temporal_data)));

spatial_data = all_sessions{2,2}.A;
spatial_data = full(bsxfun(@rdivide,spatial_data,max(spatial_data,[],1)));

trial_average = zeros(size(spatial_data,1),301);

for trial = 1:num_trials
    
    temp_event = session_struct.(['Trial_',num2str(trial)]).event_tmsp;
    temp_frames = session_struct.(['Trial_',num2str(trial)]).frame_tmsp;
    
    if ~isempty(temp_event)
        align_time = temp_event(1);
        
        [~,align_ind] = min(abs(temp_frames - align_time));
        
        trial_temporal = temporal_data(sortedIndices(:,1)==trial,:);
        for cell = 1:size(trial_temporal,2)
            trial_temporal(:,cell) = locsmooth(trial_temporal(:,cell),10,2,1);
        end
        
        trial_data = spatial_data*trial_temporal';
        
        corrected = trial_data(:,(align_ind-100):(align_ind+200));
        
        trial_average = trial_average + corrected;
        
    end
        
end

trial_average = trial_average./100;


%% create temporally-colored heatmap of activations relative to cue-onset

video = reshape(trial_average,400,400,size(trial_average,2));

filt_rad=16; % gauss filter radius
filt_alpha=6; % gauss filter alpha
lims=4; % contrast prctile limits (i.e. clipping limits lims 1-lims)
cmap= colormap(jet);%  cubehelix(200,[0.9,-1,7,1]));
per=3; % baseline percentile (0 for min)
bgcolor=[ .75 .75 .75 ]; % rgb values for axis background
time_select=0;
startT = 1;
stopT = size(video,3);

[rows,columns,frames]=size(video);

disp('Gaussian filtering the movie data...');

h=fspecial('gaussian',filt_rad,filt_alpha);
video=imfilter(video,h,'circular');

disp(['Converting to df/f using the ' num2str(per) ' percentile for the baseline...']);

baseline=repmat(prctile(video,per,3),[1 1 frames]);

dff=((video.^2-baseline.^2)./(baseline)).*100;

% take the center of mass across dim 3 (time) for each point in space
disp('Computing the center of mass...');

com_idx=zeros(1,1,frames);

for i=1:frames
	com_idx(:,:,i)=i;
end

com_idx=repmat(com_idx,[rows columns 1]);

mass=sum(dff,3);
com_dff=sum((dff.*com_idx),3)./mass;
com_dff(isnan(com_dff))=0;

max_proj=std(dff,[],3);
max_proj(isnan(max_proj))=0;

disp('Creating images...');

clims(1)=prctile(max_proj(:),lims);
clims(2)=prctile(max_proj(:),100-lims);

norm_max_proj=min(max_proj,clims(2)); % clip to max
norm_max_proj=max(norm_max_proj-clims(1),0); % clip min
norm_max_proj=norm_max_proj./(clims(2)-clims(1)); % normalize to [0,1]

% map to [0,1] for ind2rgb

% Relative scaling between [0,1]: 
clims(1)=min(com_dff(:));
clims(2)=max(com_dff(:));

norm_dff=min(com_dff,clims(2)); % clip to max
norm_dff=max(norm_dff-clims(1),0); % clip min
norm_dff=norm_dff./(clims(2)-clims(1)); % normalize to [0,1]

idx_img=round(norm_dff.*size(cmap,1));
im1_rgb=ind2rgb(idx_img,cmap);

%
 cbar_idxs=linspace(0,size(cmap,1),1e3);
  imwrite(im1_rgb,'Filename.png','Alpha',norm_max_proj);
 I = imread('Filename.png', 'BackgroundColor',[0 0 0]);
imshow(I);
 imwrite(I, 'Rat21_Reinstatement_cuelocked.jpg')

%%

for trial = 1:num_trials
    
    temp_event = session_struct.(['Trial_',num2str(trial)]).event_tmsp;
    temp_frames = session_struct.(['Trial_',num2str(trial)]).frame_tmsp;
    
    if ~isempty(temp_event)
        align_time = temp_event(1);
        
        [~,align_ind] = min(abs(temp_frames - align_time));
        
        trial_data = temporal_data(sortedIndices(:,1)==trial,:);
        
%         preEvent = trial_data((align_ind-50):(align_ind-1),:);
%         baselines = prctile(preEvent,5,1);
        
%         mu = median(preEvent,1); %baseline correction
%         sigma = mad(preEvent,1,1);
%         sigma = std(preEvent,[],1);
%         sigma(sigma==0)=1;
  
%         corrected = bsxfun(@rdivide,bsxfun(@minus,trial_data((align_ind-20):(align_ind+30),:),mu),sigma);
        
%         corrected = bsxfun(@minus,trial_data((align_ind-50):(align_ind+100),:),baselines);
        corrected = trial_data((align_ind-50):(align_ind+100),:);
        
        event_aligned = cat(3,event_aligned,corrected);
        
    end
        
end


% visualize

% event_vect = zeros(151,1);
% event_vect([51 71]) = 1;

for i = 1:size(event_aligned,2)
    
    plot(squeeze(event_aligned(:,i,:)),'b-','LineWidth',1);
%     hold on; stem(event_vect); axis tight

%     for trial = 1:100
%         temp_frames = session_struct.(['Trial_',num2str(trial)]).frame_tmsp;
%         temp_event = session_struct.(['Trial_',num2str(i)]).event_tmsp;
%         [~,align_ind] = min(abs(temp_frames - align_time));

    hold on; hist([51 71],151); axis tight;
    disp(num2str(i));
    pause; hold off;
end

%%


responsive_trial_IDs = find(any(neuron27 > 0.2,1));
responsive_trials = neuron27(:,responsive_trial_IDs);

plot(responsive_trials,'r-','LineWidth',0.5);
hold on; plot(mean(responsive_trials,2),'b-','LineWidth',2.5)

responsive_trial_IDs = find(any(neuron72 > 0.4,1));
responsive_trials = neuron72(:,responsive_trial_IDs);

plot(responsive_trials,'r-','LineWidth',0.5);
hold on; plot(mean(responsive_trials,2),'b-','LineWidth',2.5)


responsive_trial_IDs = find(any(neuron85 > 0.25,1));
responsive_trials = neuron85(:,responsive_trial_IDs);

plot(responsive_trials,'r-','LineWidth',0.5);
hold on; plot(mean(responsive_trials,2),'b-','LineWidth',2.5);


smoothed27 = zeros(size(neuron27)); 
for trial = 1:100
    smoothed27(:,trial) = locsmooth(neuron27(:,trial),10.5,2,1);
end
    
    


