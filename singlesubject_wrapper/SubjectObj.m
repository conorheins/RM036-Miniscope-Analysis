classdef SubjectObj < handle
    
    % This class is a wrapper for loading and analyzing multineuronal Ca2+ activity for a single subject (rat) from 
    % a single session
    % Author: Conor Heins, conor.heins@gmail.com
    
    %% properties
    properties
        ratname;         % ratname  (string)
        session_name;    % session name, i.e. 'SA1', 'SA2', 'Ext1', etc.. (string)
        num_trials;      % number of trials
        behav_fnam;      % behavioral data file (.mpc file)
        neural_fnam;     % neural data file (.mat file)
        tmsp_fdir;       % NeuralMapper timestamp directory (folder full of .txt files)
        C;               % raw fluorescence data (N x T, where N is number of neurons)
        sortedIndices;   % T x 2 array of trial number (first column) and trial-relative frame index (second column)
        Fs;              % imaging rate (in fps). If not given, estimated from inter-timestamp-interval distribution from NeuralMapper files
        C_smooth;        % smoothed fluorescence (N x T)
        C_decon;         % deconvolved fluorescence (N x T)
        spikes;          % estimated spike trains (N x T)
        noise_vals;      % noise-values estimated from OASIS deconvolutino
        spikes_conv;     % Gaussian-convolved spike trains
        spikes_denoised; % denoised spikes
        event_names;     % cell of strings, each representing the name of an event type (1 x numEventTypes)
        event_matrix;    % array of frame-timestamps and events (T x numEventTypes)
        event_locked;    % cell of event-locked neural data (1 x numEventTypes cell array, where each entry is a N x windowSize x numEventInstances tensor)
        stats_results;   % structure comprised of 4 matrices: significance (0s and 1s), directionality of modulation (-1, +1), associated pVal (0 < pVal < 1), and response fidelity (proportion of trials with significant firing, only in significant responding neurons)
        meta_data;       % structure containing meta-data about analysis choices
        bad_trials;      % trials that were removed (due to artifacts, etc.)
        bad_neurons;     % neurons that should be excluded (e.g. due to bad deconvolution results, etc.)
    end
    
    %% methods
    methods
        %% load data
        function load_data(obj,rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder)
            
            if ~exist('rat_id','var') || isempty(rat_id)
                error('Must provide a subject ID');
            end
                  
            if ~exist('session_name','var') || isempty(session_name)
                warning('No session type selected...defaulting to loading Self-Admin Session 1 (SA1)')
                session_name = 'SA1';
                obj.session_name = session_name;
            else
                obj.session_name = session_name;
            end
            
            if ~exist('behav_data_folder','var') || isempty(behav_data_folder)
                
                if isempty(obj.behav_fnam)
                    fprintf('Choose behavioral data file (MED-PC output file)\n')
                    [behav_fnam_within,path] = uigetfile('*');
                    obj.behav_fnam = fullfile(path,behav_fnam_within);
                end
                
            else
                                
                if rat_id < 10
                    temp = ['R0',num2str(rat_id)];
                else
                    temp = ['R',num2str(rat_id)];
                end
                
                path = fullfile(behav_data_folder,temp);
                datfile_nm = dir(path);
                datfile_nm = {datfile_nm.name};
                
                sess_fileIDX = find(cellfun(@(x) ~isempty(x),strfind(datfile_nm,session_name)));
                
                if length(sess_fileIDX) == 1
                    behav_fnam_within = datfile_nm{sess_fileIDX};
                    obj.behav_fnam = fullfile(path,behav_fnam_within);
                elseif length(sess_fileIDX) > 1
                    behav_fnam_within = datfile_nm(sess_fileIDX);
                    obj.behav_fnam = fullfile(path,behav_fnam_within);
                end
                
            end
            
            if ~exist('neural_data_folder','var') || isempty(neural_data_folder)
                
                if isempty(obj.neural_fnam)
                    fprintf('Choose neural data file (.mat file)\n')
                    [neural_fnam_within, neural_folder] = uigetfile('*.mat');
                    obj.neural_fnam = fullfile(neural_folder,neural_fnam_within);
                end
                
            else
                
                datfile_nm = dir([neural_data_folder,filesep,'*.mat']);
                datfile_nm = {datfile_nm.name};
                neural_fnam_within = datfile_nm{cellfun(@(x) ~isempty(x),strfind(datfile_nm,num2str(rat_id)))};
                obj.neural_fnam = fullfile(neural_data_folder,neural_fnam_within);
                
            end
            
            if ~exist('timestamps_folder','var') || isempty(timestamps_folder)
                
                if isempty(obj.tmsp_fdir)
                    fprintf('Choose neuralmapper timestamp folder (directory of .txt files)\n')
                    [obj.tmsp_fdir] = uigetdir();
                end
                
            else
                
                obj.tmsp_fdir = fullfile(timestamps_folder,['Rat',num2str(rat_id),'_Timestamps',filesep,session_name]);
                
            end 
                
            [obj.ratname, obj.C, obj.num_trials, obj.sortedIndices, obj.event_names, obj.event_matrix, obj.Fs] = load_data(obj.behav_fnam,obj.session_name,obj.neural_fnam,obj.tmsp_fdir);
            
        end 
        
        function bad_trials = quickView(obj)
            
            [NaN_neuron,~] = find(isnan(obj.C));
            NaN_neuron = unique(NaN_neuron);
            toshow = obj.C;
            toshow(NaN_neuron,:) = [];

            for i = 1:length(NaN_neuron)
                fprintf('This neuron has NaN entries: %d\n',NaN_neuron(i))
            end
            
            plot(sum(toshow,1));
            limz = ylim;
            
            [~,trial_starts] = ismember(unique(obj.sortedIndices(:,1)),obj.sortedIndices(:,1),'rows');
            
            x_coords = reshape([trial_starts,trial_starts,nan(obj.num_trials,1)]',3*obj.num_trials,1);
            y_coords = reshape([limz(1)*ones(obj.num_trials,1),limz(2)*ones(obj.num_trials,1),nan(obj.num_trials,1)]',3*obj.num_trials,1);
            hold on;
            plot(x_coords,y_coords);
            xticks(trial_starts); xticklabels({unique(obj.sortedIndices(:,1))});
            
            bad_trials = [];      
            keep_choosing_flag = 1;
            
            start_filtering = input('Do you want to choose some trials to throw out?(y/n)\n','s');
            if strcmp(start_filtering,'y')
                while keep_choosing_flag
                    temp = input('Enter indices of bad trials: \n');
                    bad_trials = [bad_trials,temp];
                    keep_going = input('Keep selecting?(y/n)\n','s');
                    
                    if strcmp(keep_going,'y')
                        keep_choosing_flag = 1;
                    else
                        keep_choosing_flag = 0;
                    end                   
                end
            end 
            
            close gcf;
            
        end  
                 
        function remove_trials(obj,trial_idx)
            timestamps = find(ismember(obj.sortedIndices(:,1),trial_idx));
            obj.C(:,timestamps) = [];
            obj.sortedIndices(timestamps,:) = [];
            obj.event_matrix(timestamps,:) = [];
            obj.bad_trials = trial_idx;
            obj.num_trials = obj.num_trials - length(trial_idx);
        end
        
        function filter_neurons(obj,keep_or_except_flag,neuron_idx)
            
            if strcmp(keep_or_except_flag,'keep') || isempty(keep_or_except_flag)
                obj.C = obj.C(neuron_idx,:);
            elseif strcmp(keep_or_except_flag,'remove')
                obj.C(neuron_idx,:) = [];
            end
            
        end
             
        function smooth_C(obj,window_length,win_overlap)
            
            if isempty(obj.Fs)
                obj.Fs = input('Please input frame rate: \n');
            end
            
            obj.C_smooth = smooth_dFF(obj.C,obj.sortedIndices,obj.num_trials,obj.Fs,window_length,win_overlap);
            
        end
        
        function deconvTrials(obj)
            
            [obj.C_decon,obj.spikes,obj.noise_vals] = trial_wise_deconvolve(obj.C,obj.sortedIndices,obj.num_trials);
            
        end
        
        function find_bad_neurons(obj)
            
            bad_neurons = [];
            keep_choosing_flag = 1;
            
            figure;
            imagesc(obj.spikes);
            
            start_filtering = input('Do you want to choose some neurons to throw out?(y/n)\n','s');
            if strcmp(start_filtering,'y')
                while keep_choosing_flag
                
                    temp = input('Enter rows of bad neurons: \n');
                    bad_neurons = [bad_neurons,temp];
                    keep_going = input('Keep selecting?(y/n)\n','s');
                    
                    if strcmp(keep_going,'y')
                        keep_choosing_flag = 1;
                    else
                        keep_choosing_flag = 0;
                    end
                    
                end
            end
            
            obj.bad_neurons = bad_neurons;
            
            close gcf;
            
        end

        
        function denoiseSpikes(obj,cont_data,point_process_data,denoising_params)
            
            if ~exist('cont_data','var') || isempty(cont_data)
                cont_data = obj.C;
            end
            
            if ~exist('point_process_data','var') || isempty(point_process_data)
                if isempty(obj.spikes)
                    warning('No existing spike array! Deconvolving trials now to obtain spikes\n')
                    [obj.C_decon,obj.spikes,~] = trial_wise_deconvolve(obj.C,obj.sortedIndices,obj.num_trials);
                    point_process_data = obj.spikes;
                end
            end
            
            dFF_thr = denoising_params.dFF_thr;
            look_back_f = denoising_params.look_back_f;
            look_ahead_f = denoising_params.look_ahead_f;
            counts_flag = denoising_params.counts_flag;
            
            obj.spikes_denoised = denoise_spikes(obj.sortedIndices,cont_data,point_process_data,...
                dFF_thr,look_back_f,look_ahead_f,counts_flag,obj.bad_neurons);
            
        end
        
        function overlay_spikes_calc(obj,cont_data,point_process_data,spacing_interval)
            
            if ~exist('cont_data','var') || isempty(cont_data)
                cont_data = obj.C;
            end
            
            if ~exist('point_process_data','var') || isempty(point_process_data)
                point_process_data = obj.spikes_denoised;
            end
            
            max_normed = bsxfun(@rdivide,bsxfun(@minus,cont_data,min(cont_data,[],2)),(max(cont_data,[],2) - min(cont_data,[],2)));
            max_normed(isnan(max_normed)) = 0;
            
            if ~exist('spacing_interval','var') || isempty(spacing_interval)
                spacing_interval = max(range(max_normed,2));
            end
            
            num_neurons = size(max_normed,1);
            T = size(max_normed,2);
            
            point_process_data = double(point_process_data > 0);
            cont_shifted = max_normed + repmat([spacing_interval:spacing_interval:(spacing_interval*num_neurons)]',1,T);
            plot(cont_shifted');
                        
            hold on;
            
            for i = 1:num_neurons
                
                spike_times = find(point_process_data(i,:));
                
                tick_bottom = prctile(cont_shifted(i,:),5);
                tick_top = tick_bottom + 0.5 * range(cont_shifted(i,:));
                
                xpoints = [spike_times;spike_times;nan(1,length(spike_times))];
                ypoints = repmat([tick_bottom;tick_top;NaN],1,length(spike_times));
                
                xpoints = xpoints(:);
                ypoints = ypoints(:);
                
                plot(xpoints,ypoints,'k');
                                                
            end
            
            axis tight;
            
        end
            
        
        function convSpikes(obj,data_array,transpose_flag,kernel)
            
            if ~exist('data_array','var') || isempty(data_array)
                data_array = obj.spikes;
            end
            
            if ~exist('tranpose_flag','var') || isempty(transpose_flag)
                transpose_flag = 0; % assume data doesn't need to be tranposed
            end
            
            if ~exist('kernel','var') || isempty(kernel)
                kernel = normpdf(-3:3,0,1.5);
            end

            obj.spikes_conv = convolve_timeseries(data_array,transpose_flag,kernel);
            
        end
        
        
        function lock2events(obj,data_type,before_align,after_align)
                
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            [obj.event_locked,obj.meta_data.event_tmsp] = align_neural_to_behav(obj.event_names,before_align,after_align,obj.event_matrix,data2align);
            obj.meta_data.data_type_align = data_type;
            
        end
        
        function compute_selectivity(obj,test_type,surround_time,num_shuffles,alpha)
            
            
            if strcmpi(test_type,'permutation')
                
                if ~exist('num_shuffles','var') || isempty(num_shuffles)
                    num_shuffles = 100;
                end
                
                if ~exist('alpha','var') || isempty(alpha)
                    alpha = 0.01;
                end
                
                obj.stats_results = permutation_test(obj.event_locked,surround_time,obj.meta_data.event_tmsp,num_shuffles,alpha);
                
            elseif or(strcmpi(test_type,'wilcoxon'),strcmpi(test_type,'ks'))
                
                obj.stats_results = test_firingrates(test_type,obj.event_locked,surround_time,obj.meta_data.event_tmsp);
                
            end
            
            obj.meta_data.surround_time_4stats = surround_time;
            obj.meta_data.num_shuffles = num_shuffles;
            obj.meta_data.alpha = alpha;
            
            
        end
        
        function display_aligned_responses(obj,event_of_interest,surround_time,modulation_direction)
            
            if ~exist('modulation_direction','var') || isempty(modulation_direction)
                modulation_direction = 'both';
            end
            
            event_idx = cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest));
            traces = obj.event_locked{event_idx};
            event_tmsp = obj.meta_data.event_tmsp;
            
            if strcmp(modulation_direction,'pos')
                sig_neuron_idx = obj.stats_results.Modulation_direction(:,event_idx) == 1;
            elseif strcmp(modulation_direction,'neg')
                sig_neuron_idx = obj.stats_results.Modulation_direction(:,event_idx) == -1;
            elseif strcmp(modulation_direction,'both')
                sig_neuron_idx = obj.stats_results.SigMatrix(:,event_idx) == 1;
            end
            
            % regardless of which modulation type is shown, always use non-significant responders as comparison 
            ns_neuron_idx = obj.stats_results.SigMatrix(:,event_idx) == 0;
            
            display_idx = (event_tmsp - surround_time(1)):(event_tmsp + surround_time(2));
            
            numSig = length(find(sig_neuron_idx));
            sig_responders = traces(sig_neuron_idx,display_idx,:);
            trial_average_sig = mean(sig_responders,3);
            trial_average_sig = bsxfun(@rdivide,trial_average_sig,max(trial_average_sig,[],2));
            trial_average_sig(isnan(trial_average_sig)) = 0;
            ensemble_mean_sig = mean(trial_average_sig,1);

            numNS = length(find(ns_neuron_idx));
            nonsig = traces(ns_neuron_idx,display_idx,:);
            trial_average_ns = mean(nonsig,3);
            trial_average_ns = bsxfun(@rdivide,trial_average_ns,max(trial_average_ns,[],2));
            trial_average_ns(isnan(trial_average_ns)) = 0;
            ensemble_mean_ns = mean(trial_average_ns,1);
            
            
            display_min = min([ensemble_mean_sig(:);ensemble_mean_ns(:)]);
            display_max = max([ensemble_mean_sig(:);ensemble_mean_ns(:)]);
            
            figure('Position',[200 400 1100 300])
            subplot(121)
            legend_title_sig = sprintf('Average of %d neurons aligned to %s',numSig,event_of_interest);
            plot(ensemble_mean_sig,'r-','DisplayName',legend_title_sig)
            axis tight; ylim([display_min,display_max]);
            legend('show')
            subplot(122)
            legend_title_ns = sprintf('Average of %d neurons aligned to %s',numNS,event_of_interest);
            plot(ensemble_mean_ns,'b-','DisplayName',legend_title_ns)
            axis tight; ylim([display_min,display_max]);
            legend('show')
           
        end
        
        function display_multiunit(obj,event_of_interest,surround_time,subset_idx,trial_idx,avg_flag)
            
            event_tmsp = obj.meta_data.event_tmsp;
            event_idx = cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest));
            traces = obj.event_locked{event_idx};
            display_idx = (event_tmsp - surround_time(1)):(event_tmsp + surround_time(2));
            
            % if no subset of neuron IDs provided, defaults to positively-modulated ones
            if ~exist('subset_idx','var') || isempty(subset_idx) 
                subset_idx = obj.stats_results.Modulation_direction(:,event_idx) == 1;
            end
            
            % if no subset of trial IDs provided, defaults to all trials
            if ~exist('trial_idx','var') || isempty(trial_idx)
                trial_idx = 1:obj.num_trials;
            end
            
            % if no avg_flag provided, assumes to not average and show
            % individual trials
            if ~exist('avg_flag','var') || isempty(avg_flag)
                avg_flag = false;
            end
            
            dat2display = traces(subset_idx,display_idx,trial_idx);
            colors = winter(size(dat2display,1));
            if avg_flag
                dat2display = mean(dat2display,3);
                dat2display = bsxfun(@rdivide,dat2display,max(dat2display,[],2));
                baseline = dat2display(:, 1: (surround_time(1)) );
                baseline_corrected = bsxfun(@minus,dat2display,min(baseline,[],2));
                for ii = 1:size(baseline_corrected,1)
                    plot(baseline_corrected(ii,:),'color',colors(ii,:));
                    hold on;
                end
                plot([surround_time(1), surround_time(1)],[min(baseline_corrected(:)),max(baseline_corrected(:))],'r--','LineWidth',1.5)
                axis tight; xlabel('Time ( 100 ms frames )'); ylabel('Relative fluorescence change (max normalized)');
                title(['Trial-Averaged dF/F for units modulated by event: ',event_of_interest])
                fprintf('Currently displaying average over trials: %s\n',num2str(trial_idx))
            else
                max_norm_factor = max(max(dat2display,[],3),[],2);
%                 dir_nm = fullfile('figures_02282018_2',sprintf('%s_Triallocked',event_of_interest));
%                 mkdir(dir_nm);
                for trial = 1:length(trial_idx)
                    trial_dat = squeeze(dat2display(:,:,trial));
                    trial_dat = bsxfun(@rdivide,trial_dat,max_norm_factor);
                    baseline = trial_dat(:,1: (surround_time(1)));
                    baseline_corrected = bsxfun(@minus,trial_dat,min(baseline,[],2));
                    for ii = 1:size(baseline_corrected,1)
                        plot(baseline_corrected(ii,:),'color',colors(ii,:));
                        hold on;
                    end
                    plot([surround_time(1), surround_time(1)],[min(baseline_corrected(:)),max(baseline_corrected(:))],'r--','LineWidth',1.5)
                    axis tight; xlabel('Time ( 100 ms frames )'); ylabel('Relative fluorescence change (max normalized)');
%                     saveas(gcf,fullfile(dir_nm,['Trial',num2str(trial),'.png']));
                    fprintf('Currently displaying activity on Trial: %d\n',trial_idx(trial))
%                     close gcf;
                    pause; close gcf;
                end
            end
            
        end
        
        function displayRasters(obj,event_of_interest,subset_idx,trial_idx)
            
            if ~ismember(obj.meta_data.data_type_align,{'spikes','spikes_denoised'})
                error('Must use spike raster data -- use display_multiunit() method for continuous event-locked data')
            end
            
            event_idx = cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest));
            spikes_binar = obj.event_locked{event_idx};
            event_tmsp = obj.meta_data.event_tmsp;
            
            if strcmp(obj.meta_data.data_type_align,'spikes')
                spikes_binar = double(spikes_binar > 0);
            end
            
            % if no subset of neuron IDs provided, defaults to positively-modulated ones
            if ~exist('subset_idx','var') || isempty(subset_idx) 
                subset_idx = obj.stats_results.Modulation_direction(:,event_idx) == 1;
            end
            
            % if no subset of trial IDs provided, defaults to all trials
            if ~exist('trial_idx','var') || isempty(trial_idx)
                trial_idx = 1:obj.num_trials;
            end
            
%             dir_nm = fullfile('figures_02282018_2',sprintf('%s_Triallocked_spikes',event_of_interest));
%             mkdir(dir_nm);
            
            dat2display = spikes_binar(subset_idx,:,:);
             for trial = 1:length(trial_idx)
                 trial_dat = squeeze(dat2display(:,:,trial));
                 spk_array = cell(size(trial_dat,1),1);
                 for neur = 1:length(spk_array)
                     spk_array{neur} = find(trial_dat(neur,:))./obj.Fs;
                 end
                 plotSpikeRaster(spk_array,'PlotType','vertline','Autolabel',true);
                 max_y = ylim; max_y = max_y(2);
                 hold on; h2 = plot( [event_tmsp/obj.Fs,event_tmsp/obj.Fs],[0,max_y],'r--','LineWidth',1.5);
                 title(sprintf('Neurons aligned to event: %s',event_of_interest))
                 ylabel('Unit ID');
%                  saveas(gcf,fullfile(dir_nm,['Trial',num2str(trial),'.png']));
%                  close gcf;
                 legend(h2,sprintf('Time of event: %s',event_of_interest))
                 pause; close gcf
             end
                 
        end
                       
        function stats_summary(obj)
            
            num_neurons = size(obj.C,1);
            
            sig_modulated_all = sum(obj.stats_results.SigMatrix == 1, 1);
            positively_modulated = sum(obj.stats_results.Modulation_direction == 1,1);
            negatively_modulated = sum(obj.stats_results.Modulation_direction == -1,1);
            
%             colors = winter(numel(sig_modulated_all));

            ymin_val = 0;
            ymax_val = max([sig_modulated_all,positively_modulated,negatively_modulated]);
            figure;
            subplot(311)
%             for bar_i = 1:numel(sig_modulated_all)
%                 bar(bar_i, (sig_modulated_all(bar_i)./num_neurons) * 100,'FaceColor',colors(bar_i,:));
%                 hold on;
%             end
%             hold off;
            dat2plot = (sig_modulated_all./num_neurons) * 100;
            bar(diag(dat2plot),'stacked');
            xticklabels(obj.event_names); ylim([ymin_val,ymax_val]);
            ylabel('Percent selective')
            title('Total significantly-modulated neurons')

            subplot(312)
%             for bar_i = 1:numel(positively_modulated)
%                 bar(bar_i, (positively_modulated(bar_i)./num_neurons) * 100,'FaceColor',colors(bar_i,:));
%                 hold on;
%             end
%             hold off;
            dat2plot = (positively_modulated./num_neurons) * 100;
            bar(diag(dat2plot),'stacked');
            xticklabels(obj.event_names); ylim([ymin_val,ymax_val]);
            ylabel('Percent selective')
            title('Positively-modulated neurons')

            subplot(313)
%             for bar_i = 1:numel(negatively_modulated)
%                 bar(bar_i, (negatively_modulated(bar_i)./num_neurons) * 100,'FaceColor',colors(bar_i,:));
%                 hold on;
%             end
%             hold off;     
            dat2plot = (negatively_modulated./num_neurons) * 100;
            bar(diag(dat2plot),'stacked')
            xticklabels(obj.event_names); ylim([ymin_val,ymax_val]);
            ylabel('Percent selective')
            title('Negatively-modulated neurons')

        end
        
        function compute_fidelity(obj)
            
            if isempty(obj.stats_results)
                error('Must compute single-cell statistics before performing operations like e.g. fidelity calculations')
            end
            
            event_tmsp = obj.meta_data.event_tmsp;
            surround_time = obj.meta_data.surround_time_4stats;
            
            obj.stats_results.fidelity_matrix = compute_fidelity(obj.event_locked,event_tmsp,surround_time,obj.stats_results);
                      
        end
        
        function histogram_data = fidelity_distribution(obj,event_of_interest,bin_width)
            
            if isempty(obj.stats_results)
                error('Must compute single-cell statistics before performing operations like e.g. fidelity calculations')
            end
            
            if ~isfield(obj.stats_results,'fidelity_matrix')
               obj.compute_fidelity();
            end
            
            if ~exist('event_of_interest','var') || isempty(event_of_interest)
                warning('Defaulting to computing distributions for all events, in absence of specified event_type')
                event_idx = 1:length(obj.event_names);
            else
                event_idx = find(cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest)));
            end
            
            if ~exist('bin_width','var') || isempty(bin_width)
                bin_width = 0.05;
            end
            
            fidelity_matrix = obj.stats_results.fidelity_matrix;
            
            if length(event_idx) > 1
                histogram_data = cell(length(event_idx),1);
            else
                histogram_data = [];
            end
            
            figure;
            colors = parula(length(event_idx));
            
            for ii = 1:length(event_idx)
                
                event_id = event_idx(ii);
                fidelity_vals = fidelity_matrix(~isnan(fidelity_matrix(:,event_id)),event_id);
                [probs,bin_edges] = histcounts(fidelity_vals,'BinWidth',bin_width);
                bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
                plot(bin_centers,probs,'LineWidth',1.25,'Color',colors(ii,:),'DisplayName',obj.event_names{event_id})
                hold on;
                
                if length(event_idx) > 1
                    histogram_data{ii} = [probs',bin_centers'];
                else
                    histogram_data = [probs',bin_centers'];
                end
                                
            end
            xlabel('Fidelity (proportion of all trials)')
            ylabel('Number of units')
            legend('show');
            title('Response fidelity distribution')
                    
        end
        
        function cdf_estimates = fidelity_dist_cdf(obj,histogram_data,event_of_interest)
            
            if ~exist('histogram_data','var') || isempty(histogram_data)
            
                histogram_data = obj.fidelity_distribution(event_of_interest);
                
            end
            
            figure;
            if iscell(histogram_data)
                colors = parula(length(histogram_data));
                
                for event_id = 1:length(histogram_data)
                    cd = histogram_data{event_id}(:,1);
                    cd = cumsum( cd / sum(cd) );
                    plot(histogram_data{event_id}(:,2),cd,'Color',colors(event_id,:),'DisplayName',obj.event_names{event_id})
                    hold on;
                end
                
                legend('show');
                title('Cumulative distribution of response fidelity for different behavioral events')
                xlabel('Percentage of trials');
                ylabel('Cumulative probability'); 
                hold off;
            else
                cd = histogram_data(:,1);
                cd = cumsum(cd/sum(cd));
                plot(histogram_data(:,2),cd,'DisplayName',event_of_interest)
                title(sprintf('Cumulative response fidelity for event: %s',event_of_interest))
                xlabel('Percentage of trials');
                ylabel('Cumulative probability');
            end
            
            
        end
            
        function sig_correlations = compute_signal_correlation(obj,event_of_interest,after_align)
            
            if ~exist('event_of_interest','var') || isempty(event_of_interest)
                warning('Defaulting to computing distributions for all events, in absence of specified event_type')
                event_idx = 1:length(obj.event_names);
            else
                event_idx = find(cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest)));
            end
            
            if ~exist('after_align','var') || isempty(after_align)
                after_align = obj.meta_data.surround_time_4stats(2);
            end
            
            event_tmsp = obj.meta_data.event_tmsp;
            
            if length(event_idx) > 1
                sig_correlations = cell(length(event_idx),1);
            else
                sig_correlations = [];
            end
           
            for ii = 1:length(event_idx)
                
                event_id = event_idx(ii);
                
                condition_traces = obj.event_locked{event_id};

                range = event_tmsp: min((event_tmsp+after_align),size(condition_traces,2));
                
                signal_avg = mean(condition_traces(:,range,:),3);
                
                if length(event_idx) > 1
                    correlations = corr(signal_avg','type','Spearman'); 
                    correlations(isnan(correlations)) = 0;
                    sig_correlations{ii} = correlations;
                else
                    correlations = corr(signal_avg','type','Spearman'); 
                    correlations(isnan(correlations)) = 0;
                    sig_correlations = correlations;
                end
                
            end
                                           
        end
        
        function noise_correlations = compute_noise_correlations(obj)
            
            sig_modulated = obj.stats_results.SigMatrix;
            [num_neurons,num_events] = size(sig_modulated);
            event_matrix_trialwise = permute(reshape(full(obj.event_matrix)',num_events,[],obj.num_trials),[2,1,3]);
            trial_length = size(event_matrix_trialwise,1);

            tensor_mask = false(num_neurons,trial_length,obj.num_trials);

            for neuron = 1:num_neurons
    
                sig_events = logical(sig_modulated(neuron,:));
                trial_inds2include = true(trial_length,obj.num_trials);
    
                for trial = 1:obj.num_trials
        
                    event_times = sum(event_matrix_trialwise(:,sig_events,trial),2);
                    starts = find(event_times);
        
                    for ii = 1:length(starts)
            
                        idx = starts(ii);
                        trial_inds2include( (idx:min(trial_length, idx + 100)), trial) = false;
            
                    end
        
                end
    
                tensor_mask(neuron,trial_inds2include) = true;
    
            end

            reshaped_mask = reshape(tensor_mask,num_neurons,trial_length*obj.num_trials);
            
            data2use = obj.(obj.meta_data.data_type_align);
            
            noise_correlations = zeros(num_neurons,num_neurons);

            for neuron = 1:num_neurons
                
                mask = reshaped_mask(neuron,:);
                noise_correlations(neuron,:) = corr(data2use(neuron,mask)',data2use(:,mask)','type','Spearman');
                
            end
        end   
    
        function sorted_spikes = sort_calcium_events_obj(obj,spk_array,fluorescence_array,options)
            
            if ~exist('spk_array','var') || isempty(spk_array)
                spk_array = obj.spikes;
            end
            
            if ~exist('fluorescence_array','var') || isempty(fluorescence_array)
                fluorescence_array = obj.C;
            end
            
            if ~exist('options','var') || isempty(options)
                options.ISI_thresh = 10;
                options.before_win = 20;
                options.after_win = 30;
                options.K = 2;
            end
            
            sorted_spikes = sort_calcium_events(spk_array,fluorescence_array,options);
            
        end
        
        function population_PSTH(obj,data_type,surround_time,binsize)
            
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            [event_locked,event_tmsp] = align_neural_to_behav(obj.event_names,surround_time(1),surround_time(2),obj.event_matrix,data2align);
            
            nrows = ceil(length(obj.event_names)/3);
            ncols = ceil(length(obj.event_names)/nrows);
            
            event_tmsp = event_tmsp/obj.Fs;
           
            figure('Position',[400 400 700 600])
            
            for event_i = 1:length(obj.event_names)
                
                data2show = event_locked{event_i};
                
                if and(~isempty(strfind(lower(data_type),'spikes')),and(exist('binsize','var'),~isempty(binsize)))
                    
                    max_t = ceil(size(data2show,2)./obj.Fs);
                    
                    binned = zeros(size(data2show,1),length(0:binsize:max_t),size(data2show,3));
                    
                    for trial = 1:size(data2show,3)
                        for neur = 1:size(data2show,1)
                            dat = data2show(neur,:,trial);
                            temp = find(dat)./obj.Fs; % turns binary spike pattern into spike times in seconds 
                            binned(neur,:,trial) = binspikes(temp,(1./binsize),[0 max_t]); % chronux function for binning
                        end
                    end
                    
                    pop_sum = squeeze(sum(binned,1));
                    pop_average = mean(pop_sum,2);
                    pop_average(isnan(pop_average)) = 0;
                    
                    t = 0:binsize:max_t;
                    
        
                else
                    
                    pop_average = mean(sum(data2show,1),3);
                    pop_average(isnan(pop_average)) = 0;
                    
                    t = [1:size(data2show,2)]./obj.Fs;
                    
                end
                
                subplot(nrows,ncols,event_i)
                plot(t,pop_average);
                ylimz = ylim;
                hold on; plot([event_tmsp, event_tmsp],[ylimz(1),ylimz(2)],'r--','LineWidth',1.5)
                title(sprintf('Population activity aligned to: %s',obj.event_names{event_i}))
                xlabel('Time (s)')
                ylabel('Smoothed firing rate')
                hold off; axis tight;

                
            end
                    
        end
        
        function reshaped_distances = sequence_vectors(obj,event_of_interest,data_type,surround_time)
            
            event_idx = find(cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest)));
            
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            [event_locked,event_tmsp] = align_neural_to_behav(obj.event_names,surround_time(1),surround_time(2),obj.event_matrix,data2align);
            
            data = event_locked{event_idx};
            
            [num_neurons,sequence_length,num_trials] = size(data);
            
            trialwise_mean = mean(data,3);
            
            [~,max_inds] = max(trialwise_mean,[],2);
            
            [sorted_maxtimes,srt_scheme] = sort(max_inds,'ascend');
            
            spktime_patterns = zeros(num_trials,num_neurons);
            
            for patt = 1:num_trials
                
                temp = data(:,:,patt);
                
                [maxValz,spktime_patterns(patt,:)] = max(temp,[],2);
                spktime_patterns(patt,maxValz == 0) = NaN;
                
            end
            
            pattern_distance = zeros(num_neurons,num_neurons,num_trials);
            
            
            %% below spike-time pattern distance is inspired by Raichmann & Ben-Jacob, 2008:
            % 'Identifying repeating motifs in the activation of synchronized bursts in cultured neuronal networks'
            
            for patt = 1:num_trials 
                
                currPatt = spktime_patterns(patt,:);
                
                for ii = 1:num_neurons
                    for jj = 1:num_neurons
                        if isnan(currPatt(ii) + currPatt(jj))
                            pattern_distance(ii,jj,patt) = 1; % if either neuron did not fire, put distance between spike-times is maximal
                        else
                            pattern_distance(ii,jj,patt) = (currPatt(ii) - currPatt(jj))/sequence_length; % otherwise, distance measure is difference between spike-times, normalized by total sequence length
                        end
                    end
                end
            end
            
            reshaped_distances = zeros( num_trials, (num_neurons^2 - num_neurons)/2 );
            
            for patt = 1:num_trials
                reshaped_distances(patt,:) = get_lower_tri(pattern_distance(:,:,patt));
            end
                 
        end
        
        function pattern_vectors = pattern_vectors(obj,event_of_interest,data_type,surround_time)
            
            event_idx = find(cellfun(@(x) ~isempty(x), strfind(obj.event_names,event_of_interest)));
            
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            [event_locked,event_tmsp] = align_neural_to_behav(obj.event_names,surround_time(1),surround_time(2),obj.event_matrix,data2align);
            
            data = event_locked{event_idx};
            
            pattern_vectors = squeeze(mean(data,2))';
        end
        
        function [trial_average,trial_average_CIs,pop_average,pop_average_CIs] = population_activity_wholeTrial(obj,data_type,bound_type,CI_bounds,neur_idx,trial_idx,plot_flag,throw_out_indices)
            
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            if ~exist('plot_flag','var')
                plot_flag = true;
            end
            
            if exist('throw_out_indices','var')
                data2align(:,throw_out_indices) = [];
            end
            
            if isempty(neur_idx)
                neur_idx = 1:size(data2align,1);
            end
            
            data2align = data2align(neur_idx,:);
            
            num_neurons = size(data2align,1);
            num_events = length(obj.event_names);
            trial_length = size(data2align,2)/obj.num_trials;
            reshaped = reshape(data2align,[num_neurons, trial_length, obj.num_trials]);
            reshaped_events = permute(reshape(full(obj.event_matrix),[],obj.num_trials,num_events),[1,3,2]);

            if ~exist('trial_idx','var') || isempty(trial_idx) % in case that trials to plot are not specific, use all trials
                num_trials = obj.num_trials;
            else
                reshaped = reshaped(:,:,trial_idx);
                reshaped_events = reshaped_events(:,:,trial_idx);
                num_trials = size(reshaped,3);
            end
            
            trial_average = mean(reshaped,3); % all neurons' trial-averaged activity
            pop_average = mean(trial_average,1); % within time-bin averages of trial-averaged neural activity
            trial_average_CIs = zeros(num_neurons,trial_length,2);
            pop_average_CIs = zeros(2,trial_length);
            
            switch lower(bound_type)
                case 'bootstrap'
                    
                    num_bootstraps = 100;
                    
                    lower_bound = CI_bounds(1)/100;
                    upper_bound = CI_bounds(2)/100;

                    for neur_i = 1:size(reshaped,1)
                        
                        tmp = datasample(reshaped(neur_i,:,:),num_bootstraps*num_trials,3);
                        bootstrap_dists_neur_i = reshape(squeeze(tmp),trial_length,num_bootstraps,[]);
                        sample_means = mean(bootstrap_dists_neur_i,3);
                        
                        for timebin = 1:size(sample_means,1)
                            
                            [cd,bins] = histcounts(sample_means(timebin,:),50);
                            cd = cumsum(cd/sum(cd));
                            
                            trial_average_CIs(neur_i,timebin,1) = bins(find(cd>upper_bound,1));
                            trial_average_CIs(neur_i,timebin,2) = bins(find(cd>lower_bound,1));
                            
                        end
                        
                    end
                    
                    
                    % bootstrapping to obtain deviation of cross-neuron firing rate
                    % average from the mean firing rate
                    
                    num_bootstraps = 1000;
                    bootstrap_dists = datasample(trial_average,num_neurons*num_bootstraps,1);
                    bootstrap_dists_r = reshape(bootstrap_dists,num_neurons,num_bootstraps,[]);
                    sample_means = squeeze(mean(bootstrap_dists_r,1));
                    
                    for timebin = 1:size(sample_means,2)
                        
                        [cd,bins] = histcounts(sample_means(:,timebin),50);
                        cd = cumsum(cd/sum(cd));
                        
                        pop_average_CIs(1,timebin) = bins(find(cd>upper_bound,1));
                        pop_average_CIs(2,timebin) = bins(find(cd>lower_bound,1));
                        
                    end
                case 'jackknife'
                    
                    
                    [pop_average,SDest] = jackknife(trial_average');
                    
                    pop_average_CIs(1,:) = SDest;
                    CIpop_average_CIss(2,:) = SDest;
                    
                    pop_average = pop_average'; % because of how Chronux (source of jackknife function) shapes their arrays
                    
                case 'analytic'
                    SEM = std(trial_average,0,1)./sqrt(num_neurons); % equation for SEM
                    pop_average_CIs(1,:) = SEM;
                    pop_average_CIs(2,:) = SEM;
            end
            
            if plot_flag
                event_tmspz = zeros(1,num_events);
                for event_id = 1:num_events
                    [~,event_tmspz(event_id)] = max(sum(reshaped_events(:,event_id,:),3));
                end
                
                [sorted_tmspz,srt] = sort(event_tmspz,'ascend');
                
                figure('Position',[100 400 1200 400])
                xt = (1:length(pop_average)) * (1./obj.Fs);
                plot(xt,pop_average,'b-','LineWidth',2,'DisplayName','Average firing rate')
                hold on; plot(xt,pop_average_CIs(1,:),'b-','LineWidth',0.25,'DisplayName',sprintf('%.1f CI',CI_bounds(2)));
                hold on; plot(xt,pop_average_CIs(2,:),'b-','LineWidth',0.25,'DisplayName',sprintf('%.1f CI',CI_bounds(1)));
                ylimz = ylim;
                
                event_names = obj.event_names;
                event_names{find(cellfun(@(x) ~isempty(x),strfind(event_names,'bbk')))} = 'Foodport Entry';
                
                colors = parula(length(event_names));
                
                for event_id = 1:length(event_names)
                    plot([sorted_tmspz(event_id),sorted_tmspz(event_id)]./obj.Fs,[ylimz(1),ylimz(2)],'--','LineWidth',1.5,'Color',colors(srt(event_id),:),'DisplayName',event_names{srt(event_id)})
                end
                
                xlim([1 75])
                legend('show','Location','Northeast')
                
                title('Average single-neuron activity with various experimental events time-stamped')
                xlabel('Time (seconds)')
                ylabel('Trial-averaged firing rate')
                
                hold off;
            end
        end
        
        function [trial_average,CIs] = PCA_trialplot(obj,data_type,numPCs,trial_indices,bound_type,CI_bounds,plot_flag)
            
            if ~exist('data_type','var')
                warning('No data-type selected for event-alignment, defaulting to raw fluorescence C\n')
                data_type = 'C';
            end
            
            switch lower(data_type)
                case 'c'
                    data2align = obj.C;
                case 'c_smooth'
                    data2align = obj.C_smooth;
                case 'spikes'
                    data2align = obj.spikes;
                case 'c_decon'
                    data2align = obj.C_decon;
                case 'spikes_denoised'
                    data2align = obj.spikes_denoised;
                case 'spikes_conv'
                    data2align = obj.spikes_conv;
            end
            
            num_neurons = size(data2align,1);
            num_events = length(obj.event_names);
            trial_length = size(data2align,2)/obj.num_trials;
            
            reshaped_events = permute(reshape(full(obj.event_matrix),[],obj.num_trials,num_events),[1,3,2]);

            if ~isempty(trial_indices)
                temp = reshape(data2align,num_neurons,trial_length,obj.num_trials);
                temp = temp(:,trial_indices,:);
                data2align = reshape(temp,num_neurons,obj.num_trials*length(trial_indices));
                trial_length = length(trial_indices);
                reshaped_events = reshaped_events(trial_indices,:,:);
            end
            [~, scores] = pca(data2align');
            
            reshaped = permute(reshape(scores(:,1:numPCs),trial_length,obj.num_trials,numPCs),[1,3,2]);

            trial_average = mean(reshaped,3);
            
            CIs = zeros(2,trial_length,numPCs);
            
            switch lower(bound_type)
                case 'bootstrap'
                    
                    num_bootstraps = 1000;
                    
                    sample_means = zeros(trial_length,numPCs,num_bootstraps);
                    
                    % bootstrapping by running PCA on original data with
                    % different resamplings of the original 100 trials 
                    
                    reshaped_data = reshape(data2align',trial_length,obj.num_trials,num_neurons);
                    
                    sampled_trials = reshape(datasample([1:obj.num_trials],num_bootstraps*obj.num_trials),num_bootstraps,obj.num_trials);
                    
                    tic
                    for bootstrap_i = 1:num_bootstraps
                        
                        resampled_data = reshaped_data(:,sampled_trials(bootstrap_i,:),:);
                        resampled_data = reshape(resampled_data,trial_length*obj.num_trials,num_neurons);
                        [~,rs_scores] = pca(resampled_data);
                        
                        rs_scores = permute(reshape(rs_scores(:,1:numPCs),trial_length,obj.num_trials,numPCs),[1,3,2]);
                        sample_means(:,:,bootstrap_i) = mean(rs_scores,3);
                        
                    end
                    fprintf('Time taken to estimate confidence intervals on principal component projections: %.2f minutes',toc/60)
                    
                    clear sampled_trials resampled_data rs_scores;
                    
                    % bootstrapping to obtain confidence intervals on
                    % cross-trial mean PC activation -- takes less time
                    % than option above
                    
%                     bootstrap_dists = datasample(reshaped,num_bootstraps*obj.num_trials,3);
%                     bootstrap_dists_r = reshape(bootstrap_dists,[trial_length,numPCs,num_bootstraps,obj.num_trials]);
%                     sample_means = squeeze(mean(bootstrap_dists_r,4));
                    
                    lower_bound = CI_bounds(1)/100;
                    upper_bound = CI_bounds(2)/100;
                    
                    for timebin = 1:size(sample_means,1)
                        
                        if numPCs > 1
                            for pc_i = 1:numPCs
                            
                                [cd,bins] = histcounts(sample_means(timebin,pc_i,:),50);
                                cd = cumsum(cd/sum(cd));
                                
                                CIs(1,timebin,pc_i) = bins(find(cd>upper_bound,1));
                                CIs(2,timebin,pc_i) = bins(find(cd>lower_bound,1));
                            
                            end
                        else
                            [cd,bins] = histcounts(sample_means(timebin,1,:),50);
                            cd = cumsum(cd/sum(cd));
                            CIs(1,timebin,1) = bins(find(cd>upper_bound,1));
                            CIs(2,timebin,1) = bins(find(cd>lower_bound,1));
                        end
                            
                        
                    end
                    
                case 'jackknife'
                    
                    for pc_i = 1:numPCs
                        
                        [~,SDest] = jackknife(squeeze(reshaped(:,pc_i,:)));
                    
                        CIs(1,:,pc_i) = SDest;
                        CIs(2,:,pc_i) = SDest;
                        
                    end
                                        
                case 'analytic'
                    
                    for pc_i = 1:numPCs
                        SEM = std(squeeze(reshaped(:,pc_i,:)),0,2)./sqrt(obj.num_trials);
                        CIs(1,:,pc_i) = SEM;
                        CIs(2,:,pc_i) = SEM;
                    end
                    
            end
            
            event_tmspz = zeros(1,num_events);
            for event_id = 1:num_events
                [~,event_tmspz(event_id)] = max(sum(reshaped_events(:,event_id,:),3));
            end
            
            [sorted_tmspz,srt] = sort(event_tmspz,'ascend');
            
            if plot_flag
                figure('Position',[100 400 1200 400])
                xt = (1:trial_length) * (1./obj.Fs);
                
                PC_colors = winter(numPCs);
                for pc_i = 1:numPCs
                    plot(xt,trial_average(:,pc_i),'LineWidth',2,'Color',PC_colors(pc_i,:),'DisplayName',sprintf('PC%d',pc_i))
                    hold on;
                    plot(xt,CIs(1,:,pc_i),'--','LineWidth',0.25,'Color',PC_colors(pc_i,:),'DisplayName',sprintf('%.1f CI',CI_bounds(2)));
                    plot(xt,CIs(2,:,pc_i),'--','LineWidth',0.25,'Color',PC_colors(pc_i,:),'DisplayName',sprintf('%.1f CI',CI_bounds(1)));
                end
                xlim([1 75])
                ylimz = ylim;
                
                event_names = obj.event_names;
                event_names{8} = 'Foodport Entry';
                
                colors = parula(length(event_names));
                
                for event_id = 1:length(event_names)
                    plot([sorted_tmspz(event_id),sorted_tmspz(event_id)]./obj.Fs,[ylimz(1),ylimz(2)],'--','LineWidth',1.5,'Color',colors(srt(event_id),:),'DisplayName',event_names{srt(event_id)})
                end
                
                legend('show')
                ylim(ylimz)
                
                title('Principal Component Activity aligned to various experimental events')
                xlabel('Time (seconds)')
                ylabel('PC Activation')
                
                
                hold off;
            end
            
        end
        
                
    end
                               
end


