%% Script 12.13.2018

% edit 12.14.2018

% Goal: concatenate activity vectors across multiple bins, generate linear/quadratic discriminant vectors
% based on training Fisher discriminators against response vs. non-response
% data (maybe try other classes as well).
% Then look at the geometry (e.g. correlations, etc.) between data vectors
% projected onto the discriminant subspace over time


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


%% parameters for plotting/rats to use/etc.
RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

data_type = 'spikes_conv'; 

sessions2track = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

events_of_interest = {'HLON','leverOUT','leverIN'}; 

temporal_bins = {[-5,5],[-5,5],[-5,5],[-5,5]};
Fs = 10.49; % sampling rate, in frames/second


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
    
    % find neurons that show up in the intersection of Self-Admin through
    % end of Extinction
    
    
    sess2track_ids = [];
    for kk = 1:length(sessions2track)
        sess2track_ids = [sess2track_ids,find(strcmp(sessions2track{kk},session_names(sess_ids)))];
    end
    
    absolute_idx_neurons = find(sum(cell_map(:,sess2track_ids)>0,2) == length(sess2track_ids));
    
    accum_data = cell(length(events_of_interest),1);

    
    %% accumulate data from appropriate sessions
    
    all_trial_flags = [];
  
    for kk = 1:length(sess2track_ids)
        
        sess_name = session_names{sess_ids(sess2track_ids(kk))};
        
        sess_specific_idx = cell_map(absolute_idx_neurons,sess_ids(sess2track_ids(kk)));
        
        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        reshaped_events = reshape(full(Sess_object.event_matrix),799,Sess_object.num_trials,size(Sess_object.event_matrix,2));
        
        R_trial_bool = logical(sum(squeeze(reshaped_events(:,:,strcmp('cueON',Sess_object.event_names))),1));
        
        all_trial_flags = [all_trial_flags;[sess_ids(sess2track_ids(kk))*ones(length(R_trial_bool),1), R_trial_bool']];
                
        reshaped_neur = reshape(Sess_object.(data_type)(sess_specific_idx,1:(799*Sess_object.num_trials)),[length(sess_specific_idx), 799,Sess_object.num_trials]);

        for event_i = 1:length(events_of_interest)
            
            [event_tmsp,trial_id] = find(squeeze(reshaped_events(:,:,strcmp(events_of_interest{event_i},Sess_object.event_names))));
            
            if ~isempty(event_tmsp)
                [trial_id,temp] = unique(trial_id,'first');
                event_tmsp = event_tmsp(temp);
                event_i_data = zeros(length(sess_specific_idx),length((Fs*temporal_bins{event_i}(1)):(Fs*temporal_bins{event_i}(2))),length(trial_id));
                
                for jj = 1:length(event_tmsp)
                    temp = event_tmsp(jj);
                    win_edges = round([temp + Fs*temporal_bins{event_i}(1), temp + Fs*temporal_bins{event_i}(2)]);
                    event_i_data(:,:,jj) = reshaped_neur(:,win_edges(1):win_edges(2),trial_id(jj));
                end
                
                accum_data{event_i} = cat(3,accum_data{event_i},event_i_data);
                
            end
                
        end
        
        clear Sess_object
    end
    
    %% look at decay of pattern similarity from self-administration average vector 
    
    time_win = 40:80;
%     time_win = 45:70;
        
    for event_i = 1:length(events_of_interest)
        
        %%
        
        close gcf;

        all_event_i_data = squeeze(mean(accum_data{event_i}(:,time_win,:),2));
%         all_event_i_data = squeeze(prctile(accum_data{event_i}(:,time_win,:),95,2));
 
%         tf_idf_d = calcTFIDF(all_event_i_data);
        
        data2use = all_event_i_data;
        
        training_sess = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4'}; % which sessions to fit linear discriminator on

        train_ids = [];
        for kk = 1:length(training_sess)
            train_ids = [train_ids,find(strcmp(training_sess{kk},session_names(sess_ids)))];
        end
        
        training_filter_vector = ismember(all_trial_flags(:,1),train_ids);
        
        training_data = data2use(:,training_filter_vector);     
        response_var_train = repmat({'Non Response'},length(find(training_filter_vector)),1);
        response_var_train(all_trial_flags(training_filter_vector,2) == 1) = {'Response'};
        
        
        
        % train LDA to distinguish response v. non-response trials from
        % last two extinction sessions
        
%         training_sess = {'Ext3','Ext4'}; % which sessions to fit linear discriminator on
%         train_ids = [];
%         for kk = 1:length(training_sess)
%             train_ids = [train_ids,find(strcmp(training_sess{kk},session_names(sess_ids)))];
%         end
%         
%         training_filter_vector = ismember(all_trial_flags(:,1),train_ids);
%         
%         training_data = data2use(:,training_filter_vector);
        
%         r_samples = find(all_trial_flags(training_filter_vector,2) == 1);
%         nr_samples = find(all_trial_flags(training_filter_vector,2) == 0);
%         
% %         rand_idx = randperm(length(r_samples));
%         rand_idx = randperm(length(nr_samples),length(r_samples));
%         nr_samples_to_use = nr_samples(rand_idx);
        
%         training_indices_to_use = [r_samples;nr_samples_to_use];
%         rand_idx2 = randperm(length(training_indices_to_use));
%         training_indices_to_use = training_indices_to_use(rand_idx2);
%         
%         training_data_subsampled = training_data(:,training_indices_to_use);
                

%         response_var_train = repmat({'Non Response'},length(find(training_filter_vector)),1);
%         response_var_train(all_trial_flags(training_filter_vector,2) == 1) = {'Response'};
        
%         response_var_train_subsampled = response_var_train(training_indices_to_use);
        
%         response_var = repmat({'Non Response'},size(data2use,2),1);
%         response_var(all_trial_flags(:,2) == 1) = {'Response'};

        %% train 10-fold cross-validated LDA on chosen sessions (90/10
        % split of train/validation sets) to get good hyperparameters (delta and
        % gamma -- just simple grid search) and then try to classify activity patterns from previous
        % sessions' response vs. non-response trials 
        
        score_transform = 'none';
        
        num_folds = 25;
        
        delta = exp(linspace(log(0.01),log(1),25));
        gamma = exp(linspace(log(0.1),log(1),25));
               
%         P = 0.5; % proportion of test samples, relative to (1-P) proportion of training samples
        
%         C = cv_partition(response_var,P,num_folds,{'Non Response','Response'});
        
%         C = cvpartition(response_var,'KFold',num_folds);
        C = cvpartition(response_var_train,'KFold',num_folds);

        % do cross_validation on model
%         all_accurs = crossValidate_LDA(data2use,response_var,C,delta,gamma,score_transform);
        all_accurs = crossValidate_LDA(training_data,response_var_train,C,delta,gamma,score_transform);
        
        avg_train_accur = squeeze(mean(all_accurs(:,:,1,:),4));
        avg_test_accur = squeeze(mean(all_accurs(:,:,2,:),4));
        
        [~,best_parameter_idx] = max(avg_test_accur(:));
        [best_delta_idx,best_gamma_idx] = ind2sub(size(avg_test_accur),best_parameter_idx);
        
        discr = fitcdiscr(training_data',response_var_train,'ClassNames',{'Non Response','Response'},'ScoreTransform',score_transform,...
                'Gamma',gamma(best_gamma_idx),'Delta',delta(best_delta_idx));
            
        M = mahal(discr,data2use');
        
        all_Mahals_NR = squareform(pdist([discr.Mu(1,:);data2use'],'mahalanobis',nancov(data2use')));
        all_Mahals_R = squareform(pdist([discr.Mu(2,:);data2use'],'mahalanobis',nancov(data2use')));
         
        response_var = repmat({'Non Response'},size(data2use,2),1);
        response_var(all_trial_flags(:,2) == 1) = {'Response'};
        
        R_distances_R_trials = all_Mahals_R(1,find(strcmp(response_var,'Response'))+1);
        R_times = find(strcmp(response_var,'Response'));
        R_distances_NR_trials = all_Mahals_R(1,find(strcmp(response_var,'Non Response'))+1);
        NR_times = find(strcmp(response_var,'Non Response'));
        bin_width = 50;
        
        [binned_R_Rtrials,binned_R_NRtrials,bin_centers] = bin_distances(R_distances_R_trials,R_times,R_distances_NR_trials,NR_times,bin_width);
        figure(1);
        plot(bin_centers,binned_R_Rtrials,'b-','LineWidth',2);
        hold on; plot(bin_centers,binned_R_NRtrials,'r-','LineWidth',2);
        title(sprintf('Mahalanobis Distance from mean %s-locked activity for Response Trials',events_of_interest{event_i}));
        ylabel('Mahanalogibs Distance from Response Trials Pattern')
        xlabel('Trials')
        
        set(gca,'FontSize',16)

        
        figure(1);
        plot(find(strcmp(response_var,'Response')),all_Mahals_R(1,find(strcmp(response_var,'Response'))+1),'b.','MarkerSize',20,'DisplayName','Response Trials')
        hold on; plot(find(strcmp(response_var,'Non Response')),all_Mahals_R(1,find(strcmp(response_var,'Non Response'))+1),'r.','MarkerSize',20,'DisplayName','Non Response Trials')
        title(sprintf('Mahalanobis Distance from mean %s-locked activity for Response Trials',events_of_interest{event_i}));
        legend('show')
        
        figure(2);
        plot(find(strcmp(response_var,'Response')),all_Mahals_NR(1,find(strcmp(response_var,'Response'))+1),'b.','MarkerSize',20,'DisplayName','Response Trials')
        hold on; plot(find(strcmp(response_var,'Non Response')),all_Mahals_NR(1,find(strcmp(response_var,'Non Response'))+1),'r.','MarkerSize',20,'DisplayName','Non Response Trials')
        title(sprintf('Mahalanobis Distance from mean %s-locked activity for Non-Response Trials',events_of_interest{event_i}));
        legend('show')
        
        test_sess = {'Reinstatement'}; % which sessions to test LDA performance on
        test_ids = [];
        for kk = 1:length(test_sess)
            test_ids = [test_ids,find(strcmp(test_sess{kk},session_names(sess_ids)))];
        end
        
        test_filter_vector = ismember(all_trial_flags(:,1),test_ids);
        
        test_data = data2use(:,test_filter_vector);
        
        response_var_test = repmat({'Non Response'},length(find(test_filter_vector)),1);
        response_var_test(all_trial_flags(test_filter_vector,2) == 1) = {'Response'};
        
        labels = discr.predict(test_data');
        
        accur = sum(strcmp(labels,'Non Response') == strcmp(response_var_test,'Non Response'))./length(labels); 
        
        M = mahal(discr,test_data');
        
%         distance_ratio = M(:,1)./M(:,2);
%         plot(find(strcmp(response_var_test,'Response')),distance_ratio(strcmp(response_var_test,'Response')),'b.','MarkerSize',20)
%         hold on; plot(find(strcmp(response_var_test,'Non Response')),distance_ratio(strcmp(response_var_test,'Non Response')),'r.','MarkerSize',20)
%         
%         all_Mahals_NR = squareform(pdist([discr.Mu(1,:);data2use'],'mahalanobis',nancov(data2use')));
%         all_Mahals_R = squareform(pdist([discr.Mu(2,:);data2use'],'mahalanobis',nancov(data2use')));
%         
        %%
        
%         all_models = cell(1,num_folds);
%         
%         for cv_i = 1:num_folds           
%             all_idx = [C{cv_i}{1};
%                         C{cv_i}{2}];   
%             all_models{cv_i} = fitcdiscr(training_data(:,all_idx)',response_var_train(all_idx),...
%                 'ClassNames',{'Non Response','Response'},'ScoreTransform',score_transform,...
%                 'Gamma',gamma(best_gamma_idx),'Delta',delta(best_delta_idx));           
%         end
%         
%         all_mu = zeros(length(unique(response_var_train)),size(training_data,1),length(all_models));
%         all_coefs = zeros(size(training_data,1),length(all_models));
%         all_constants = zeros(1,length(all_models));
%         
%         for cv_i = 1:length(all_models)
%             
%             all_mu(:,:,cv_i) = all_models{cv_i}.Mu;
%             
%             all_coefs(:,cv_i) = all_models{cv_i}.Coeffs(1,2).Linear;
%             
%             all_constants(cv_i) = all_models{cv_i}.Coeffs(1,2).Const;
%             
%         end
%         
%         true_mu = zeros(length(unique(response_var_train)),size(training_data,1));
%         true_mu(1,:) = mean(training_data(:,strcmp(response_var_train,'Non Response')),2); 
%         true_mu(2,:) = mean(training_data(:,strcmp(response_var_train,'Response')),2);
%         
%         average_beta = mean(all_coefs,2);
%         
%         response_var = repmat({'Non Response'},size(data2use,2),1);
%         response_var(all_trial_flags(:,2) == 1) = {'Response'};
%    
%         projections = (data2use' - mean(true_mu,1))*average_beta + mean(all_constants);
%         binar_labels = strcmp(response_var,'Response');
%         
%         x_ticks_r = find(binar_labels);
%         x_ticks_nr = find(~binar_labels);
%         
%         figure(1);
%         
%         plot(x_ticks_r,projections(binar_labels),'b.','DisplayName','Response Trials')
%         hold on; plot(x_ticks_nr,projections(~binar_labels),'r.','DisplayName','Non-response Trials');
%         xlabel('Trials')
%         ylabel('Alignment along maximally-discriminative dimension (projection magnitude)');
%         legend('show')
%         axis tight;
%         
%         discr = fitcdiscr(training_data',response_var_train,...
%                 'ClassNames',{'Non Response','Response'},'ScoreTransform',score_transform,...
%                 'Gamma',gamma(best_gamma_idx),'Delta',delta(best_delta_idx));
%             
%         pred_labels = discr.predict(data2use');
%         x_ticks_r = find(strcmp(pred_labels,'Response'));
%         x_ticks_nr = find(strcmp(pred_labels,'Non Response'));
%         
%         plot(x_ticks_r,projections(strcmp(pred_labels,'Response')),'b.')
%         hold on; plot(x_ticks_nr,projections(strcmp(pred_labels,'Non Response')),'r.')
        
        
        

        
                    
%         discr = fitcdiscr(training_data',response_var_train,'ClassNames',{'Non Response','Response'},...
%             'ScoreTransform',score_transform,'Delta',delta(best_delta_idx),'Gamma',gamma(best_gamma_idx));
        
                        
%         
%         discr = fitcdiscr(training_data',response_var_train,'ClassNames',{'Non Response','Response'},...,
%             'ScoreTransform',score_transform,'CrossVal','on');
% %         discr = fitcdiscr(training_data_subsampled',response_var_train_subsampled,'ClassNames',...
% %             {'Non Response','Response'},'ScoreTransform',score_transform);
%         discr = fitcdiscr(training_data_subsampled',response_var_train_subsampled,'ClassNames',...
%             {'Non Response','Response'},'ScoreTransform',score_transform,'OptimizeHyperparameters','auto');
%         
        % now test the classifier on neural activity from earlier, to see
        % how prediction accuracy fares during self-administration and the
        % extinction sessions
%         
%         test_sess = {'SA1','SA2','Ext1','Ext2'}; % which sessions to fit linear discriminator on
%         test_ids = [];
%         for kk = 1:length(test_sess)
%             test_ids = [test_ids,find(strcmp(test_sess{kk},session_names(sess_ids)))];
%         end
%         
%         test_filter_vector = ismember(all_trial_flags(:,1),test_ids);
%         
%         test_data = data2use(:,test_filter_vector);
%         
%         response_var_test = repmat({'Non Response'},length(find(test_filter_vector)),1);
%         response_var_test(all_trial_flags(test_filter_vector,2) == 1) = {'Response'};
%         
        % add in the training data that wasn't used to train the
        % classifier (e.g. the rest of the training data that is not in
        % 'subsampled')
        
%         training_data_extra = training_data;
%         training_data_extra(:,training_indices_to_use) = [];
%          
%         response_var_train_extra = response_var_train;
%         response_var_train_extra(training_indices_to_use)= [];
        
%         full_test_data = [test_data,training_data_extra];
%         full_test_data = test_data;
%         full_test_data = [test_data,training_data];
% 
%         
% %         full_test_response_var = [response_var_test;response_var_train_extra];
% %         full_test_response_var = response_var_test;
%         full_test_response_var = [response_var_test;response_var_train];
%         
%         
%         bins = 1:20:size(full_test_data,2);
%         bins(end) = size(full_test_data,2);
%         
%         accur = zeros(length(bins)-1,1);
%         for bin_i = 1:length(bins)-1
%             
%             [~,posterior_probs,cost] = discr.predict(full_test_data(:,bins(bin_i):bins(bin_i+1))');
%             %              true_labels = all_trial_flags(bins(bin_i):bins(bin_i+1),2);
%             true_labels = strcmp(full_test_response_var(bins(bin_i):bins(bin_i+1)),'Response');
%             pred_labels = posterior_probs(:,1) < 0.5;
%             accur(bin_i) = sum(true_labels == pred_labels)/size(posterior_probs,1);
%             
%         end
%         
%         figure(1);
%         title('Accuracy as a function of trials');
%         plot(bins(2:end),accur);
%         xlabel('Trials')
%         ylabel('Proportion correct classification');
%         axis tight;
        
        
        % look at Mahalanobis distance between individual activity vectors
        % and class means
        
%         M = discr.mahal(full_test_data');
        
        % project each data vector onto the discriminant axis and plot
        % projection value over time
%         beta_vec = discr.Coeffs(1,2).Linear;
%         
%         projections = (full_test_data' - mean(discr.Mu,1))*beta_vec;
%         binar_labels = strcmp(full_test_response_var,'Response');
%         
%         x_ticks_r = find(binar_labels);
%         x_ticks_nr = find(~binar_labels);
%         
%         figure(2);
%         
%         plot(x_ticks_r,projections(binar_labels),'b.','DisplayName','Response Trials')
%         hold on; plot(x_ticks_nr,projections(~binar_labels),'r.','DisplayName','Non-response Trials');
%         xlabel('Trials')
%         ylabel('Alignment along maximally-discriminative dimension (projection magnitude)');
%         legend('show')
%         axis tight;
%         pause;
        
    end
    
end
        
    