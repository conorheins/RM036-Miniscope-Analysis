%% Script September 19, 2018
% last edit October 7th 2018
% edit October 13th 2018
% edit October 14th 2018
% decoding analysis -- using neural activity to decode behavioral data

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

% decoding analysis parameters
event_names = {'HLON','leverOUT','leverIN','HLOFF'}; 
use_null_label = 'n'; % include a 'null' label in the model as a control for 'non-event' activity?
window_size = 10; % length of future (in frames) to include when aggregating neural activity into prediction vectors
data_type = 'spikes_conv'; 
train_proportion = 0.7; % proportion of the data to use as training data -- 1 - train_proportion will be used as test data
lambda_vec = logspace(-6,1,20); % regularization parameters (optimize with cross-validation)
num_folds = 5;

%% main-loop through rats/neurons/sessions

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    all_coefs = cell(length(session_names),2); % one column for the thetas, another column for the lambda value
    all_models_acc = zeros(length(event_names),length(session_names));

    %% TWO THINGS TO TRY
    % 1. BUILD MODELS USING ALL NEURONS FOR THAT SESSION, BUT THEN ONLY TRACK COEFFICIENTS OF CONSISTENT NEURONS
  
    for kk = 1:length(session_names)
        
        sess_id = session_names{kk};
        
        load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
        
        % set up accuracy array to store accuracy (per event as well as for
        % full multi-class classification), for both training and test
        % sets, for every value of lambda
        
        all_accuracies = zeros(length(event_names)+1,2,num_folds,length(lambda_vec));
        
        for fold_i = 1:num_folds
            [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,[],train_proportion);
            for jj = 1:length(lambda_vec)
                
                lambda = lambda_vec(jj);
                fprintf('\nTraining One-vs-All Logistic Regression...\n')
                [all_theta] = oneVsAll(trainX, trainY, length(label_num_list), lambda);
                [all_accuracies(:,1,fold_i,jj)] = accuracyEval(all_theta,trainX,trainY,label_num_list);
                [all_accuracies(:,2,fold_i,jj)] = accuracyEval(all_theta,testX,testY,label_num_list);
                
            end
        end
        
        %       Plot classification accuracy (averaged over folds) on train and test sets
        %         plot_CV_curves(all_accuracies,label_num_list,label_names,lambda_vec);
        
        for event_i = 1:length(event_names)
            event_i_acc = squeeze(mean(all_accuracies(strcmp(label_names,event_names{event_i}),2,:,:),3));
            all_models_acc(label_num_list(event_i),kk) = max(event_i_acc);
        end
            
        fullModel_testAcc = squeeze(mean(all_accuracies(end,2,:,:),3));
        [~,idx] = max(fullModel_testAcc);
        best_lambda = lambda_vec(idx);
        
        [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,[],train_proportion);
        allX = [trainX;testX]; allY = [trainY;testY];
        final_coefs = oneVsAll(allX,allY,length(label_num_list),best_lambda);
        all_coefs{kk,1} = final_coefs; all_coefs{kk,2} = best_lambda;
        
    end
        
    for kk = 1:length(session_names)
        
        sess2track = kk:(min(kk+1,length(session_names))); % which adjacent sessions to compare
        common_neurons = sum(cell_map(:,sess2track)>0,2) == length(sess2track);
        
        coefs2compare = zeros(length(event_names),length(find(common_neurons)),length(sess2track));
        
        for sess_i = 1:length(sess2track)
            
            neurons_of_interest = cell_map(common_neurons,sess2track(sess_i));
            
            coefs2compare(:,:,sess_i) = all_coefs{sess2track(sess_i),1}(:,neurons_of_interest+1);
            
        end
        
        for event_i = 1:length(label_names)
            figure;
            scatter(squeeze(coefs2compare(strcmp(label_names,event_names{event_i}),:,1)),squeeze(coefs2compare(2,:,2)));
            xlabel(sprintf('Regression coefficients from session: %s',session_names{sess2track(1)}));
            ylabel(sprintf('Regression coefficients from session: %s',session_names{sess2track(2)}));
            title(sprintf('Adjacent session predictor comparison for event %s',event_names{event_i}));
        end
        
        pause;
        close gcf;
        
    end
    
    save(sprintf('Rat%d_LRModels_%s',rat_id,datestr(now,'mmddyy')),'all_coefs','all_models_acc','label_names')

               
end


% 2. BUILD MODELS USING ONLY NEURONS THAT PERSIST, AND SEE HOW THEY
%    COMPARE TO MODELS WITH ANY NEURONS FROM THAT SESSION

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    all_coefs = cell(length(session_names),2); % one column for the thetas, another column for the lambda value
    
    for kk = 1:length(session_names) % loop over all sessions
        
  
        sess2track = kk:(min(kk+1,length(session_names))); % which adjacent sessions to compare
        common_neurons = sum(cell_map(:,sess2track)>0,2) == length(sess2track);
        persist_models = zeros(length(sess2track),length(event_names),length(find(common_neurons)));
        persist_models_acc = zeros(length(sess2track),1);
        random_models = zeros(length(sess2track),length(event_names),length(find(common_neurons)));
        random_models_acc = zeros(length(sess2track),1);
        
        for sess_i = 1:length(sess2track)
            
            sess_id = session_names{sess2track(sess_i)};
            
            neurons_of_interest = cell_map(common_neurons,sess2track(sess_i));
            
            %% first, train a model using the persistently-tracked neurons
            load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
            
            all_accuracies = zeros(length(event_names)+1,2,num_folds,length(lambda_vec));
            
            for fold_i = 1:num_folds
                [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,neurons_of_interest,train_proportion);
                for jj = 1:length(lambda_vec)
                    
                    lambda = lambda_vec(jj);
                    fprintf('\nTraining One-vs-All Logistic Regression...\n')
                    [all_theta] = oneVsAll(trainX, trainY, length(label_num_list), lambda);
                    [all_accuracies(:,1,fold_i,jj)] = accuracyEval(all_theta,trainX,trainY,label_num_list);
                    [all_accuracies(:,2,fold_i,jj)] = accuracyEval(all_theta,testX,testY,label_num_list);
                    
                end
            end
            
            leverOUT_testAcc = squeeze(mean(all_accuracies(strcmp(label_names,'leverOUT'),2,:,:),3));
            [best_acc,idx] = max(leverOUT_testAcc);
            best_lambda = lambda_vec(idx);
            
            [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,neurons_of_interest,train_proportion);
            allX = [trainX;testX]; allY = [trainY;testY];
            final_coefs = oneVsAll(allX,allY,length(label_num_list),best_lambda);
            persist_models(sess_i,:,:) = final_coefs;
            persist_models_acc(sess_i) = best_acc;
            
            %% next, train a model using a random subset of neurons (equal in
            %number to the persistent neurons for this pair of sessions)

            random_neurons = randperm(size(Sess_object.C,1));
            random_neurons = random_neurons(1:length(neurons_of_interest));
            
            all_accuracies = zeros(length(event_names)+1,2,num_folds,length(lambda_vec));
            
            for fold_i = 1:num_folds
                [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,random_neurons,train_proportion);
                for jj = 1:length(lambda_vec)
                    
                    lambda = lambda_vec(jj);
                    fprintf('\nTraining One-vs-All Logistic Regression...\n')
                    [all_theta] = oneVsAll(trainX, trainY, length(label_num_list), lambda);
                    [all_accuracies(:,1,fold_i,jj)] = accuracyEval(all_theta,trainX,trainY,label_num_list);
                    [all_accuracies(:,2,fold_i,jj)] = accuracyEval(all_theta,testX,testY,label_num_list);
                    
                end
            end
            
            leverOUT_testAcc = squeeze(mean(all_accuracies(strcmp(label_names,'leverOUT'),2,:,:),3));
            [best_acc,idx] = max(leverOUT_testAcc);
            best_lambda = lambda_vec(idx);
            
            [trainX,trainY,testX,testY,label_names,label_num_list] = traintest_divide(event_names,use_null_label,window_size,Sess_object,data_type,random_neurons,train_proportion);
            allX = [trainX;testX]; allY = [trainY;testY];
            final_coefs = oneVsAll(allX,allY,length(label_num_list),best_lambda);
            random_models(sess_i,:,:) = final_coefs;
            random_models_acc(sess_i) = best_acc;

                      
        end
        
    end
    
end
            
            
            
            


%% OLD STUFF
        
%         %if you wanna do multinomial logistic regression, use the code
%         %below -- but beware of errors like the following:
%         %         Error using linsolve
%         % Matrix must be positive definite.
%         % Error in mnrfit (line 248)
%         %             bcov = linsolve(hess,eye(size(hess)),struct('SYM',true,'POSDEF',true));
% %         [B,DEV,STATS] = mnrfit(fullX,fullY);
%         
%         % check out the following: https://de.mathworks.com/help/stats/fitcecoc.html
%         % multi-class classification using an SVM
%         
%         class_labels = cell(length(fullY),1);
%         for label_i = 1:length(label_names)
%             class_labels(fullY == label_i) = label_names(label_i);         
%         end
%         
%         % some rudimentary 'feature selection'
%         highestPred = find(sum(fullX,1) >= prctile(sum(fullX,1),25));
%         
%         
%         t = templateSVM('Standardize',1);
%         Mdl = fitcecoc(fullX(:,highestPred),class_labels,'Learners',t,'ClassNames',label_names);
%         CVMdl = crossval(Mdl);
%         loss = kfoldLoss(CVMdl);
%         disp(loss)
%         
%         % do some binary classification by turning label of everything except 
%         % class of interest into 'null' labels
%         
%         event2pred = 'leverOUT';
%         pos_inds = strcmp(class_labels,event2pred);
%         neg_inds = ~pos_inds;
%         
%         new_labels = class_labels;
%         new_labels(neg_inds) = {'NULL'};
%         
%         SVMModel = fitcsvm(fullX(:,highestPred),new_labels,'Standardize',true,'KernelFunction','linear',...
%             'KernelScale','auto');
%         CVSVMModel = crossval(SVMModel);
%         classLoss = kfoldLoss(CVSVMModel);
%         
%         Mdl = fitcsvm(fullX(:,highestPred),new_labels,'OptimizeHyperparameters','auto',...
%             'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%             'expected-improvement-plus'));
%         
%         CVSVMModel = crossval(Mdl);
%         classLoss = kfoldLoss(CVSVMModel);
%         
%         
%         Lambda = logspace(-6,1,20);
%         CVMdl = fitclinear(fullX(:,highestPred),new_labels,'KFold',5,...
%             'Learner','logistic','Solver','sparsa','Regularization','lasso',...
%             'Lambda',Lambda,'GradientTolerance',1e-8);
%         ce = kfoldLoss(CVMdl);
%         
%         Mdl = fitclinear(fullX(:,highestPred),new_labels,...
%             'Learner','logistic','Solver','sparsa','Regularization','lasso',...
%             'Lambda',Lambda,'GradientTolerance',1e-8);
%         numNZCoeff = sum(Mdl.Beta~=0);
%         
%         figure;
%         [h,hL1,hL2] = plotyy(log10(Lambda),log10(ce),...
%             log10(Lambda),log10(numNZCoeff));
%         hL1.Marker = 'o';
%         hL2.Marker = 'o';
%         ylabel(h(1),'log_{10} classification error')
%         ylabel(h(2),'log_{10} nonzero-coefficient frequency')
%         xlabel('log_{10} Lambda')
%         title('Test-Sample Statistics')
%         hold off
    



            
            
            
            
            
            

            
            
            
            
        
        
        
        
    
    
    

    