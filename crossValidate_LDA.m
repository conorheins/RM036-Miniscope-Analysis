function [ all_accurs ] = crossValidate_LDA(data,response_var,cross_val_partition,delta,gamma,score_transform)
%crossValidate_LDA:
% INPUTS: data -- training data, P x N
%         response_var -- labels, cell array of strings 'Response v. Non Response'
%         cross_val_partition -- cross-validation partition object
%         delta -- vector of delta parameters to fit with (regularizer on covariance matrix for LDA)
%         gamma -- vector of gamma parameters t fit with (threshold on linear coefficients for LDA)
%         score_transform -- score transformation for LDA scores
% OUTPUTS: all_accurs -- length(delta) x length(gamma) x 2 x
%          cross_val_partition.NumTestSets array, containing training and test
%          accuracies for every parameter combination and every fold

all_accurs = zeros(length(delta),length(gamma),2,cross_val_partition.NumTestSets);

for d_i = 1:length(delta)
    
    for g_i = 1:length(gamma)
        
        for cv_i = 1:cross_val_partition.NumTestSets
            
            %  trIdx = C{cv_i}{1};
            %  teIdx = C{cv_i}{2};
            trIdx = cross_val_partition.training(cv_i);
            teIdx = cross_val_partition.test(cv_i);
            
            % discr = fitcdiscr(data(:,trIdx)',response_var(trIdx),...
            %  'ClassNames',{'Non Response','Response'},'ScoreTransform',score_transform,...
            %  'Gamma',gamma(g_i),'Delta',delta(d_i));
            
            discr = fitcdiscr(data(:,trIdx)',response_var(trIdx),...
                'ClassNames',{'Non Response','Response'},...
                'ScoreTransform',score_transform,'Gamma',gamma(g_i),'Delta',delta(d_i));
            
            pred_labels = discr.predict(data(:,trIdx)');
            pred_labels = strcmp(pred_labels,'Response');
            true_labels = strcmp(response_var(trIdx),'Response');
            all_accurs(d_i,g_i,1,cv_i) = sum(true_labels == pred_labels)/length(true_labels);
            
            pred_labels = discr.predict(data(:,teIdx)');
            pred_labels = strcmp(pred_labels,'Response');
            true_labels = strcmp(response_var(teIdx),'Response');
            all_accurs(d_i,g_i,2,cv_i) = sum(true_labels == pred_labels)/length(true_labels);
            
        end
    end
end

end

