function plot_CV_curves( all_accuracies,label_num_list,label_names,lambda_vec )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num_labels = length(label_num_list);

figure(1);

for lab_i = 1:num_labels
    
    subplot(ceil(num_labels/2),2,lab_i)
    [h,hL1,hL2] = plotyy(log10(lambda_vec),squeeze(mean(all_accuracies(lab_i,1,:,:),3)),...
        log10(lambda_vec),squeeze(mean(all_accuracies(lab_i,2,:,:),3)));
    hL1.Marker = 'o';
    hL2.Marker = 'o';
    ylabel(h(1),'Training Set Accuracy')
    ylabel(h(2),'Test Set Accuracy')
    xlabel('Lambda')
    xlabel('log_{10} Lambda')
    title(sprintf('Cross-Validation Curve for predicting %s',label_names{label_num_list(lab_i)}));
    hold off
    
end

figure(2);

total_idx = length(label_num_list)+1;
[h,hL1,hL2] = plotyy(log10(lambda_vec),squeeze(mean(all_accuracies(total_idx,1,:,:),3)),...
    log10(lambda_vec),squeeze(mean(all_accuracies(total_idx,2,:,:),3)));
hL1.Marker = 'o';
hL2.Marker = 'o';
ylabel(h(1),'Training Set Accuracy')
ylabel(h(2),'Test Set Accuracy')
xlabel('Lambda')
xlabel('log_{10} Lambda')
title('Cross-Validation Curve for multi-class prediction');
hold off;

end

