function [ accuracy ] = accuracyEval( all_theta,X,Y,label_num_list )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pred = predictOneVsAll(all_theta, X);
accuracy = zeros(length(label_num_list)+1,1);
        
for lab_i = 1:length(label_num_list)
    
    label_num = label_num_list(lab_i);
    label_specific_pred = pred;
    label_specific_pred(pred ~= label_num) = 0;
    label_specific_y = Y;
    label_specific_y(Y ~= label_num) = 0;    
    accuracy(lab_i) = mean(double(label_specific_pred == label_specific_y)) * 100;
    
end

accuracy(length(label_num_list)+1,1) = mean(double(pred == Y)) * 100;


end

