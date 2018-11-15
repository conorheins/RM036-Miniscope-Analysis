function [ rate_array ] = compare_rates_acrossdays(session_objects,session_names,data_type,plot_flag)
%compare_rates_acrossdays Summary of this function goes here
%   Detailed explanation goes here

num_neurons = size(session_objects{1}.C,1);

rate_array = zeros(num_neurons,length(session_objects));

if ~exist('data_type','var') || isempty(data_type)    
    data_type = 'spikes_denoised';
end

if ~exist('plot_flag','var') || isempty(plot_flag) 
    plot_flag = false;
end

for sess = 1:length(session_objects)
    
    data_array = session_objects{sess}.(data_type);
    
    rate_array(:,sess) = sum(data_array,2)/length(data_array);

end

if plot_flag
    
    % remove outliers for plotting
    
    rate_array_2plot = rate_array;
    
    [neuron,~] = find(bsxfun(@gt,rate_array,mean(rate_array,1) + 5*std(rate_array,0,1)));
    
    rate_array_2plot(neuron,:) = [];
    
    scatter(rate_array_2plot(:,1),rate_array_2plot(:,2),'ro')
    
end

xlabel(sprintf('Firing rates on session: %s',session_names{1}))
ylabel(sprintf('Firing rates on session: %s',session_names{2}))

