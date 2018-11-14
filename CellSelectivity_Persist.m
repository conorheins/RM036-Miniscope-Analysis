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

RatIDs = [8,9,10,11,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names
sess2align = [3:6]; % indices of which sessions to look at persistence across
surround_time = [52 52]; % ~5 seconds before and after
alpha_level = 0.05; % test at the 0.05 significance level

num_days_sig = cell(1,length(RatIDs));
num_sig_on_day = zeros(length(sess2align),length(RatIDs));
day_relative_sig = zeros(length(sess2align),length(sess2align),length(RatIDs));
day_relative_sig_null = zeros([size(day_relative_sig),1000]);

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    cell_map = load_cell_map(registration_folder,rat_id);
    
    neuron_ids = find(sum(cell_map(:,sess2align) > 0,2) == length(sess2align));
    
    neuron_sig = zeros(length(neuron_ids),2,length(sess2align));
    
    for jj = 1:length(sess2align)
        
        sess_id = session_names{sess2align(jj)};
        
        load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
        
        Sess_object.lock2events('spikes_denoised',52,52)
        Sess_object.compute_selectivity('permutation',surround_time,1000,alpha_level)
        Sess_object.compute_fidelity();
        
%         neuron_sig(:,1,jj) = Sess_object.stats_results.pVals(cell_map(neuron_ids,sess2align(jj)),2); % second column is leverIN
        neuron_sig(:,1,jj) = Sess_object.stats_results.Modulation_direction(cell_map(neuron_ids,sess2align(jj)),2); % second column is leverIN
        neuron_sig(:,2,jj) = Sess_object.stats_results.fidelity_matrix(cell_map(neuron_ids,sess2align(jj)),2); % second column is leverIN
        
    end
    
    num_days_sig_rat = zeros(length(neuron_ids),1);
    
    for neur = 1:length(neuron_ids)
        num_days_sig_rat(neur) = length(find(~isnan(squeeze(neuron_sig(neur,1,:)))));
    end
    
    num_days_sig{ii} = num_days_sig_rat;
    
    for day = 1:size(neuron_sig,3)
        num_sig_on_day(day,ii) = length(find(~isnan(squeeze(neuron_sig(:,1,day)))));
    end
    
    for day_i = 1:size(neuron_sig,3)
%         assymmetric measure
%         sig_neurons_day_i = find(~isnan(squeeze(neuron_sig(:,1,day_i))));
%         symmetric measure
         sig_neurons_day_i = ~isnan(squeeze(neuron_sig(:,1,day_i)));
        
         for day_j = 1:size(neuron_sig,3)
             
%            assymmetric measure 
%            sig_neurons_day_j = find(~isnan(squeeze(neuron_sig(:,1,day_j))));
%            day_relative_sig(day_i,day_j,ii) = length(find(ismember(sig_neurons_day_j,sig_neurons_day_i)))./length(sig_neurons_day_i);
%            symmetric measure
             sig_neurons_day_j = ~isnan(squeeze(neuron_sig(:,1,day_j)));
             day_relative_sig(day_i,day_j,ii) = pdist([sig_neurons_day_i';sig_neurons_day_j'],'jaccard');
             
             for shuff_i = 1:1000
                sig_neurons_day_j_shuffle = rand(length(sig_neurons_day_j),1) > sum(sig_neurons_day_j)/length(sig_neurons_day_j);
                day_relative_sig_shuffle(day_i,day_j,ii,shuff_i) = pdist([sig_neurons_day_i';sig_neurons_day_j_shuffle'],'jaccard');
             end
        end
            
    end
            
    
end


%% do some analysis on the accumulated significance values

day_relative_sig(isnan(day_relative_sig)) = 0;
identity_cube = logical(repmat(eye(size(day_relative_sig,1),size(day_relative_sig,2)),1,1,size(day_relative_sig,3)));
day_relative_sig(identity_cube) = 0;


upper_bound = 0.95;
lower_bound = 0.05;

CIs = zeros([size(day_relative_sig),2]);
for rat_i = 1:length(RatIDs)
    for day_i = 1:size(day_relative_sig,1)
        
        for day_j = 1:size(day_relative_sig,2)
            
            [cd,bins] = histcounts(day_relative_sig_shuffle(day_i,day_j,rat_i,:),50);
            
            cd = cumsum(cd/sum(cd));
            
            CIs(day_i,day_j,rat_i,1) = bins(find(cd>upper_bound,1));
            CIs(day_i,day_j,rat_i,2) = bins(find(cd>lower_bound,1));
            
        end
    end
end

%% plotting bullshit


for rat_i = 1:length(RatIDs)
    
    anchor_day = 1; % day to track overlap relative to
    
    figure('rend','painters','pos',[10 10 900 600])
    main_line = plot(1-day_relative_sig(anchor_day,:,rat_i),'LineWidth',2);
    hold on;
    upper_bound_line = 1-CIs(anchor_day,:,rat_i,1);
    lower_bound_line = 1-CIs(anchor_day,:,rat_i,2);
    upper_bound_line(anchor_day) = 1;
    lower_bound_line(anchor_day) = 1;
    CI_line = plot(upper_bound_line,'b--','LineWidth',1.2);
    plot(lower_bound_line,'b--','LineWidth',1.2);
    leg = legend([main_line,CI_line],{sprintf('Overlap of responsive neurons relative to %s',session_names{sess2align(anchor_day)}),'+3/-3 SEM of a shuffled overlaps'},'FontSize',14);
    legend('show')
    title(sprintf('%s - relative neuron overlap: Rat %d',session_names{sess2align(anchor_day)},RatIDs(rat_i)),'FontSize',16)
    ylabel('Jaccard similarity','FontSize',16)
    xlabel('Sesssion name','FontSize',16)
    
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', session_names(sess2align))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 14)
    
    dir_name = fullfile('Significance_over_time_23062018_ExtinctionSessions',session_names{sess2align(anchor_day)});
    if exist(dir_name) ~= 7
        mkdir(dir_name)
    end
    
    saveName = fullfile(dir_name,sprintf('Rat%d_SigNeuronOverlap_%s',RatIDs(rat_i),session_names{sess2align(anchor_day)}));
    saveas(gcf,[saveName,'.png'])
    
end
    
    

    
    

    
    


    



