%% Script 10.16.2018
% plot binned lever press data for each rat, at different bin-widths

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

%% initialize variables

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};
RatIDs = [8,9,10,11,15,21,23];
bin_widths = [1,2:2:20]; % how many trials do we bin by
trial_length = 799;  % 799 frames per trial

all_rats_behav =  cell(1,length(RatIDs));

binar_flag = false;

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
    
    sess_ids = [];
    for kk = 1:length(session_names)
        if sum(strcmp(exist_sessions,session_names{kk})) > 0
            sess_ids = [sess_ids,kk];
        end
    end
    
    currRat_behav = cell(2,length(sess_ids));
    
    for kk = 1:length(sess_ids)
             
        sess_name = session_names{sess_ids(kk)};
        
        currRat_behav{2,kk} = sess_name;

        load(fullfile(base_directory,fullfile(rat_folder,sprintf('Rat%d_%s.mat',rat_id,sess_name))));
        
        if binar_flag
            all_responses = full(Sess_object.event_matrix(:,strcmp(Sess_object.event_names,'cueON'))); % use cueON flags if you want to just count whether a trial is successful or not
        else
            all_responses = full(Sess_object.event_matrix(:,strcmp(Sess_object.event_names,'press')));
        end
        
        all_bin_widths = cell(1,length(bin_widths));
        
        for bin_i = 1:length(bin_widths)
            
            segment_size = trial_length*bin_widths(bin_i);
            
            num_segments = floor(length(all_responses)/segment_size);

            press_segs = [];
            for seg_i = 1:num_segments
                
                idx = (segment_size*(seg_i-1) + 1): (segment_size*seg_i);
                press_segs = [press_segs;sum(all_responses(idx))];
                
            end
            
            all_bin_widths{1,bin_i} = press_segs;
            
        end
        
        currRat_behav{1,kk} = all_bin_widths;
        
    end
    
    all_rats_behav{1,ii} = currRat_behav;
    
end
                

% for ii = 1:length(RatIDs)
%     
%     figure(ii);
%     for ext_sess = 3:6
%         plot(all_rats_behav{ii}{1,ext_sess}{2},'DisplayName',session_names{ext_sess})
%         hold on;
%     end
%     legend('show')
%     
%     title(sprintf('Extinction data for Rat %d',RatIDs(ii)));
%     
%     pause;
% end

%% 

all_rats_behav_stacked =  cell(1,length(RatIDs));

for ii = 1:length(RatIDs)
        
    current_rat_stacks = cell(1,length(bin_widths));
    
    sessions = all_rats_behav{ii}(2,:);
    
    for bin_i = 1:length(bin_widths)
        
        this_bin_data = [];
        
        for kk = 1:length(sessions)
%         for kk = 3:6 % extinction only
            
            temp_data = all_rats_behav{ii}{1,kk}{bin_i};
            this_bin_data = [this_bin_data; [kk*ones(length(temp_data),1),temp_data]];
            
        end
        
        current_rat_stacks{bin_i} = this_bin_data;
        
    end
    all_rats_behav_stacked{ii} = current_rat_stacks;
end
        

%%
for ii = 1:length(RatIDs)
    
    for bin_i = 1:length(bin_widths)
        fprintf('Bin Width: %d\n',bin_widths(bin_i))
        plot(all_rats_behav_stacked{ii}{bin_i}(:,2));
        pause;
    end
    
end

%%

writedir = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_11132018/noSmooth';

for bin2use = 1:length(bin_widths)

    for ii = 1:length(RatIDs)
        
        figure(ii)
        rat_id = RatIDs(ii);
        this_rat_dat = all_rats_behav_stacked{ii}{bin2use};
        time_x = linspace(1,400,length(this_rat_dat));
        plot(time_x,this_rat_dat(:,2),'LineWidth',0.75,'DisplayName','Lever Presses')
        hold on;
        for sess = 3:6
            next_sess_idx = find(this_rat_dat(:,1) == sess,1);
            plot([time_x(next_sess_idx), time_x(next_sess_idx)],[0, max(this_rat_dat(:,2))],'--','LineWidth',1.5,...
                'DisplayName',sprintf('Ext%d',sess-2));
        end
        legend('show')

        xlabel(sprintf('Time (%d trial bins)',bin_widths(bin2use)))
        ylabel('Lever Presses')
        title(sprintf('Extinction data for rat %d',rat_id));
        saveas(gcf,fullfile(writedir,sprintf('%dTrialBins/Rat%d_Extinction_%dTrialBins',bin_widths(bin2use),rat_id,bin_widths(bin2use))),'png')
        close gcf;
        
    end
    
end

%% now do a convolved version (with a square filter)

writedir = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_11132018/Smooth';

smoothing_widths = 2:10;

for win_width_i = 1:length(smoothing_widths)
    
    this_width = smoothing_widths(win_width_i);
    avg_kernel = ones(1,this_width)./this_width;

    for ii = 1:length(RatIDs)
        
        figure(ii)
        rat_id = RatIDs(ii);
        this_rat_dat = all_rats_behav_stacked{ii}{1};
        
        press_dat = conv(this_rat_dat(:,2)',avg_kernel,'same');
        time_x = linspace(1,400,length(press_dat));
        plot(time_x,press_dat','LineWidth',0.75,'DisplayName','Lever Presses')
        hold on;
        for sess = 3:6
            next_sess_idx = find(this_rat_dat(:,1) == sess,1);
            plot([time_x(next_sess_idx), time_x(next_sess_idx)],[0, max(press_dat)],'--','LineWidth',1.5,...
                'DisplayName',sprintf('Ext%d',sess-2));
        end
        legend('show')
        
        xlabel('Time (trials)')
        ylabel('Smoothed Lever Presses')
        title(sprintf('Extinction data for rat %d',rat_id));
        saveas(gcf,fullfile(writedir,sprintf('%dTrialBins/Rat%d_Extinction_%dTrialBins',this_width,rat_id,this_width)),'png')
        close gcf;
        
    end
    
end
        

        
        
        
        
        
        
        
                
                
            
            
        
        

    
    
        
        
    
    



