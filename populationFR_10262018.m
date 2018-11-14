%% Script October 3rd, 2018
% gets population average firing rate for all rats/sessions

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
dat_type = 'spikes_conv';
est_type = 'bootstrap';
CI_bounds = [16 84];

%% main-loop through rats/neurons/sessions

%Following arrays are empty
% ii = 1, kk = 4, trial_type = 1
% ii = 1, kk = [1,2], trial_type = 2
% ii = 2, kk = 4, trial_type = 1
% ii = 2, kk = [1,2], trial_type = 2
% ii = 3, kk = [1,2], trial_type = 2
% ii = 4, kk = [1,2,7], trial_type = 2
% ii = 5, kk = [4,6,7], trial_type = 1
% ii = 5, kk = [1,2,7], trial_type = 2
% ii = 6, kk = [6], trial_type = 1
% ii = 6, kk = [1,2], trial_type = 2
% ii = 7, kk = [1,2], trial_type = 2


all_popFRs = zeros(length(session_names),length(RatIDs),3,3,799); % added a third dimension for average of all trials, not separated by response vs. non-response

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    for kk = 1:length(session_names)
        
        if rat_id == 15 && kk == 7
            break
        end
        
        sess_id = session_names{kk};

        load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
        
        if exist('non_recorded_indices','var')           
            [~,~,all_popFRs(kk,ii,1,3,:),all_popFRs(kk,ii,2:end,3,:)] = Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],[],false,non_recorded_indices(1:(length(non_recorded_indices)-1)));
        else
            [~,~,all_popFRs(kk,ii,1,3,:),all_popFRs(kk,ii,2:end,3,:)] = Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],[],false);
        end
            
        % use event_matrix structure to find all the trials where the
        % Cue light went off, signalling a successful trial
        R_trial_bool = logical(sum(full(reshape(Sess_object.event_matrix(:,find(strcmp('cueON',Sess_object.event_names))),799,Sess_object.num_trials)),1));
        NR_trial_bool = ~R_trial_bool;
        
        if length(find(R_trial_bool)) >= 10
            
            if exist('non_recorded_indices','var')
                [~,~,all_popFRs(kk,ii,1,1,:),all_popFRs(kk,ii,2:end,1,:)] = ...
                    Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],R_trial_bool,false,non_recorded_indices(1:(length(non_recorded_indices)-1)));
            else
                [~,~,all_popFRs(kk,ii,1,1,:),all_popFRs(kk,ii,2:end,1,:)] = ...
                    Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],R_trial_bool,false);
            end
            
        end
        
        if length(find(NR_trial_bool)) >= 10
            
            if exist('non_recorded_indices','var')
                [~,~,all_popFRs(kk,ii,1,2,:),all_popFRs(kk,ii,2:end,2,:)] = ...
                    Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],NR_trial_bool,false,non_recorded_indices(1:(length(non_recorded_indices)-1)));
            else
                [~,~,all_popFRs(kk,ii,1,2,:),all_popFRs(kk,ii,2:end,2,:)] = ...
                    Sess_object.population_activity_wholeTrial(dat_type,est_type,CI_bounds,[],NR_trial_bool,false);
            end
            
        end
        
        clearvars -except base_directory RatIDs rat_id session_names dat_type est_type CI_bounds all_popFRs ii kk

    end
     
end

%% plotting

RatIDs = [8,9,10,11,15,21,23]; % which rats do you want to analyze
session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'}; % all session names

display_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_10040218';

time_x = (10:780)/10.5;
plot_interval = 0.05;
shifts = flipud(repmat([0:plot_interval:((length(session_names)-1)*plot_interval)]',1,length(time_x)));

load(['rats_processed',filesep,'all_FRs.mat'])

for ii = 1:length(RatIDs)
    
    figure('Position',[100 100 1000 700]);
    
    rat_folder = sprintf('Rat%d',RatIDs(ii));
    
    sess_names_R = session_names;
    sess_names_NR = session_names;

    if ~isdir(fullfile(display_folder,rat_folder))
        mkdir(fullfile(display_folder,rat_folder))
    end
    
    all_days_R = squeeze(all_popFRs(:,ii,:,1,10:780));
    y_R = squeeze(all_days_R(:,1,:)) + shifts;
    ymax_R = squeeze(all_days_R(:,2,:)) + shifts;
    ymin_R = squeeze(all_days_R(:,3,:)) + shifts;
    
    all_days_NR = squeeze(all_popFRs(:,ii,:,2,10:780));
    y_NR = squeeze(all_days_NR(:,1,:)) + shifts;
    ymax_NR = squeeze(all_days_NR(:,2,:)) + shifts;
    ymin_NR = squeeze(all_days_NR(:,3,:)) + shifts;
    
    y_upper_lim = max( max(ymax_NR(1,:)) + 0.75 * var(ymax_NR(1,:)), max(ymax_R(1,:)) + 0.75 * var(ymax_R(1,:)) );
    y_lower_lim = min( min(ymin_NR(end,:)) - 0.75 * var(ymin_NR(end,:)), min(ymin_R(end,:)) - 0.75 * var(ymin_R(end,:)) );
    
    if y_upper_lim == y_lower_lim
        y_upper_lim = y_lower_lim + 1e-10;
    end
    
    if ~isempty(find(sum(y_R - shifts,2) == 0))
        zero_rows = find(sum(y_R-shifts,2) == 0);
        y_R(zero_rows,:) = [];
        ymin_R(zero_rows,:) =[];
        ymax_R(zero_rows,:) = [];
        sess_names_R(zero_rows) = [];
    end
        
       
    g(1,1) = gramm('x',time_x,'y',y_R,'ymin',ymin_R,'ymax',ymax_R,...
        'color',sess_names_R');
    g(1,1).set_order_options('color',0);
    g(1,1).set_layout_options('Position',[0.02 0.05 0.5 0.8],... %Set the position in the figure (as in standard 'Position' axe property)
        'legend_pos',[0.1 0.8 0.9 0.125]) %,... % No need to display legend for side histograms
    %             'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    %             'margin_width',[0.1 0.02],...
    %             'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
    g(1,1).axe_property('ylim',[y_lower_lim y_upper_lim]);
    g(1,1).set_names('x','Time (seconds)','y','Session');
    g(1,1).set_title('Response trials')
    g(1,1).geom_interval('geom','area');
    
    
    if ~isempty(find(sum(y_NR - shifts,2) == 0))
        zero_rows = find(sum(y_NR-shifts,2) == 0);
        y_NR(zero_rows,:) = [];
        ymin_NR(zero_rows,:) =[];
        ymax_NR(zero_rows,:) = [];
        sess_names_NR(zero_rows) = [];
    end

    g(2,1) =gramm('x',time_x,'y',y_NR,'ymin',ymin_NR,'ymax',ymax_NR,...
        'color',sess_names_NR');
    g(2,1).set_order_options('color',0);
    g(2,1).geom_interval('geom','area');
    g(2,1).set_layout_options('Position',[0.5 0.05 0.5 0.8],...
        'legend_pos',[0.6 0.8 0.9 0.125]) %,... %We detach the legend from the plot and move it to the top right
    %            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
    %             'margin_width',[0.1 0.02],...
    %             'redraw',false);
    g(2,1).axe_property('ylim',[y_lower_lim y_upper_lim]);
    g(2,1).set_names('x','Time (seconds)','y','Session');
    g(2,1).set_title('Non-response trials')
    
    g.set_title(sprintf('Rat %d, Trial-Averaged Population Activity',RatIDs(ii)),'FontSize',16);
    g.draw();
    
    %% save image and close
    
    fig_name = [rat_folder,filesep,sprintf('Rat%d_PopAverage',RatIDs(ii)),'.png'];
    
    saveas(gcf,fullfile(display_folder,fig_name));
    
    close gcf; clear g
    
end
    
    
    
    
  

