%% Script September 11th, 2018


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
surround_time = [52 52]; % ~5 seconds before and after
event_nam = 'leverOUT';

display_folder = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/analysis_09112018';

%% main-loop through rats/neurons/sessions

for ii = 1:length(RatIDs)
    
    rat_id = RatIDs(ii);
    
    cell_map = load_cell_map(registration_folder,rat_id);
    
    rat_folder = sprintf('Rat%d',rat_id);
    
    if ~isdir(fullfile(display_folder,rat_folder))
        mkdir(fullfile(display_folder,rat_folder))
    end
    
    for jj = 1:size(cell_map,1)
        
        %% get relevant sessions for this neuron and set up arrays
        sess_nums = find(cell_map(jj,:));
        
        all_days_R = zeros(3,sum(surround_time)+1,length(sess_nums));
        all_days_NR = zeros(3,sum(surround_time)+1,length(sess_nums));
        
        %% run through sessions
        
        for kk = 1:length(sess_nums)
            
            sess_id = session_names{sess_nums(kk)};
            
            load(fullfile(base_directory,sprintf('rats_processed/Rat%d/Rat%d_%s.mat',rat_id,rat_id,sess_id)));
            
            event_idx = find(strcmp(event_nam,Sess_object.event_names));
            
            Sess_object.lock2events('spikes_conv',52,52)
            
            % use event_matrix structure to find all the trials where the
            % Cue light went off, signalling a successful trial
            R_trial_bool = logical(sum(full(reshape(Sess_object.event_matrix(:,find(strcmp('cueON',Sess_object.event_names))),799,Sess_object.num_trials)),1));
            NR_trial_bool = ~R_trial_bool;
            
            R_trials = Sess_object.event_locked{1}(cell_map(jj,sess_nums(kk)),:,R_trial_bool);
            if ismatrix(R_trials) || size(R_trials,3) < 10
                fprintf('Less than 10 response trials for Rat %d on session: %s\n',rat_id,sess_id);
                all_days_R(:,:,kk) = NaN;
            else
                all_days_R(1,:,kk) = mean(R_trials,3);
                
                % use bootstrap samples to compute 5% - 95% confidence
                % intervals (for the bootci_array function, use the same
                % dimension over which the trials are concatenated in order to
                % get intervals on the trial-average)
           
                all_days_R(2:end,:,kk) = bootci_array(1000,@(x)mean(x),squeeze(R_trials),2)'; % with the transposition of the output & assignment: 2nd row is lower interval, 3rd row is upper interval
            end
               
            
            NR_trials = Sess_object.event_locked{1}(cell_map(jj,sess_nums(kk)),:,NR_trial_bool);
            if ismatrix(NR_trials) || size(NR_trials,3) < 10 % don't plot if there's less than 10 trials 
                fprintf('Less than 10 non-response trials for Rat %d on session: %s\n',rat_id,sess_id);
                all_days_NR(:,:,kk) = NaN;
            else
                all_days_NR(1,:,kk) = mean(NR_trials,3);
                all_days_NR(2:end,:,kk) = bootci_array(1000,@(x)mean(x),squeeze(NR_trials),2)';
            end
            
            
        end
        
        
        %% plotting 
        
        figure('Position',[100 100 1000 700]);
        
        time_x = (-surround_time(1):surround_time(2))/10.5;
        plot_interval = 0.25;
        shifts = flipud(repmat([0:plot_interval:((length(sess_nums)-1)*plot_interval)]',1,length(time_x)));
        
        if ismatrix(all_days_R)
            y_R = all_days_R(1,:) + shifts;
            ymin_R = all_days_R(2,:) + shifts;
            ymax_R = all_days_R(3,:) + shifts;
           
        else
            y_R = squeeze(all_days_R(1,:,:))' + shifts;
            ymin_R = squeeze(all_days_R(2,:,:))' + shifts;
            ymax_R = squeeze(all_days_R(3,:,:))' + shifts;
           
        end
        
        if ismatrix(all_days_NR)
            
            y_NR = all_days_NR(1,:) + shifts;
            ymin_NR = all_days_NR(2,:) + shifts;
            ymax_NR = all_days_NR(3,:) + shifts;
        else
            
            y_NR = squeeze(all_days_NR(1,:,:))' + shifts;
            ymin_NR = squeeze(all_days_NR(2,:,:))' + shifts;
            ymax_NR = squeeze(all_days_NR(3,:,:))' + shifts;
        end
        
        
        y_upper_lim = max( max(ymax_NR(1,:)) + 0.75 * var(ymax_NR(1,:)), max(ymax_R(1,:)) + 0.75 * var(ymax_R(1,:)) );
        y_lower_lim = min( min(ymin_NR(end,:)) - 0.75 * var(ymin_NR(end,:)), min(ymin_R(end,:)) - 0.75 * var(ymin_R(end,:)) );
        
        if y_upper_lim == y_lower_lim
            y_upper_lim = y_lower_lim + 1e-10;
        end
        
        if isvector(y_R)
            g(1,1) = gramm('x',time_x,'y',y_R,'ymin',ymin_R,'ymax',ymax_R);
        else
            g(1,1) = gramm('x',time_x,'y',y_R,'ymin',ymin_R,'ymax',ymax_R,...
                'color',session_names(sess_nums)');
            g(1,1).set_order_options('color',0);
        end
        g(1,1).set_layout_options('Position',[0.02 0.05 0.5 0.8],... %Set the position in the figure (as in standard 'Position' axe property)
            'legend',false) %,... % No need to display legend for side histograms
%             'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
%             'margin_width',[0.1 0.02],...
%             'redraw',false); %We deactivate automatic redrawing/resizing so that the axes stay aligned according to the margin options
        g(1,1).axe_property('ylim',[y_lower_lim y_upper_lim]);
        g(1,1).set_names('x',sprintf('Time relative to %s (seconds)',event_nam),'y','Session');
        g(1,1).set_title('Response trials')
        g(1,1).geom_interval('geom','area');
        
        
        if isvector(y_NR)
            g(2,1) = gramm('x',time_x,'y',y_NR,'ymin',ymin_NR,'ymax',ymax_NR);
        else
            g(2,1) =gramm('x',time_x,'y',y_NR,'ymin',ymin_NR,'ymax',ymax_NR,...
                'color',session_names(sess_nums)');
            g(2,1).set_order_options('color',0);
        end

        g(2,1).geom_interval('geom','area');
        g(2,1).set_layout_options('Position',[0.5 0.05 0.5 0.8],...
            'legend_pos',[0.1 0.8 0.9 0.125]) %,... %We detach the legend from the plot and move it to the top right
%            'margin_height',[0.02 0.05],... %We set custom margins, values must be coordinated between the different elements so that alignment is maintained
%             'margin_width',[0.1 0.02],...
%             'redraw',false);
        g(2,1).axe_property('ylim',[y_lower_lim y_upper_lim]);
        g(2,1).set_names('x',sprintf('Time relative to %s (seconds)',event_nam),'y','Session');
        g(2,1).set_title('Non-response trials')

        g.set_title(sprintf('Rat %d, Neuron %d, aligned to %s',rat_id,jj,event_nam),'FontSize',16);
        g.draw();
        
        %% save image and close
        
        fig_name = [rat_folder,filesep,sprintf('Rat%d_Neuron%d_%s',rat_id,jj,event_nam),'.png'];
        
        saveas(gcf,fullfile(display_folder,fig_name));
        
        close gcf; clear g
        
        %% save data before next neuron
        
        dat_name = fullfile(rat_folder,sprintf('Rat%d_Neuron%d_%s.mat',rat_id,jj,event_nam));
        session_names_indexed = session_names{sess_nums};
        save(fullfile(display_folder,dat_name),'all_days_R','all_days_NR','surround_time','session_names_indexed','event_nam');

    end
    
end

            
            
            
            
            
            

            
            
            
            
        
        
        
        
    
    
    

    