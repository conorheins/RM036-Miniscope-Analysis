
fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

%%

clear;

[fnam,fdir] = uigetfile('*.mat');

load(fullfile(fdir,fnam))

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

% allPCs = zeros(length(sess_names),799);
% allCIs = zeros(length(sess_names),2,799);

which_PC = 2;
% time_idx = 10:790;
time_idx = 1:771; % when the trial-indices have already been cut out earlier

allPCs = trialavg_PC(:,time_idx,which_PC);
allCIs = CI_PC(:,:,time_idx,which_PC);

allPCs_shift = allPCs - allPCs(:,1);
allCIs_shift = allCIs - allPCs(:,1);

shift_width = 0.65;

shifts = repmat([shift_width:shift_width:(length(session_names)*shift_width)]',1,799);
shifts = shifts(:,time_idx);

colors = hsv(length(event_matrices));

time_axis = [10:780]/10.5;

figure('Position',[400 400 600 600]);

for i = 1:length(session_names)
    plot(time_axis,allPCs_shift(i,:) + shifts(i,:),'k-','LineWidth',1.5,'Color',colors(i,:));
    hold on; plot(time_axis,squeeze(allCIs_shift(i,1,:)) + shifts(i,:),'k--','LineWidth',0.5,'Color',colors(i,:));
    plot(time_axis,squeeze(allCIs_shift(i,2,:)) + shifts(i,:),'k--','LineWidth',0.5,'Color',colors(i,:));
end

xlabel('Time (seconds)','FontSize',12)
ylim([0.1 5])
xlim([time_axis(1) time_axis(end)])
yticks([shifts(:,1)])
yticklabels(session_names)

xt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

title(sprintf('Principal component %d across different sessions',which_PC))

writenam = sprintf('PC%d_allsessions_persist_1to78SECONDS.png',which_PC);
write_dir = sprintf('/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/Rat21');

saveas(gcf,fullfile(write_dir,writenam))

% for i = 1:length(sess_names)
%     
% %     temp_name = sprintf('PC1_%s_BS_1SEM.mat',sess_names{i});
%     
% %     load(temp_name);
%     
%     allPCs(i,:) = trialavg_PC - trial_average(time_idx(1));
%     allCIs(i,:,:) = CIs - trial_average(time_idx(1));
%     
% end

% shift_width = 0.75;
% 
% shifts = repmat([shift_width:shift_width:(length(sess_names)*shift_width)]',1,799);
% 
% time_axis = time_idx/10.5;
% 
% figure(1);
% 
% for i = 1:length(sess_names)
%     plot(time_axis,allPCs(i,time_idx) + shifts(i,time_idx),'k-');
%     hold on; plot(time_axis,squeeze(allCIs(i,1,time_idx)) + shifts(i,time_idx),'k--','LineWidth',0.5);
%     plot(time_axis,squeeze(allCIs(i,2,time_idx)) + shifts(i,time_idx),'k--','LineWidth',0.5);
% end
% xlabel('Time (seconds)','FontSize',12)
% ylabel('Session','FontSize',12)
% ylim([0 6.0])
% xlim([time_axis(1) time_axis(end)])
% yticks([shifts(:,1)])
% yticklabels(sess_names)
% 
% xt = get(gca, 'YTick');
% set(gca, 'FontSize', 16)
% 
% title('Principal component activities across different sessions')