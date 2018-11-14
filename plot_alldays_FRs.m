fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

%%
clear;

[fnam,fdir] = uigetfile('*.mat'); % choose file with all the firing rates aligned over days

load(fullfile(fdir,fnam))

session_names = {'SA1','SA2','Ext1','Ext2','Ext3','Ext4','Reinstatement'};

time_idx = 30:750;

allFRs = trialavg_FR(:,time_idx);
allCIs = CI_FR(:,:,time_idx);

allFRs_shift = allFRs - mean(allFRs,2);
allCIs_shift = allCIs - mean(allFRs,2);

shift_width = 0.02;

shifts = repmat([shift_width:shift_width:(length(session_names)*shift_width)]',1,799);
shifts = shifts(:,time_idx);

colors = hsv(length(event_matrices));

time_axis = time_idx/10.5;

figure('Position',[400 400 600 600]);

for i = 1:length(session_names)
    plot(time_axis,allFRs_shift(i,:) + shifts(i,:),'k-','LineWidth',1.5,'Color',colors(i,:));
    hold on; plot(time_axis,squeeze(allCIs_shift(i,1,:)) + shifts(i,:),'k--','LineWidth',0.5,'Color',colors(i,:));
    plot(time_axis,squeeze(allCIs_shift(i,2,:)) + shifts(i,:),'k--','LineWidth',0.5,'Color',colors(i,:));
end

xlabel('Time (seconds)','FontSize',12)
ylim([0.001 0.155])
xlim([time_axis(1) time_axis(end)])
yticks([shifts(:,1)])
yticklabels(session_names)

xt = get(gca, 'YTick');
set(gca, 'FontSize', 16)

title('Average firing rates across different sessions')

writenam = 'AvgFR_allsessions_persist.png';
write_dir = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Displays/Figures/Rat23';

saveas(gcf,fullfile(write_dir,writenam))
