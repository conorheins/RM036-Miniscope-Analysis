function visualize_all_sess(rat_data,session_names,trial_indices,Fs)
%visualize_all_sess Uses rat_data (cell array of session-specific
%structures) to plot daily population averages for response vs.
%non-response trials. Uses trial_indices and Fs to generate x-axis


all_sess_R_popavg = zeros(length(trial_indices),3,length(session_names));
all_sess_NR_popavg = zeros(length(trial_indices),3,length(session_names));

NR_sess_idx = [];
    
for jj = 1:length(session_names)
    
    all_sess_R_popavg(:,:,jj) = rat_data{jj}.R_trials_popavg;
    
    if ~isempty(rat_data{jj}.NR_trials_popavg)
        all_sess_NR_popavg(:,:,jj) = rat_data{jj}.NR_trials_popavg;
        NR_sess_idx = [NR_sess_idx,jj];
    end
    
end

figure(1);

R_average = squeeze(all_sess_R_popavg(:,1,:));
R_CI_high = squeeze(all_sess_R_popavg(:,2,:));
R_CI_low = squeeze(all_sess_R_popavg(:,3,:));

offset = max(range(R_average,1));

shifts = repmat([offset:offset: (offset * length(session_names))],length(trial_indices),1);

visualize_array = R_average + shifts;
ytick_locs = visualize_array(1,:);

xt = trial_indices./Fs;

plot(xt,visualize_array,'b','LineWidth',1.25);
hold on;
yticks(ytick_locs);
yticklabels(session_names);
xlabel('Time (seconds)');
   
hold on;
plot(xt,R_CI_high + shifts,'b--','LineWidth',0.75);
plot(xt,R_CI_low + shifts,'b--','LineWidth',0.75);

title('Trial-average of response trials')
    
axis tight;

figure(2);

NR_average = squeeze(all_sess_NR_popavg(:,1,NR_sess_idx));
NR_CI_high = squeeze(all_sess_NR_popavg(:,2,NR_sess_idx));
NR_CI_low = squeeze(all_sess_NR_popavg(:,3,NR_sess_idx));

offset = max(range(NR_average,1));

shifts = repmat([offset:offset: (offset * length(session_names(NR_sess_idx)))],length(trial_indices),1);

visualize_array = NR_average + shifts;
ytick_locs = visualize_array(1,:);

xt = trial_indices./Fs;

plot(xt,visualize_array,'r','LineWidth',1.25);
hold on;
yticks(ytick_locs);
yticklabels(session_names(NR_sess_idx));
xlabel('Time (seconds)');

hold on;
plot(xt,NR_CI_high + shifts,'r--','LineWidth',0.75);
plot(xt,NR_CI_low + shifts,'r--','LineWidth',0.75);

title('Trial-average of no-response trials')

axis tight;

end

