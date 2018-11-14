% experimenting with using gramm library to plot trial-averaged firing
% rates

addpath(genpath('/Users/conorheins/Documents/MATLAB/piermoral_gramm_bbeffcb/'))
load('050218_TrialAvgs_allneurons.mat');

Rat_ID = 21; %rat number
sess_num = 1; % session number
Fs = 10.5;

example_spike_array = persistent_trialavg_FRs{find(RatIDs==Rat_ID)}{sess_num}.R_trials_allneur; 

[num_neurons,T,~] = size(example_spike_array);
time_x = (1:T)/Fs;

all_averages = squeeze(example_spike_array(:,:,1));
all_ymin = squeeze(example_spike_array(:,:,2));
all_ymax = squeeze(example_spike_array(:,:,3));

row_subtractor = min(all_averages,[],2);
row_divisor = max(all_averages,[],2) - min(all_averages,[],2);

all_averages = (all_averages - row_subtractor)./row_divisor;
all_ymin = (all_ymin - row_subtractor)./row_divisor;
all_ymax = (all_ymax - row_subtractor)./row_divisor;

plot_interval = 1.5;
shifts = repmat([0:plot_interval:((num_neurons-1)*plot_interval)]',1,771);

g = gramm('x',time_x,'y',all_averages + shifts,'ymin',all_ymin+shifts,'ymax',all_ymax+shifts);
g.geom_interval('geom','area');
figure('Position',[100 100 800 450]);
g.draw(); axis tight;

