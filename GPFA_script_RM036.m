%% RM036 -- Fitting of a Gaussian-process factor analyzer to the single trial data

clear all; clc;

fprintf('Choose base directory\n')
base_directory = uigetdir();
cd(base_directory)

addpath(genpath('singlesubject_wrapper'));
addpath(genpath('toolboxes'));

cd([base_directory,filesep,'toolboxes/OASIS_matlab-master'])
setup;
cd(base_directory);

%% Initialize instance of SubjectObj for user-chosen rat and load data

rat_id = 21;
behav_data_folder = 'Behavioral_Data';
neural_data_folder = 'NeuralData/PostMerge';
timestamps_folder = 'Timestamps';
session_name = 'SA1';

currRat = SubjectObj;

currRat.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);

%% remove bad trials, if any

% currRat.remove_trials(40); % for instance, Rat21, Ext3 had a big artifact at Trial 40

%% Process fluorescence traces and estimate spiking (deconvolution + spike-thresholding, etc.)
% deconvolve data

currRat.deconvTrials()

% denoise spikes estimated via deconvolution with a local/global
% thresholding process

dFF_thr = 6;
look_back_f = 10;
look_ahead_f = 50;
counts_flag = true;
currRat.denoiseSpikes(currRat.C,threshold_timeseries(currRat.spikes,'mean',1.5),dFF_thr,look_back_f,look_ahead_f,counts_flag)

% visualize results of spike-denoising

currRat.overlay_spikes_calc(currRat.C,double(currRat.spikes_denoised > 0))

%% GPFA fitting

most_active_id = find(sum(currRat.spikes_denoised,2) > prctile(sum(currRat.spikes_denoised,2),25));
spike_matrix = currRat.spikes_denoised(most_active_id,:);

trial_length = 799;
trial_indices = [10:780];
seq = GPFA_preprocess(spike_matrix,currRat.num_trials,trial_length,trial_indices,true);

gpfa_path = '/Users/conorheins/Desktop/CALCIUM_IMAGING/RM036/RM036_Analysis/NeuralTraj_v0300';
cd(gpfa_path);
startup;
cd(base_directory);
addpath(genpath(gpfa_path));

runIdx = 1;

result = neuralTraj(runIdx, seq,'datFormat', 'seq', ...
    'method', 'gpfa', 'binWidth',100,'xDims',8,'kernSDList', [100,200,300]);

[estParams, seqTrain] = postprocess(result, 'kernSD', [100,200,300]);

% Plot neural trajectories in 3D space
plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
% NOTES:
% - This figure shows the time-evolution of neural population
%   activity on a single-trial basis.  Each trajectory is extracted from
%   the activity of all units on a single trial.
% - This particular example is based on multi-electrode recordings
%   in premotor and motor cortices within a 400 ms period starting 300 ms 
%   before movement onset.  The extracted trajectories appear to
%   follow the same general path, but there are clear trial-to-trial
%   differences that can be related to the physical arm movement. 
% - Analogous to Figure 8 in Yu et al., J Neurophysiol, 2009.
% WARNING:
% - If the optimal dimensionality (as assessed by cross-validation in 
%   Section 2) is greater than 3, then this plot may mask important 
%   features of the neural trajectories in the dimensions not plotted.  
%   This motivates looking at the next plot, which shows all latent 
%   dimensions.

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime(seqTrain, 'xorth', 95);
% NOTES:
% - These are the same neural trajectories as in the previous figure.
%   The advantage of this figure is that we can see all latent
%   dimensions (one per panel), not just three selected dimensions.  
%   As with the previous figure, each trajectory is extracted from the 
%   population activity on a single trial.  The activity of each unit 
%   is some linear combination of each of the panels.  The panels are
%   ordered, starting with the dimension of greatest covariance
%   (in the case of 'gpfa' and 'fa') or variance (in the case of
%   'ppca' and 'pca').
% - From this figure, we can roughly estimate the optimal
%   dimensionality by counting the number of top dimensions that have
%   'meaningful' temporal structure.   In this example, the optimal 
%   dimensionality appears to be about 5.  This can be assessed
%   quantitatively using cross-validation in Section 2.
% - Analogous to Figure 7 in Yu et al., J Neurophysiol, 2009.

% ========================================================
% 2) Full cross-validation to find:
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================

% Select number of cross-validation folds
numFolds = 4;

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
xDims = [2 5 8];

% If 'parallelize' is true, all folds will be run in parallel using 
% Matlab's parfor construct. If you have access to multiple cores, this 
% provides significant speedup. 

datFormat = 'seq';

parallelize = true;
neuralTraj(runIdx, seq, 'datFormat', datFormat, 'method',  'pca', ...
    'xDims', xDims, 'numFolds', numFolds, 'parallelize', parallelize);
neuralTraj(runIdx, seq, 'datFormat', datFormat, 'method', 'ppca', ...
    'xDims', xDims, 'numFolds', numFolds, 'parallelize', parallelize);
neuralTraj(runIdx, seq, 'datFormat', datFormat, 'method',   'fa', ...
    'xDims', xDims, 'numFolds', numFolds, 'parallelize', parallelize);
neuralTraj(runIdx, seq, 'datFormat', datFormat, 'method', 'gpfa', ...
    'xDims', xDims, 'numFolds', numFolds, 'parallelize', parallelize);

fprintf('\n');
% NOTES:
% - These function calls are computationally demanding.  Cross-validation 
%   takes a long time because a separate model has to be fit for each 
%   state dimensionality and each cross-validation fold.

% Plot prediction error versus state dimensionality.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
kernSD = 100; % select kernSD for two-stage methods
plotPredErrorVsDim(runIdx, kernSD);
% NOTES:
% - Using this figure, we i) compare the performance (i.e,,
%   predictive ability) of different methods for extracting neural
%   trajectories, and ii) find the optimal latent dimensionality for
%   each method.  The optimal dimensionality is that which gives the
%   lowest prediction error.  For the two-stage methods, the latent
%   dimensionality and smoothing kernel width must be jointly
%   optimized, which requires looking at the next figure.
% - In this particular example, the optimal dimensionality is 5. This
%   implies that, even though the raw data are evolving in a
%   53-dimensional space (i.e., there are 53 units), the system
%   appears to be using only 5 degrees of freedom due to firing rate
%   correlations across the neural population.
% - Analogous to Figure 5A in Yu et al., J Neurophysiol, 2009.

% Plot prediction error versus kernelSD.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
xDim = 5; % select state dimensionality
plotPredErrorVsKernSD(runIdx, xDim);
% NOTES:
% - This figure is used to find the optimal smoothing kernel for the
%   two-stage methods.  The same smoothing kernel is used for all units.
% - In this particular example, the optimal standard deviation of a
%   Gaussian smoothing kernel with FA is 30 ms.
% - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.
