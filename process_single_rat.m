function [ daily_rat_obj ] = process_single_rat(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder,params)
%process_single_rat Wrapper function for running a bunch of methods from the SubjectObj class
%   Includes loading data from a specific RatID and session, smoothing the
%   data, deconvolving it, and doing spike detection on it, as well
%   re-convolving it with a specific kernel to generate instantaneous
%   firing rates


%% load data

daily_rat_obj = SubjectObj;
daily_rat_obj.load_data(rat_id,session_name,behav_data_folder,neural_data_folder,timestamps_folder);

%% remove bad trials, if any

if params.quick_view
    trial_ids = daily_rat_obj.quickView();
end

if ~isempty(trial_ids)
    daily_rat_obj.remove_trials(trial_ids);
end

%% Process fluorescence traces and estimate spiking (deconvolution + spike-thresholding, etc.)

% smooth data with LLR
if params.smooth_flag
    daily_rat_obj.smooth_C(params.win_len,params.win_overlap);
end

% deconvolve data with a trial-by-trial procedure

daily_rat_obj.deconvTrials()

daily_rat_obj.find_bad_neurons(); % view all deconvolved traces and select neurons with bad results

% denoise spikes estimated via deconvolution with a local/global
% thresholding process

daily_rat_obj.denoiseSpikes(daily_rat_obj.C,threshold_timeseries(daily_rat_obj.spikes,params.threshold_params.stat_type,params.threshold_params.thr_factor),...
    params.denoising_params)

% convolve spikes
daily_rat_obj.convSpikes(daily_rat_obj.spikes_denoised,0,params.conv_kernel)

end

