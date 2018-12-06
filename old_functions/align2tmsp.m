function aligned = align2tmsp(data,trialsframes,eventTmsp,summaryType,preT,postT,baselineCorr)
% align2tmsp: Aligns temporal traces from multiple channels to a
% user-provided event timestamp (relative to trial-start), and either returns all traces aligned, 
% or a summary statistic of the aligned traces for each channel.
%
% Conor Heins, 05172017, conor.heins@gmail.com
%
% INPUTS:
%   -data: matrix of T x N, where N is number of channels/neurons, T is number of
%   temporal samples.
%   -trialsframes: Trial and Frame array of size T x 2, where the first
%   column contains the trial number and the second column contains the
%   timestep number.
%   -eventTmsp: the event timestamp of interest (positive integer with
%   maximum value that's the length of the trial)
%   -summaryType: a string variable determining the nature of the output, time-aligned trace(s). 
%   Possible values: {[] (default), 'avg', 'stdev', 'median'}
%   -preT: number of frames before eventTmsp to use for display/for baseline correction (if
%   baselineCorr = 1)
%   -postT: numbers of frames after eventTmsp to use for display
%   baselineCorr: set to 1 to normalize traces by activity over preT frames, 0 for no
%   normalization
%
% OUTPUTS:
%   -aligned: a 2D ( {length(preT:postT)+1 x N} in the case of summaryType being
%   non-empty) or 3D ( {length(preT:postT)+1 x N x Trials} in the case of
%   summaryType being empty) matrix of event-aligned traces

T = size(data,1);
N = size(data,2);
numTrials = length(unique(trialsframes(:,1)));

if ~exist('baselineCorr', 'var') || isempty(baselineCorr) % default: no baseline normalization
    baselineCorr = 0; 
end

if ~exist('postT','var') || isempty(postT) % default: display goes until the end of the trial
    postT = max(trialsframes(:,2))-eventTmsp;
end

if ~exist('preT','var') || isempty(preT) % default: display starts from the beginning of trial
    preT = eventTmsp - 1;
end

if ~exist('summaryType','var') || isempty(summaryType) % default: no summary statistic, aligns all traces and returns all individual traces 
    summaryType = 'None';   
end

% max normalize data
data = bsxfun(@rdivide,bsxfun(@minus,data,min(data)),(max(data)-min(data)));

eventTimes = find(trialsframes(:,2)==eventTmsp);

temp = zeros(preT+postT+1,N,numTrials);

if logical(baselineCorr)
    for i = 1:numTrials
        preChunk = eventTimes(i)-preT;
        postChunk = eventTimes(i)+postT;
        
        mu = mean(data(preChunk:eventTimes(i)-1,:),1); %baseline correction
        sigma = std(data(preChunk:eventTimes(i)-1,:),[],1);
        sigma(sigma==0)=1;
  
        temp(:,:,i) = bsxfun(@rdivide,bsxfun(@minus,data(preChunk:postChunk,:),mu),sigma);
    end
else
    for i = 1:length(eventTimes)
        preChunk = eventTimes(i)-preT;
        postChunk = eventTimes(i)+postT;
        temp(:,:,i) = data(preChunk:postChunk,:);
    end
end

switch summaryType
    case 'None'
        aligned = temp;
    case 'avg'
        aligned = mean(temp,3);
    case 'stdev'
        aligned = std(temp,[],3);
    case 'median'
        aligned = median(temp,[],3);
end
        
end

