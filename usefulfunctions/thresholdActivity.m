function [thresholded,frameIndices] = thresholdActivity(x,pVal,precision,shuffleMethod)
% Use shuffled version of data to determine high activity frames
% at user-given statistical threshold and with user-given precision
%   Detailed explanation goes here

if nargin < 4
    shuffleMethod = 'time';
end

if ~any(strcmp(shuffleMethod, {'frames','time','time_shift','isi','cell','exchange'}))
    warning('No shuffle method provided: defaulting to shuffling time to preserve firing rates but eliminate correlations')
    method = 'time';
end

if nargin < 3
    precision = 0.1;
end

if nargin < 2
    pVal = 0.001;
end

if size(x,1) > size(x,2) % always assumes less neurons than timesteps (less variables than observations)
    data = x';
else
    data = x;
end

numIterations = 1/precision; % number of permutations depends on desired precision/conservativeness

criticalVals = zeros(numIterations,1);
for i = 1:numIterations
    shuffled = shuffle(data,shuffleMethod);
    distribution = sort(sum(shuffled));
    chance = ceil((1-pVal)*length(distribution));
    criticalVals(i) = distribution(chance);
end

threshold = mean(criticalVals);
frameIndices = find(sum(data)>threshold);

if size(x,1) > size(x,2)
    thresholded = data(:,frameIndices)';
else
    thresholded = data(:,frameIndices);
end

end

