% Restricted Boltzmann Machine for Rat PFC

sprintf('Add source extraction path')
addpath(genpath(uigetdir()));

sprintf('Add deconvolution path')
addpath(genpath(uigetdir()));

sprintf('Any other paths to add?')
addpath(genpath(uigetdir()));

fnam = [uigetdir(),filesep,uigetfile()];
load(fnam);

vars = whos('-file',fnam);
vars = {vars(:).name};
if any(strcmp(vars,'S'))
    already_decon = true;
    sprintf('Data already deconvolved')
else
    already_decon = false;
    sprintf('Data not yet deconvolved -- deconvolving now')
end

if ~already_decon

    C_decon = zeros(size(C_init));
    S = zeros(size(C_init));
    sn = zeros(numCells,1);

    deconvOptions.bas_nonneg = 0;
    for i = 1:numCells
        fprintf('Deconvolved Neuron No. %d of %d total\n',i,numCells);
        [C_decon(i,:), ~, ~, ~, ~, S(i,:)] = constrained_foopsi(C_init(i,:),[],[],[],[],deconvOptions);
        sn(i) = GetSn(C_init(i,:));
    end
    
end

tic
deconvOptions.bas_nonneg = 0;
for i = 1:numCells
    fprintf('Deconvolved Neuron No. %d of %d total\n',i,numCells);
    [C_decon(i,:), ~, ~, ~, ~, S(i,:)] = constrained_foopsi(C(i,:),[],[],[],[],deconvOptions);
    sn(i) = GetSn(C(i,:));
end
fprintf('Time taken to deconvolve all traces: %.2f minutes\n',(toc/60))
    
dFdT = diff(C_decon,1,2);
dFdT = [dFdT,zeros(numCells,1)];
thrStatistic1 = 5*sn;
thrStatistic2 = median(dFdT,2,'omitnan')+3.5*mad(dFdT,1,2);
aboveThr = (bsxfun(@gt,dFdT,thrStatistic2) & bsxfun(@gt,C_init,thrStatistic1))';
% aboveThr = bsxfun(@gt,dFdT,thrStatistic2)';

Fs = 10;
tooCloseWin = 1; %exclusion window (in seconds)
crossingTimes = double((diff(aboveThr)==1));
[frame, unit] = find(crossingTimes==1);
for cell = unit(1):unit(end)
    cellXTimes = frame(unit==cell);
    tooCloseInds = find(diff(cellXTimes)<=ceil(tooCloseWin*Fs))+1;
    if ~isempty(tooCloseInds)
        crossingTimes(cellXTimes(tooCloseInds),cell)=0;
    end
end

temp_in = input('Would you like to review the traces and event time-stamps? (y/n, default n)\n',s);
if strcmp(temp_in,'y')
    vis = true;
else
    vis = false;
end

if vis
    for i = 1:numCells
        plot(C_decon(i,:)); hold on; plot(.75*max(C_decon(i,:)).*crossingTimes(:,i),'r.');
        pause; hold off;
    end
end







