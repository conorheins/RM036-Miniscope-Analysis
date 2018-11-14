[fnam,fdir] = uigetfile();
filename = fullfile(fdir,fnam);
data = load(filename);
ratname = fieldnames(data);
ratname = ratname{1};
A = data.(ratname){2,1}.A;
C = data.(ratname){2,1}.C;

numCells = size(C,1);

trialsframes = [reshape(repmat(1:100,799,1),799*100,1), repmat([1:799]',100,1)];

Csmooth = zeros(size(C));
aboveThr = zeros(size(C));

alltrials = unique(trialsframes(:,1));
for i = 1:numCells
    for trial = 1:length(alltrials)
        thisTrial = locsmooth(C(i,trialsframes(:,1)==alltrials(trial)),10,3);
        Csmooth(i,trialsframes(:,1)==alltrials(trial)) = thisTrial;
        thisTrial_deriv = diff(thisTrial,1);
        thrStatistic = median(thisTrial_deriv,2,'omitnan') + 4*mad(thisTrial_deriv,1,2);
        aboveThr(i,trialsframes(:,1)==alltrials(trial)) = [bsxfun(@gt,thisTrial_deriv,thrStatistic),0];
    end
end


aboveThr = aboveThr';
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
        plot(Csmooth(i,:)); hold on; plot(.75*max(Csmooth(i,:)).*crossingTimes(:,i),'r.');
        pause; hold off;
    end
end