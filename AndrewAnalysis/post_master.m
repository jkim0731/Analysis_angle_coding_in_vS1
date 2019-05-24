clear all; close all; clc;

%% loads
expert = load('Y:/Whiskernas/JK/suite2p/AndrewAnalysis/expert_analysis.mat');
expert = expert.processedResults;
naive = load('Y:/Whiskernas/JK/suite2p/AndrewAnalysis/naive_analysis.mat');
naive = naive.processedResults;

%% general variables
nCellsExpert = max([expert.nCells]);
nCellsNaive = max([naive.nCells]);
minCellsExpert = min([expert.nCells]);
minCellsNaive = min([naive.nCells]);
nAnimalsExpert = length(expert);
nAnimalsNaive = length(naive);

%% average classifier performance
expertPerformance = nan(nAnimalsExpert, nCellsExpert);
naivePerformance = nan(nAnimalsNaive, nCellsNaive);
fitType = 'power2';
start = [-0.1, -0.5, 0.4];

for i = 1:nAnimalsExpert
   thisAccuracy = expert(i).accuracy';
   x = repmat(1:size(thisAccuracy, 1), 1, size(thisAccuracy, 2));
   y = thisAccuracy(:);
   
   f = fit(x', y, fitType, 'StartPoint', start);
   expertPerformance(i,:) = f(1:nCellsExpert);
end

for i = 1:nAnimalsNaive
   thisAccuracy = naive(i).accuracy';
   x = repmat(1:size(thisAccuracy, 1), 1, size(thisAccuracy, 2));
   y = thisAccuracy(:);
   
   f = fit(x', y, fitType, 'StartPoint', start);
   naivePerformance(i,:) = f(1:nCellsNaive);
end

figure; hold on;
shadedErrorBar(1:size(expertPerformance,2), mean(expertPerformance), std(expertPerformance), 'b');
shadedErrorBar(1:size(naivePerformance,2), mean(naivePerformance), std(naivePerformance), 'm');
ylim([0, 1]);

%% required ROIs to reach performance level
performanceLevels = 0.1:0.05:0.7;
roiRequiredExpert = nan(nAnimalsExpert, length(performanceLevels));
roiRequiredNaive = nan(nAnimalsNaive, length(performanceLevels));

for i = 1:nAnimalsExpert
   thisAccuracy = mean(expert(i).accuracy);
   for p = 1:length(performanceLevels)
       roiRequired = find(thisAccuracy >= performanceLevels(p), 1);
       if ~isempty(roiRequired)
        roiRequiredExpert(i, p) = roiRequired;
       end
   end
end

for i = 1:nAnimalsNaive
   thisAccuracy = mean(naive(i).accuracy);
   for p = 1:length(performanceLevels)
       roiRequired = find(thisAccuracy >= performanceLevels(p), 1);
       if ~isempty(roiRequired)
        roiRequiredNaive(i, p) = roiRequired;
       end
   end
end

figure; hold on;
errorbar(performanceLevels, nanmean(roiRequiredExpert), nanstd(roiRequiredExpert), 'b.');
errorbar((performanceLevels)+0.01, nanmean(roiRequiredNaive), nanstd(roiRequiredNaive), 'm.');

%% max true accuracy reached

maxAccuracyExpert = nan(1, nAnimalsExpert);
maxAccuracyNaive = nan(1, nAnimalsNaive);

for i = 1:nAnimalsExpert
   thisAccuracy = expert(i).accuracy; 
   maxAccuracyExpert(i) = max(max(thisAccuracy));
end

for i = 1:nAnimalsNaive
   thisAccuracy = naive(i).accuracy; 
   maxAccuracyNaive(i) = max(max(thisAccuracy));
end

figure; hold on;
bar(1:2, [mean(maxAccuracyNaive) mean(maxAccuracyExpert)]);
errorbar(1:2, [mean(maxAccuracyNaive) mean(maxAccuracyExpert)], [std(maxAccuracyNaive) std(maxAccuracyExpert)], '.k');
xlim([-1 4])

