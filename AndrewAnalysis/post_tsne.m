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

%% distance metrics
xFit = linspace(0, 90, 90);

yFitExpert = nan(nAnimalsExpert, length(xFit));
for i = 1:length(expert)
    angleDistance = expert(i).angleDistance; distance = expert(i).distance;
    
    f = fit(angleDistance(:), distance(:), 'poly1');
    
    yFitExpert(i, :) = f(xFit);
end
figure; hold on;
shadedErrorBar(xFit, mean(yFitExpert), std(yFitExpert), 'b', 1);

yFitNaive = nan(nAnimalsNaive, length(xFit));
for i = 1:length(naive)
    angleDistance = naive(i).angleDistance; distance = naive(i).distance;
    
    f = fit(angleDistance(:), distance(:), 'poly1');
    
    yFitNaive(i, :) = f(xFit);
end
shadedErrorBar(xFit, mean(yFitNaive), std(yFitNaive), 'm', 1);

xlabel('Angle Difference');
ylabel('Distance')




