clear all; close all; clc;

%% loads
expert = load('Y:/Whiskernas/JK/suite2p/AndrewAnalysis/expert_analysis.mat');
expert = expert.processedResults;
naive = load('Y:/Whiskernas/JK/suite2p/AndrewAnalysis/naive_analysis.mat');
naive = naive.processedResults;

%% matched animal t-SNE distances

nonlearnerNum = [5,6,8,10:12];
naiveNum = [1:4,6,9];
expertNum = 1:6;
cmap = jet(7);
angles = 45:15:135;

currN = 6;
naiveY = naive(naiveNum(currN)).Y;
naiveAV = naive(naiveNum(currN)).angleVals;
expertY = expert(currN).Y;
expertAV = expert(currN).angleVals;
figure, 
subplot(121), hold on
for i = 1 : length(angles)
    scatter3(naiveY(naiveAV==angles(i),1), naiveY(naiveAV==angles(i),2), naiveY(naiveAV==angles(i),3), ...
        25, cmap(i,:), 'filled');
end
view(-135, 35);
title('Naive')
subplot(122), hold on
for i = 1 : length(angles)
    scatter3(expertY(expertAV==angles(i),1), expertY(expertAV==angles(i),2), expertY(expertAV==angles(i),3), ...
        25, cmap(i,:), 'filled');
end
view(-135, 35);
title('Expert')
%
% distance t-sne space (change Y to tableData to look in response space)
figure
subplot(121)
distance = squareform(pdist(naiveY));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);
imagesc(distance);
set(gca,'YDir','normal');
colormap('hot');
title('Naive')

subplot(122)
distance = squareform(pdist(expertY));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);
imagesc(distance);
set(gca,'YDir','normal');
colormap('hot');
title('Expert')
%%
% angle space distance
figure,
subplot(221)
angleDistance = squareform(pdist(naiveAV));
distance = squareform(pdist(naiveY));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);

scatter(angleDistance(:) + rand(length(angleDistance(:)), 1)*10-5, distance(:), 0.1,'k');
xticks([0:15:90])
xlim([-10, 100])
ylim([0, 1]);
title('Naive')
subplot(223);
boxplot(distance(:), angleDistance(:))
ylim([0, 1]);


subplot(222)
angleDistance = squareform(pdist(expertAV));
distance = squareform(pdist(expertY));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);

scatter(angleDistance(:) + rand(length(angleDistance(:)), 1)*10-5, distance(:), 0.1,'k');
xticks([0:15:90])
xlim([-10, 100])
ylim([0, 1]);
title('Expert')
subplot(224);

boxplot(distance(:), angleDistance(:))
ylim([0, 1]);

%%

%% general variables
% nCellsExpert = max([expert.nCells]);
% nCellsNaive = max([naive.nCells]);
% minCellsExpert = min([expert.nCells]);
% minCellsNaive = min([naive.nCells]);
% nAnimalsExpert = length(expert);
% nAnimalsNaive = length(naive);

nCells = min([expert.nCells, naive.nCells]);
%% average classifier performance
expertPerformance = nan(6, nCells);
naivePerformance = nan(6, nCells);
nonlearnerPerformance = nan(6, nCells);
% fitType = 'power2';
start = [-0.1, -0.5, 0.4];

for i = 1:length(expert)
   thisAccuracy = mean(expert(i).accuracy);
%    x = repmat(1:size(thisAccuracy, 1), 1, size(thisAccuracy, 2));
%    y = thisAccuracy(:);
%    
%    f = fit(x', y, fitType, 'StartPoint', start);
%    expertPerformance(i,:) = f(1:nCellsExpert);
expertPerformance(i,:) = thisAccuracy(1:nCells);
end

for i = 1:length(naiveNum)
   thisAccuracy = mean(naive(naiveNum(i)).accuracy);
%    x = repmat(1:size(thisAccuracy, 1), 1, size(thisAccuracy, 2));
%    y = thisAccuracy(:);
%    
%    f = fit(x', y, fitType, 'StartPoint', start);
%    naivePerformance(i,:) = f(1:nCellsNaive);
naivePerformance(i,:) = thisAccuracy(1:nCells);
end

for i = 1 : length(nonlearnerNum)
    thisAccuracy = mean(naive(nonlearnerNum(i)).accuracy);
    nonlearnerPerformance(i,:) = thisAccuracy(1:nCells);
end

figure; hold on;
plot(1:nCells, mean(nonlearnerPerformance), 'c')
plot(1:nCells, mean(naivePerformance), 'b')
plot(1:nCells, mean(expertPerformance), 'r')
shadedErrorBar(1:size(nonlearnerPerformance,2), mean(nonlearnerPerformance), std(nonlearnerPerformance)/sqrt(6), 'lineprops','c');
shadedErrorBar(1:size(naivePerformance,2), mean(naivePerformance), std(naivePerformance)/sqrt(6), 'lineprops','b');
shadedErrorBar(1:size(expertPerformance,2), mean(expertPerformance), std(expertPerformance)/sqrt(6), 'lineprops','r');

plot(1:nCells, mean(nonlearnerPerformance), 'c')
plot(1:nCells, mean(naivePerformance), 'b')
plot(1:nCells, mean(expertPerformance), 'r')

ylim([0.2, 0.7]);
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'southeast')
xlabel('# of cells')
ylabel('Performance')

%% example classifier performance
i = 6;
figure
subplot(121), hold on
naiveEgAcc = naive(naiveNum(i)).accuracy;
naiveEgAccShuffled = naive(naiveNum(i)).accuracyShuffled;
plot(1:size(naiveEgAcc,2), mean(naiveEgAcc), 'b')
plot(1:size(naiveEgAccShuffled,2), mean(naiveEgAccShuffled), 'k')
shadedErrorBar(1:size(naiveEgAcc,2), mean(naiveEgAcc), std(naiveEgAcc)/sqrt(10), 'lineprops','b');
shadedErrorBar(1:size(naiveEgAccShuffled,2), mean(naiveEgAccShuffled), std(naiveEgAccShuffled)/sqrt(10), 'lineprops','k');

ylim([0 0.8])
legend({'LDA', 'Shuffled'}, 'box', 'off')
xlabel('# of cells')
ylabel('Performance')
title('Naive')

subplot(122), hold on
expertEgAcc = expert(i).accuracy;
expertEgAccShuffled = expert(i).accuracyShuffled;
plot(1:size(expertEgAcc,2), mean(expertEgAcc), 'b')
plot(1:size(expertEgAccShuffled,2), mean(expertEgAccShuffled), 'k')
shadedErrorBar(1:size(expertEgAcc,2), mean(expertEgAcc), std(expertEgAcc)/sqrt(10), 'lineprops','b');
shadedErrorBar(1:size(expertEgAccShuffled,2), mean(expertEgAccShuffled), std(expertEgAccShuffled)/sqrt(10), 'lineprops','k');

ylim([0 0.8])
% legend({'LDA', 'Shuffled'}, 'box', 'off')
xlabel('# of cells')
ylabel('Performance')
title('Expert')
