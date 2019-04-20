%% for fig3e (or f or supple whatever)
%% and fig5 (learning stuff)
%% Distribution of sharpness values
baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
load('angle_tuning_summary');

%% preliminary exploration
%% all naive
% for i = 1 : 12
a = cell(12,1);
for i = 1 : 12
    a{i} = naive(i).sharpness(find(naive(i).tuned));
end

figure, histogram(cell2mat(a))

%% matching naive

matchingInd = [1:4,7,9];
b = cell(6,1);
for i = 1 : 6
    b{i} = naive(matchingInd(i)).sharpness(find(naive(matchingInd(i)).tuned));
end
figure, histogram(cell2mat(b))

%% expert
c = cell(6,1);
for i = 1 : 6
    c{i} = expert(i).sharpness(find(expert(i).tuned));
end
figure, histogram(cell2mat(c))


%% compare between nonlearner, naive, and expert
range = -0.1:0.1:1.6;

nonlearnerInd = [5,6,8,10:12];

nonlearnerHist = zeros(length(nonlearnerInd), length(range)-1);
naiveHist = zeros(length(matchingInd), length(range)-1);
expertHist = zeros(length(matchingInd), length(range)-1);

for i = 1 : length(nonlearnerInd)
    nonlearnerHist(i,:) = histcounts(naive(nonlearnerInd(i)).sharpness, range, 'normalization', 'cdf');
end
for i = 1 : length(matchingInd)
    naiveHist(i,:) = histcounts(naive(matchingInd(i)).sharpness, range, 'normalization', 'cdf');
end
for i = 1 : length(matchingInd)
    expertHist(i,:) = histcounts(expert(i).sharpness, range, 'normalization', 'cdf');
end


figure,
plot(range(2:end), mean(nonlearnerHist), 'c-'), hold on
plot(range(2:end), mean(naiveHist), 'b')
plot(range(2:end), mean(expertHist), 'r')
boundedline(range(2:end), mean(nonlearnerHist), std(nonlearnerHist)/sqrt(length(nonlearnerInd)), 'c-')
boundedline(range(2:end), mean(naiveHist), std(naiveHist)/sqrt(length(matchingInd)), 'b-')
boundedline(range(2:end), mean(expertHist), std(expertHist)/sqrt(length(matchingInd)), 'r-')
plot(range(2:end), mean(nonlearnerHist), 'c-'), hold on
plot(range(2:end), mean(naiveHist), 'b')
plot(range(2:end), mean(expertHist), 'r')
xlabel('Tuning sharpness')
ylabel('Proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off')
set(gca, 'box', 'off', 'fontname', 'myriadpro')


%% same thing, with modulation (max-min)
range = -0.1:0.1:2;

nonlearnerInd = [5,6,8,10:12];

nonlearnerHist = zeros(length(nonlearnerInd), length(range)-1);
naiveHist = zeros(length(matchingInd), length(range)-1);
expertHist = zeros(length(matchingInd), length(range)-1);

for i = 1 : length(nonlearnerInd)
    nonlearnerHist(i,:) = histcounts(naive(nonlearnerInd(i)).modulation, range, 'normalization', 'cdf');
end
for i = 1 : length(matchingInd)
    naiveHist(i,:) = histcounts(naive(matchingInd(i)).modulation, range, 'normalization', 'cdf');
end
for i = 1 : length(matchingInd)
    expertHist(i,:) = histcounts(expert(i).modulation, range, 'normalization', 'cdf');
end


figure,
plot(range(2:end), mean(nonlearnerHist), 'c-'), hold on
plot(range(2:end), mean(naiveHist), 'b')
plot(range(2:end), mean(expertHist), 'r')
boundedline(range(2:end), mean(nonlearnerHist), std(nonlearnerHist)/sqrt(length(nonlearnerInd)), 'c-')
boundedline(range(2:end), mean(naiveHist), std(naiveHist)/sqrt(length(matchingInd)), 'b-')
boundedline(range(2:end), mean(expertHist), std(expertHist)/sqrt(length(matchingInd)), 'r-')
plot(range(2:end), mean(nonlearnerHist), 'c-'), hold on
plot(range(2:end), mean(naiveHist), 'b')
plot(range(2:end), mean(expertHist), 'r')
xlabel('Tuning modulation')
ylabel('Proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off')
set(gca, 'box', 'off', 'fontname', 'myriadpro')


%% modulation vs sharpness
figure, hold on
for i = 1 : 12
    plot(naive(i).sharpness(find(naive(i).tuned)), naive(i).modulation(find(naive(i).tuned)), 'k.')
end

for i = 1 : 6
    plot(expert(i).sharpness(find(expert(i).tuned)), expert(i).modulation(find(expert(i).tuned)), 'r.')
end
axis equal
xlabel('Sharpness'), ylabel('Modulation')
xlim([0 2]), ylim([0 2])

%% sharpness median comparison between L2/3 L4 C2 non-C2
%% all naive
cd(baseDir)
info = load('cellFunctionRidgeDE010.mat');
sharpness = zeros(12,4);
for i = 1 : 12
    touchInd = find(ismember(info.naive(i).cellNums, naive(i).touchID));
    L23ind = find(info.naive(i).cellDepths < 350);
    L4ind = find(info.naive(i).cellDepths >= 350);
    C2ind = find(info.naive(i).isC2);
    nonC2ind = find(info.naive(i).isC2==0);
    touchL23C2ID = info.naive(i).cellNums(intersect(touchInd, intersect(L23ind, C2ind)));
    touchL23nonC2ID = info.naive(i).cellNums(intersect(touchInd, intersect(L23ind, nonC2ind)));
    touchL4C2ID = info.naive(i).cellNums(intersect(touchInd, intersect(L4ind, C2ind)));
    touchL4nonC2ID = info.naive(i).cellNums(intersect(touchInd, intersect(L4ind, nonC2ind)));
    sharpness(i,1) = mean(naive(i).sharpness(find(ismember(naive(i).touchID, touchL23C2ID))));
    sharpness(i,2) = mean(naive(i).sharpness(find(ismember(naive(i).touchID, touchL23nonC2ID))));
    sharpness(i,3) = mean(naive(i).sharpness(find(ismember(naive(i).touchID, touchL4C2ID))));
    sharpness(i,4) = mean(naive(i).sharpness(find(ismember(naive(i).touchID, touchL4nonC2ID))));
end
% %%
figure, 
bar(nanmean(sharpness), 'facecolor', 'w'), hold on
errorbar(nanmean(sharpness), nanstd(sharpness)/sqrt(12), 'k.')
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
ylabel('Mean sharpness')
set(gca, 'box', 'off')

%% comparison between before and after learning
naiveSharpness = zeros(6,4);
expertSharpness = zeros(6,4);

for i = 1 : 6
    ii = matchingInd(i);
    touchInd = find(ismember(info.naive(ii).cellNums, naive(ii).touchID));
    L23ind = find(info.naive(ii).cellDepths < 350);
    L4ind = find(info.naive(ii).cellDepths >= 350);
    C2ind = find(info.naive(ii).isC2);
    nonC2ind = find(info.naive(ii).isC2==0);
    touchL23C2ID = info.naive(ii).cellNums(intersect(touchInd, intersect(L23ind, C2ind)));
    touchL23nonC2ID = info.naive(ii).cellNums(intersect(touchInd, intersect(L23ind, nonC2ind)));
    touchL4C2ID = info.naive(ii).cellNums(intersect(touchInd, intersect(L4ind, C2ind)));
    touchL4nonC2ID = info.naive(ii).cellNums(intersect(touchInd, intersect(L4ind, nonC2ind)));
    naiveSharpness(i,1) = mean(naive(ii).sharpness(find(ismember(naive(ii).touchID, touchL23C2ID))));
    naiveSharpness(i,2) = mean(naive(ii).sharpness(find(ismember(naive(ii).touchID, touchL23nonC2ID))));
    naiveSharpness(i,3) = mean(naive(ii).sharpness(find(ismember(naive(ii).touchID, touchL4C2ID))));
    naiveSharpness(i,4) = mean(naive(ii).sharpness(find(ismember(naive(ii).touchID, touchL4nonC2ID))));
    
    touchInd = find(ismember(info.expert(i).cellNums, expert(i).touchID));
    L23ind = find(info.expert(i).cellDepths < 350);
    L4ind = find(info.expert(i).cellDepths >= 350);
    C2ind = find(info.expert(i).isC2);
    nonC2ind = find(info.expert(i).isC2==0);
    touchL23C2ID = info.expert(i).cellNums(intersect(touchInd, intersect(L23ind, C2ind)));
    touchL23nonC2ID = info.expert(i).cellNums(intersect(touchInd, intersect(L23ind, nonC2ind)));
    touchL4C2ID = info.expert(i).cellNums(intersect(touchInd, intersect(L4ind, C2ind)));
    touchL4nonC2ID = info.expert(i).cellNums(intersect(touchInd, intersect(L4ind, nonC2ind)));
    expertSharpness(i,1) = mean(expert(i).sharpness(find(ismember(expert(i).touchID, touchL23C2ID))));
    expertSharpness(i,2) = mean(expert(i).sharpness(find(ismember(expert(i).touchID, touchL23nonC2ID))));
    expertSharpness(i,3) = mean(expert(i).sharpness(find(ismember(expert(i).touchID, touchL4C2ID))));
    expertSharpness(i,4) = mean(expert(i).sharpness(find(ismember(expert(i).touchID, touchL4nonC2ID))));
end

figure,
bar(0.8:1:3.8, nanmean(naiveSharpness), 0.3, 'facecolor', 'w'), hold on
bar(1.2:1:4.2, nanmean(expertSharpness), 0.3, 'facecolor', 'k')
errorbar(0.8:3.8, nanmean(naiveSharpness), nanstd(naiveSharpness)/sqrt(6), 'k.')
errorbar(1.2:4.2, nanmean(expertSharpness), nanstd(expertSharpness)/sqrt(6), 'k.')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
ylabel('Mean sharpness')
legend({'Naive', 'Expert'}, 'box', 'off')
set(gca, 'box', 'off')
