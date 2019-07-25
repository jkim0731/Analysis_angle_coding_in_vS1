%%
% Check proportion of 90 degrees pole response touch cells and whisking
% cells, in or outside of C2, in L2/3. Compare with peron et al.

% Check consistency between glm, calcium anova, and spike anova data.

% See proportion of cell functions

% Compare between naive and expert


%%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3], [4], [4], [4], [4]}; 
% for mi = 1 : length(mice)
%     for si = 1 : length(sessions{mi})
%         Uber.buildUberArray(mice(mi), sessions{mi}(si))
%     end
% end



%%
% clear
tic
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);


for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);
    load(sprintf('JK%03dS%02dglm_cell_function_lasso_NC',mouse,session))
    naive(ni) = glm;
end

for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
%     expert(ei) = glm_results_cell_function(mouse, session, baseDir);
    load(sprintf('JK%03dS%02dglm_cell_function_lasso_NC',mouse,session))
    expert(ei) = glm;
end
% 
% % L4mice = [70,74,75,76];
% % L4sessions = [6,4,4,4];
% L4mice = [70];
% L4sessions = [6];
% 
% % L4 = struct;
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_cell_function(mouse, session, baseDir);
% end

save('Y:\Whiskernas\JK\suite2p\cellFunctionLasso_NC.mat', 'naive', 'expert')
toc


%% changing one of glm results

baseDir = 'D:\TPM\JK\suite2p\';

load([baseDir, 'cellFunctionRidgeDE010'], 'naive', 'expert')
expert(2) = glm_results_cell_function(27,10,baseDir);
save('Y:\Whiskernas\JK\suite2p\cellFunctionRidgeDE010_JK027_S10.mat', 'naive', 'expert')

expert(2) = glm_results_cell_function(27,9,baseDir);
save('Y:\Whiskernas\JK\suite2p\cellFunctionRidgeDE010_JK027_S09.mat', 'naive', 'expert')


%% Scnn1a mice

baseDir = 'D:\TPM\JK\suite2p\';
L4 = struct;
L4(1) = glm_results_cell_function(75,4,baseDir);
L4(2) = glm_results_cell_function(76,4,baseDir);

save([baseDir, 'Scnn1a_ridgeDE010'], 'L4')

%% Comparing between glm methods

lasso010 = load('cellFunctionLassoDE010.mat');
lasso005 = load('cellFunctionLassoDE005.mat');

ridge010 = load('cellFunctionRidgeDE010.mat');
ridge005 = load('cellFunctionRidgeDE005.mat');
%%
whiskingPropNaive = zeros(12,4); % 1 lasso 010, 2 lasso 005, 3 ridge 010, 4 ridge 005
touchPropNaive = zeros(12,4);
whiskingPropExpert = zeros(6,4); % 1 lasso 010, 2 lasso 005, 3 ridge 010, 4 ridge 005
touchPropExpert = zeros(6,4);
whiskingPropMatchingNaive = zeros(6,4); % 1 lasso 010, 2 lasso 005, 3 ridge 010, 4 ridge 005
touchPropMatchingNaive = zeros(6,4);

for i = 1 : 12
    whiskingPropNaive(i,1) = length(lasso010.naive(i).whiskingID) / length(lasso010.naive(i).cellNums);
    whiskingPropNaive(i,2) = length(lasso005.naive(i).whiskingID) / length(lasso005.naive(i).cellNums);
    whiskingPropNaive(i,3) = length(ridge010.naive(i).whiskingID) / length(ridge010.naive(i).cellNums);
    whiskingPropNaive(i,4) = length(ridge005.naive(i).whiskingID) / length(ridge005.naive(i).cellNums);
    
    touchPropNaive(i,1) = length(lasso010.naive(i).touchID) / length(lasso010.naive(i).cellNums);
    touchPropNaive(i,2) = length(lasso005.naive(i).touchID) / length(lasso005.naive(i).cellNums);
    touchPropNaive(i,3) = length(ridge010.naive(i).touchID) / length(ridge010.naive(i).cellNums);
    touchPropNaive(i,4) = length(ridge005.naive(i).touchID) / length(ridge005.naive(i).cellNums);
    
end

for i = 1 : 6
    whiskingPropExpert(i,1) = length(lasso010.expert(i).whiskingID) / length(lasso010.expert(i).cellNums);
    whiskingPropExpert(i,2) = length(lasso005.expert(i).whiskingID) / length(lasso005.expert(i).cellNums);
    whiskingPropExpert(i,3) = length(ridge010.expert(i).whiskingID) / length(ridge010.expert(i).cellNums);
    whiskingPropExpert(i,4) = length(ridge005.expert(i).whiskingID) / length(ridge005.expert(i).cellNums);
    
    touchPropExpert(i,1) = length(lasso010.expert(i).touchID) / length(lasso010.expert(i).cellNums);
    touchPropExpert(i,2) = length(lasso005.expert(i).touchID) / length(lasso005.expert(i).cellNums);
    touchPropExpert(i,3) = length(ridge010.expert(i).touchID) / length(ridge010.expert(i).cellNums);
    touchPropExpert(i,4) = length(ridge005.expert(i).touchID) / length(ridge005.expert(i).cellNums);    
    
end
%%
mean(whiskingPropExpert)
mean(touchPropExpert)
%%
matchingNi = [1,2,3,4,7,9];
mean(whiskingPropNaive(matchingNi,:))
mean(touchPropNaive(matchingNi,:))


%%
figure, hold all

bar(1:1:4,mean(whiskingPropNaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:1:4.3,mean(whiskingPropNaive(matchingNi,:)), 0.3, 'g', 'linestyle', 'none')
bar(1.6:1:4.6,mean(whiskingPropExpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1.6:1:4.6,mean(whiskingPropExpert), std(whiskingPropExpert)/sqrt(size(whiskingPropExpert,1)), 'r.', 'linewidth', 2)
errorbar(1:1:4,mean(whiskingPropNaive), std(whiskingPropNaive)/sqrt(size(whiskingPropNaive,1)), 'b.', 'linewidth', 2)
errorbar(1.3:1:4.3,mean(whiskingPropNaive(matchingNi,:)), std(whiskingPropNaive(matchingNi,:))/sqrt(size(whiskingPropNaive(matchingNi,:),1)), 'g.', 'linewidth', 2)    

xticks(1.3:1:4.3)
xticklabels({'Lasso 0.1', 'Lasso 0.05', 'Ridge 0.1', 'Ridge 0.05'})
ylabel('Proportion')
title('Whisking cells')
legend({'All naive', 'Matching naive', 'Expert'})

%%
figure, hold all

bar(1:1:4,mean(touchPropNaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:1:4.3,mean(touchPropNaive(matchingNi,:)), 0.3, 'g', 'linestyle', 'none')
bar(1.6:1:4.6,mean(touchPropExpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1.6:1:4.6,mean(touchPropExpert), std(touchPropExpert)/sqrt(size(touchPropExpert,1)), 'r.', 'linewidth', 2)
errorbar(1:1:4,mean(touchPropNaive), std(touchPropNaive)/sqrt(size(touchPropNaive,1)), 'b.', 'linewidth', 2)
errorbar(1.3:1:4.3,mean(touchPropNaive(matchingNi,:)), std(touchPropNaive(matchingNi,:))/sqrt(size(touchPropNaive(matchingNi,:),1)), 'g.', 'linewidth', 2)    

xticks(1.3:1:4.3)
xticklabels({'Lasso 0.1', 'Lasso 0.05', 'Ridge 0.1', 'Ridge 0.05'})
ylabel('Proportion')
title('Touch cells')
legend({'All naive', 'Matching naive', 'Expert'})

%%
touchCellIDnaive = zeros(12,2); % col 1 between 0.1, 2 between 0.05
whiskingCellIDnaive = zeros(12,2); % col 1 between 0.1, 2 between 0.05
touchCellIDexpert = zeros(6,2); % col 1 between 0.1, 2 between 0.05
whiskingCellIDexpert = zeros(6,2); % col 1 between 0.1, 2 between 0.05

for i = 1 : 12
    touchCellIDnaive(i,1) = length(intersect(lasso010.naive(i).touchID, ridge010.naive(i).touchID)) / length(union(lasso010.naive(i).touchID, ridge010.naive(i).touchID));
    touchCellIDnaive(i,2) = length(intersect(lasso005.naive(i).touchID, ridge005.naive(i).touchID)) / length(union(lasso005.naive(i).touchID, ridge005.naive(i).touchID));

    whiskingCellIDnaive(i,1) = length(intersect(lasso010.naive(i).whiskingID, ridge010.naive(i).whiskingID)) / length(union(lasso010.naive(i).whiskingID, ridge010.naive(i).whiskingID));
    whiskingCellIDnaive(i,2) = length(intersect(lasso005.naive(i).whiskingID, ridge005.naive(i).whiskingID)) / length(union(lasso005.naive(i).whiskingID, ridge005.naive(i).whiskingID));
end

for i = 1 : 6
    touchCellIDexpert(i,1) = length(intersect(lasso010.expert(i).touchID, ridge010.expert(i).touchID)) / length(union(lasso010.expert(i).touchID, ridge010.expert(i).touchID));
    touchCellIDexpert(i,2) = length(intersect(lasso005.expert(i).touchID, ridge005.expert(i).touchID)) / length(union(lasso005.expert(i).touchID, ridge005.expert(i).touchID));

    whiskingCellIDexpert(i,1) = length(intersect(lasso010.expert(i).whiskingID, ridge010.expert(i).whiskingID)) / length(union(lasso010.expert(i).whiskingID, ridge010.expert(i).whiskingID));
    whiskingCellIDexpert(i,2) = length(intersect(lasso005.expert(i).whiskingID, ridge005.expert(i).whiskingID)) / length(union(lasso005.expert(i).whiskingID, ridge005.expert(i).whiskingID));
end

%%
figure, hold all
bar(1:1:2, mean(touchCellIDnaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:1:2.3, mean(touchCellIDexpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1:1:2, mean(touchCellIDnaive), std(touchCellIDnaive)/sqrt(size(touchCellIDnaive,1)), 'b.', 'linewidth', 2)
errorbar(1.3:1:2.3, mean(touchCellIDexpert), std(touchCellIDexpert)/sqrt(size(touchCellIDexpert,1)), 'r.', 'linewidth', 2)
xticks(1.15:1:4.15)
xticklabels({'DE 0.1', 'DE 0.05'})
ylabel('Proportion')
title('Touch cells')
legend({'Naive', 'Expert'})

figure, hold all
bar(1:1:2, mean(whiskingCellIDnaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:1:2.3, mean(whiskingCellIDexpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1:1:2, mean(whiskingCellIDnaive), std(whiskingCellIDnaive)/sqrt(size(whiskingCellIDnaive,1)), 'b.', 'linewidth', 2)
errorbar(1.3:1:2.3, mean(whiskingCellIDexpert), std(whiskingCellIDexpert)/sqrt(size(whiskingCellIDexpert,1)), 'r.', 'linewidth', 2)

xticks(1.15:1:4.15)
xticklabels({'DE 0.1', 'DE 0.05'})
ylabel('Proportion')
title('Whisking cells')
legend({'Naive', 'Expert'})


%% lasso VS ridge, which is similar to 1st touch response in spike ANOVA?
baseDir = 'C:\Data\suite2p\';
mice = [25,27,30,36,39,52];
sessions = [19,16,21,17,22,21];

propAnovaLasso = zeros(length(mice), 3); % 1 touch, 2 tuned, 3 not-tuned touch
propAnovaRidge = zeros(length(mice), 3); % 1 touch, 2 tuned, 3 not-tuned touch
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions(mi);
    cd(sprintf('%s%03d',baseDir, mouse));
    anovafn = sprintf('JK%03dS%02dsingleCell_anova_spk_final.mat', mouse, session);
    load(anovafn, 'cellsTuned', 'cellsNTResponse');
    anovaTouch = union(cellsTuned, cellsNTResponse);
    propAnovaLasso(mi,1) = length(find(ismember(lasso010.expert(mi).touchID, anovaTouch))) / length(lasso010.expert(mi).touchID);
    propAnovaLasso(mi,2) = length(find(ismember(lasso010.expert(mi).touchID, cellsTuned))) / length(lasso010.expert(mi).touchID);
    propAnovaLasso(mi,3) = length(find(ismember(lasso010.expert(mi).touchID, cellsNTResponse))) / length(lasso010.expert(mi).touchID);

    propAnovaRidge(mi,1) = length(find(ismember(ridge010.expert(mi).touchID, anovaTouch))) / length(ridge010.expert(mi).touchID);
    propAnovaRidge(mi,2) = length(find(ismember(ridge010.expert(mi).touchID, cellsTuned))) / length(ridge010.expert(mi).touchID);
    propAnovaRidge(mi,3) = length(find(ismember(ridge010.expert(mi).touchID, cellsNTResponse))) / length(ridge010.expert(mi).touchID);
end
%%
figure, hold all
bar(1:1:3, mean(propAnovaLasso), 0.3, 'b', 'linestyle', 'none')
bar(1.3:1:3.3, mean(propAnovaRidge), 0.3, 'r', 'linestyle', 'none')
errorbar(1:1:3, mean(propAnovaLasso), std(propAnovaLasso)/sqrt(size(propAnovaLasso,1)), 'b.', 'linewidth', 2)
errorbar(1.3:1:3.3, mean(propAnovaRidge), std(propAnovaRidge)/sqrt(size(propAnovaRidge,1)), 'r.', 'linewidth', 2)

xticks(1.15:1:3.15)
xticklabels({'All touch', 'Tuned', 'Not-tuned'})
ylabel('Proportion')
title('Matching touch cells')
legend({'Lasso 0.1', 'Ridge 0.1'})

%% c2, non-c2, L2/3, L4. Proportions. Compare lasso 0.1 and ridge 0.1
propTouchLasso = zeros(length(mice), 4); % 1 L2/3 C2, 2 L4 C2, 3 L2/3 non-C2, 4 L4 non-C2
propTouchRidge = zeros(length(mice), 4);
propWhiskingLasso = zeros(length(mice), 4);
propWhiskingRidge = zeros(length(mice), 4);

for mi = 1 : length(mice)
    temp = lasso010.expert(mi);
    idL23C2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(temp.isC2)));
    idL4C2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(temp.isC2)));
    idL23nonC2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(1 - temp.isC2)));
    idL4nonC2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(1 - temp.isC2)));
    
    propTouchLasso(mi,1) = length(intersect(temp.touchID, idL23C2)) / length(idL23C2);
    propTouchLasso(mi,2) = length(intersect(temp.touchID, idL4C2)) / length(idL4C2);
    propTouchLasso(mi,3) = length(intersect(temp.touchID, idL23nonC2)) / length(idL23nonC2);
    propTouchLasso(mi,4) = length(intersect(temp.touchID, idL4nonC2)) / length(idL4nonC2);
    
    propWhiskingLasso(mi,1) = length(intersect(temp.whiskingID, idL23C2)) / length(idL23C2);
    propWhiskingLasso(mi,2) = length(intersect(temp.whiskingID, idL4C2)) / length(idL4C2);
    propWhiskingLasso(mi,3) = length(intersect(temp.whiskingID, idL23nonC2)) / length(idL23nonC2);
    propWhiskingLasso(mi,4) = length(intersect(temp.whiskingID, idL4nonC2)) / length(idL4nonC2);
    
    temp = ridge010.expert(mi);
    idL23C2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(temp.isC2)));
    idL4C2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(temp.isC2)));
    idL23nonC2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(1 - temp.isC2)));
    idL4nonC2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(1 - temp.isC2)));
    
    propTouchRidge(mi,1) = length(intersect(temp.touchID, idL23C2)) / length(idL23C2);
    propTouchRidge(mi,2) = length(intersect(temp.touchID, idL4C2)) / length(idL4C2);
    propTouchRidge(mi,3) = length(intersect(temp.touchID, idL23nonC2)) / length(idL23nonC2);
    propTouchRidge(mi,4) = length(intersect(temp.touchID, idL4nonC2)) / length(idL4nonC2);
    
    propWhiskingRidge(mi,1) = length(intersect(temp.whiskingID, idL23C2)) / length(idL23C2);
    propWhiskingRidge(mi,2) = length(intersect(temp.whiskingID, idL4C2)) / length(idL4C2);
    propWhiskingRidge(mi,3) = length(intersect(temp.whiskingID, idL23nonC2)) / length(idL23nonC2);
    propWhiskingRidge(mi,4) = length(intersect(temp.whiskingID, idL4nonC2)) / length(idL4nonC2);
    
    
end
%%
figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(1, mean(propTouchLasso(:,i)), 0.3, 'b', 'linestyle', 'none')
    bar(1.3, mean(propTouchRidge(:,i)), 0.3, 'r', 'linestyle', 'none')
    bar(2, mean(propWhiskingLasso(:,i)), 0.3, 'c', 'linestyle', 'none')
    bar(2.3, mean(propWhiskingRidge(:,i)), 0.3, 'm', 'linestyle', 'none')
    
    errorbar(1, mean(propTouchLasso(:,i)), std(propTouchLasso(:,i))/sqrt(length(mice)), 'b', 'lineWidth', 2)
    errorbar(1.3, mean(propTouchRidge(:,i)), std(propTouchRidge(:,i))/sqrt(length(mice)), 'r', 'lineWidth', 2)
    errorbar(2, mean(propWhiskingLasso(:,i)), std(propWhiskingLasso(:,i))/sqrt(length(mice)), 'c', 'lineWidth', 2)
    errorbar(2.3, mean(propWhiskingRidge(:,i)), std(propWhiskingRidge(:,i))/sqrt(length(mice)), 'm', 'lineWidth', 2)
    if mod(i,2)
        ylim([0 0.4])
    else
        ylim([0 0.2])
    end
    xticks([1.15, 2.15])
    xticklabels({'Touch', 'Whisking'})
end

legend({'Lasso', 'Ridge'})

%% same thing from naive mice
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
matchingMi = [1,2,3,4,7,9];
propTouchLasso = zeros(length(mice), 4); % 1 L2/3 C2, 2 L4 C2, 3 L2/3 non-C2, 4 L4 non-C2
propTouchRidge = zeros(length(mice), 4);
propWhiskingLasso = zeros(length(mice), 4);
propWhiskingRidge = zeros(length(mice), 4);

for mi = 1 : length(matchingMi)
    temp = lasso010.naive(matchingMi(mi));
    idL23C2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(temp.isC2)));
    idL4C2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(temp.isC2)));
    idL23nonC2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(1 - temp.isC2)));
    idL4nonC2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(1 - temp.isC2)));
    
    propTouchLasso(mi,1) = length(intersect(temp.touchID, idL23C2)) / length(idL23C2);
    propTouchLasso(mi,2) = length(intersect(temp.touchID, idL4C2)) / length(idL4C2);
    propTouchLasso(mi,3) = length(intersect(temp.touchID, idL23nonC2)) / length(idL23nonC2);
    propTouchLasso(mi,4) = length(intersect(temp.touchID, idL4nonC2)) / length(idL4nonC2);
    
    propWhiskingLasso(mi,1) = length(intersect(temp.whiskingID, idL23C2)) / length(idL23C2);
    propWhiskingLasso(mi,2) = length(intersect(temp.whiskingID, idL4C2)) / length(idL4C2);
    propWhiskingLasso(mi,3) = length(intersect(temp.whiskingID, idL23nonC2)) / length(idL23nonC2);
    propWhiskingLasso(mi,4) = length(intersect(temp.whiskingID, idL4nonC2)) / length(idL4nonC2);
    
    temp = ridge010.naive(matchingMi(mi));
    idL23C2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(temp.isC2)));
    idL4C2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(temp.isC2)));
    idL23nonC2 = temp.cellNums(intersect(find(temp.cellDepths < 350), find(1 - temp.isC2)));
    idL4nonC2 = temp.cellNums(intersect(find(temp.cellDepths >= 350), find(1 - temp.isC2)));
    
    propTouchRidge(mi,1) = length(intersect(temp.touchID, idL23C2)) / length(idL23C2);
    propTouchRidge(mi,2) = length(intersect(temp.touchID, idL4C2)) / length(idL4C2);
    propTouchRidge(mi,3) = length(intersect(temp.touchID, idL23nonC2)) / length(idL23nonC2);
    propTouchRidge(mi,4) = length(intersect(temp.touchID, idL4nonC2)) / length(idL4nonC2);
    
    propWhiskingRidge(mi,1) = length(intersect(temp.whiskingID, idL23C2)) / length(idL23C2);
    propWhiskingRidge(mi,2) = length(intersect(temp.whiskingID, idL4C2)) / length(idL4C2);
    propWhiskingRidge(mi,3) = length(intersect(temp.whiskingID, idL23nonC2)) / length(idL23nonC2);
    propWhiskingRidge(mi,4) = length(intersect(temp.whiskingID, idL4nonC2)) / length(idL4nonC2);
    
    
end

figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(1, mean(propTouchLasso(:,i)), 0.3, 'b', 'linestyle', 'none')
    bar(1.3, mean(propTouchRidge(:,i)), 0.3, 'r', 'linestyle', 'none')
    bar(2, mean(propWhiskingLasso(:,i)), 0.3, 'c', 'linestyle', 'none')
    bar(2.3, mean(propWhiskingRidge(:,i)), 0.3, 'm', 'linestyle', 'none')
    
    errorbar(1, mean(propTouchLasso(:,i)), std(propTouchLasso(:,i))/sqrt(length(mice)), 'b', 'lineWidth', 2)
    errorbar(1.3, mean(propTouchRidge(:,i)), std(propTouchRidge(:,i))/sqrt(length(mice)), 'r', 'lineWidth', 2)
    errorbar(2, mean(propWhiskingLasso(:,i)), std(propWhiskingLasso(:,i))/sqrt(length(mice)), 'c', 'lineWidth', 2)
    errorbar(2.3, mean(propWhiskingRidge(:,i)), std(propWhiskingRidge(:,i))/sqrt(length(mice)), 'm', 'lineWidth', 2)
    if mod(i,2)
        ylim([0 0.4])
    else
        ylim([0 0.2])
    end
    xticks([1.15, 2.15])
    xticklabels({'Touch', 'Whisking'})
end

legend({'Lasso', 'Ridge'})


%% Decided to use Ridge 0.1

%% Proportion of functions
%% From now on, focus on matching naive and experts
% sound, reward, licking.
% how many of union(touch, whisking) have those 3 functions?
matchingMi = [1,2,3,4,7,9];

propSoundNaive = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))
propSoundExpert = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))
propRewardNaive = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))
propRewardExpert = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))
propLickingNaive = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))
propLickingExpert = zeros(length(matchingMi), 4); % 1 / all, 2 / L4, 3 / union(touch, whisking), 4 / intersect(L4, union(touch, whisking))

for mi = 1 : length(matchingMi)
    temp = ridge010.naive(matchingMi(mi));
    idL4 = temp.cellNums(temp.cellDepths >= 350);
    idTW = union(temp.touchID, temp.whiskingID);
    idL4TW = intersect(idL4, idTW);
    propSoundNaive(mi,1) = length(temp.soundID) / length(temp.cellNums);
    propSoundNaive(mi,2) = length(intersect(temp.soundID, idL4)) / length(idL4);
    propSoundNaive(mi,3) = length(intersect(temp.soundID, idTW)) / length(idTW);
    propSoundNaive(mi,4) = length(intersect(temp.soundID, idL4TW)) / length(idL4TW);
    
    propRewardNaive(mi,1) = length(temp.rewardID) / length(temp.cellNums);
    propRewardNaive(mi,2) = length(intersect(temp.rewardID, idL4)) / length(idL4);
    propRewardNaive(mi,3) = length(intersect(temp.rewardID, idTW)) / length(idTW);
    propRewardNaive(mi,4) = length(intersect(temp.rewardID, idL4TW)) / length(idL4TW);
    
    propLickingNaive(mi,1) = length(temp.lickingID) / length(temp.cellNums);
    propLickingNaive(mi,2) = length(intersect(temp.lickingID, idL4)) / length(idL4);
    propLickingNaive(mi,3) = length(intersect(temp.lickingID, idTW)) / length(idTW);
    propLickingNaive(mi,4) = length(intersect(temp.lickingID, idL4TW)) / length(idL4TW);
    
    
    temp = ridge010.expert(mi);
    idL4 = temp.cellNums(temp.cellDepths >= 350);
    idTW = union(temp.touchID, temp.whiskingID);
    idL4TW = intersect(idL4, idTW);
    propSoundExpert(mi,1) = length(temp.soundID) / length(temp.cellNums);
    propSoundExpert(mi,2) = length(intersect(temp.soundID, idL4)) / length(idL4);
    propSoundExpert(mi,3) = length(intersect(temp.soundID, idTW)) / length(idTW);
    propSoundExpert(mi,4) = length(intersect(temp.soundID, idL4TW)) / length(idL4TW);
    
    propRewardExpert(mi,1) = length(temp.rewardID) / length(temp.cellNums);
    propRewardExpert(mi,2) = length(intersect(temp.rewardID, idL4)) / length(idL4);
    propRewardExpert(mi,3) = length(intersect(temp.rewardID, idTW)) / length(idTW);
    propRewardExpert(mi,4) = length(intersect(temp.rewardID, idL4TW)) / length(idL4TW);
    
    propLickingExpert(mi,1) = length(temp.lickingID) / length(temp.cellNums);
    propLickingExpert(mi,2) = length(intersect(temp.lickingID, idL4)) / length(idL4);
    propLickingExpert(mi,3) = length(intersect(temp.lickingID, idTW)) / length(idTW);
    propLickingExpert(mi,4) = length(intersect(temp.lickingID, idL4TW)) / length(idL4TW);
end
%%
figure,
subplot(131), hold all
bar(1:4, mean(propSoundNaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:4.3, mean(propSoundExpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1:4, mean(propSoundNaive), std(propSoundNaive)/sqrt(length(matchingMi)), 'b.', 'linewidth', 2)
errorbar(1.3:4.3, mean(propSoundExpert), std(propSoundExpert)/sqrt(length(matchingMi)), 'r.', 'linewidth', 2)
xticks([1:4]), xticklabels({'All', 'L4', 'Touch & Whisking', 'L4 T & W'}), xtickangle(45)
title('Sound')
ylim([0 0.03])
ylabel('Proportion')

subplot(132), hold all
bar(1:4, mean(propRewardNaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:4.3, mean(propRewardExpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1:4, mean(propRewardNaive), std(propRewardNaive)/sqrt(length(matchingMi)), 'b.', 'linewidth', 2)
errorbar(1.3:4.3, mean(propRewardExpert), std(propRewardExpert)/sqrt(length(matchingMi)), 'r.', 'linewidth', 2)
xticks([1:4]), xticklabels({'All', 'L4', 'Touch & Whisking', 'L4 T & W'}), xtickangle(45)
ylim([0 0.03])
title('Reward')

subplot(133), hold all
bar(1:4, mean(propLickingNaive), 0.3, 'b', 'linestyle', 'none')
bar(1.3:4.3, mean(propLickingExpert), 0.3, 'r', 'linestyle', 'none')
errorbar(1:4, mean(propLickingNaive), std(propLickingNaive)/sqrt(length(matchingMi)), 'b.', 'linewidth', 2)
errorbar(1.3:4.3, mean(propLickingExpert), std(propLickingExpert)/sqrt(length(matchingMi)), 'r.', 'linewidth', 2)
xticks([1:4]), xticklabels({'All', 'L4', 'Touch & Whisking', 'L4 T & W'}), xtickangle(45)
ylim([0 0.03])
title('Licking')
legend('Naive', 'Expert')


    
    
    
    
    
    
%% Radial distance tasks

clear
tic
baseDir = 'D:\JK\suite2p\';

mice = [25,27,30,36,39,52];
sessions = [22,17,22,18,24,26]; 

for i = 1 : length(mice)
    mouse = mice(i);
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions(i);
    distance(i) = glm_results_cell_function(mouse, session, baseDir);
end

save('Y:\Whiskernas\JK\suite2p\cellFunctionRidgeDE010Distance.mat', 'distance')
toc

    
    