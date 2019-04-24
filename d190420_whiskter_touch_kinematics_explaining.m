%% The influence of whisker kinematics
%% Compare between drop-out fitting and drop-out prediction.
%% DE diff comparison, DE comparison, ratio of error ratio comparison, etc.
%% Try finding out anything that makes both of the method be correlated with each other.
%% If not, try picking out the ones that are correlated, saying the others are hard to assign.
%% The question is, how many of them can be effectively explained by whisker touch kinematics?
%% Which touch kinematics best describe those cells?

%% first, look at difference in DE diff in each one
clear
baseDir = 'C:\JK\';
cd(baseDir)
fullModel = load('glm_cell_function_error_ratio_withWTV_shuffling', 'naive', 'expert');
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');
touchModel = load('glm_results_responseType', 'naive', 'expert');

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

%% Are the two methods correlated?

figure, 
subplot(121), hold on
for i = 1 : 12

    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) - DEdiffFromFull(:,2), DEdiffFromPartial(:,1) - DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch - wtv)')
ylabel('DE diff from partial models (touch - wtv)')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) - DEdiffFromFull(:,2), DEdiffFromPartial(:,1) - DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch - wtv)')
ylabel('DE diff from partial models (touch - wtv)')
title('Expert (n=6)')

%% Results: both naive and expert have correlated relationship between two

%% ratio of DE diff

figure, 
subplot(121), hold on
% for i = 1 : 12
for i = 1 : 8
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) ./ DEdiffFromFull(:,2), DEdiffFromPartial(:,1) ./ DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch / wtv)')
ylabel('DE diff from partial models (touch / wtv)')
title('All naive (n=8*)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) ./ DEdiffFromFull(:,2), DEdiffFromPartial(:,1) ./ DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch / wtv)')
ylabel('DE diff from partial models (touch / wtv)')
title('Expert (n=6)')


%% Results: There are some negative values. Mostly (large ones) from partial wvt model (touch ONLY model)
% Which means that in these cases adding touch was making the model worse
% There are some (very few) very high values in partial model, which means partial wvt model did not explain any.
% But majority of the data are hard to explain.

%% Compare between touch and wtv in each method
% is there negative correlation? how are the values distributed?
% (1) partial fitting method

figure, 
subplot(121), hold on
for i = 1 : 12
% for i = 1 : 8
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')


%% (2) partial prediction method

figure, 
subplot(121), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')

%% Results:
% They are not negatively correlated in both ways. 
% partial fitting method is more difficult to interpret because of negative values.


%% What are the values in partial prediction method for those of negative values is partial fittin method?
touchNegvalInd = cell(12,1);
wtvNegvalInd = cell(12,1);
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    touchNegvalInd{i} = find(DEdiffFromPartial(:,1) < 0);
    wtvNegvalInd{i} = find(DEdiffFromPartial(:,2) < 0);
end

figure,
subplot(221), hold on
for i = 1 : 12
    
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(touchNegvalInd{i},1), DEdiffFromFull(touchNegvalInd{i},2), 'r.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Naive touch only better than full model')


subplot(222), hold on
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(wtvNegvalInd{i},1), DEdiffFromFull(wtvNegvalInd{i},2), 'b.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Naive WTV only better than full model')
    

% now for expert
touchNegvalInd = cell(6,1);
wtvNegvalInd = cell(6,1);
for i = 1 : 6
    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.expert(i).cellID, cID));

    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    touchNegvalInd{i} = find(DEdiffFromPartial(:,1) < 0);
    wtvNegvalInd{i} = find(DEdiffFromPartial(:,2) < 0);
end


subplot(223), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(touchNegvalInd{i},1), DEdiffFromFull(touchNegvalInd{i},2), 'r.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert touch only better than full model')


subplot(224), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(wtvNegvalInd{i},1), DEdiffFromFull(wtvNegvalInd{i},2), 'b.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert WTV only better than full model')

%% Results: scattered a lot. Ther eis a tendency, but very weak.

%% How does DE diff < 0.1 from partial prediction look like in partial fitting?
% naive
touchLowvalInd = cell(12,1);
wtvLowvalInd = cell(12,1);
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    touchLowvalInd{i} = find(DEdiffFromFull(:,1) < 0.1);
    wtvLowvalInd{i} = find(DEdiffFromFull(:,2) < 0.1);
end

figure,
subplot(221), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(touchLowvalInd{i},1), DEdiffFromPartial(touchLowvalInd{i},2), 'r.')
end

xlabel('DE diff touch from partial fitting')
ylabel('DE diff WTV from partial fitting')
title('Naive without touch similar to full model')

subplot(222), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(wtvLowvalInd{i},1), DEdiffFromPartial(wtvLowvalInd{i},2), 'b.')
end

xlabel('DE diff touch from partial fitting')
ylabel('DE diff WTV from partial fitting')
title('Naive without WTV similar to full model')

%% Results: again, there is a tendency, but very weak and scattered or covers most of the points.
%% Negative values from partial fitting are not necessarily low values in partial prediction (although there is a tendency).
%% It makes it difficult to connect both of them, so just focus on partial prediction method. I have shuffling for it too.
%% It is good enough to know that there is a correlation.

%% Try finding the right threshold from shuffling
% std of each DE diff, calculated back from error ratio (stupid...)
% plot each partial DE diff value to each std (one plot for touch another
% for WTV)
% clear
baseDir = 'C:\JK\';
cd(baseDir)
load('glm_cell_function_error_ratio_withWTV_shuffling', 'naive', 'expert');
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
matchingInd = find(cellfun(@length, sessions)>1);
nonlearnerInd = setdiff(1:length(mice), matchingInd);

touchSTD = cell(12,1);
WTVSTD = cell(12,1);
for i = 1 : 12
    touchDiffShuffling = cell(length(naive(i).devExp),1);
    wtvDiffShuffling = cell(length(naive(i).devExp),1);
    for j = 1 : length(naive(i).devExp)
        touchDiffShuffling{j} = (naive(i).errorRatio{j}(1,:)-1) * naive(i).devExp(j);
        wtvDiffShuffling{j} = (naive(i).errorRatio{j}(2,:)-1) * naive(i).devExp(j);
    end
    
    touchSTD{i} = cellfun(@std, touchDiffShuffling);
    WTVSTD{i} = cellfun(@std, wtvDiffShuffling);
end
figure,
subplot(121), hold on
for i = 1 : 12
    tempDEdiff = cell2mat(naive(i).DEdiff);
    plot(tempDEdiff(:,1), touchSTD{i}, 'k.')
end
xlabel('Touch DE diff')
ylabel('Shuffling DE diff std')

subplot(122), hold on
for i = 1 : 12
    tempDEdiff = cell2mat(naive(i).DEdiff);
    plot(tempDEdiff(:,6), WTVSTD{i}, 'k.')
end
xlabel('WTV DE diff')
ylabel('Shuffling DE diff std')


%% How are they distributed?

figure, histogram(cell2mat(touchSTD)), hold on, histogram(cell2mat(WTVSTD))

%% Results: They are correlated, but 99th percentile of std is about 0.01. 
% 3 times of it is still 0.03, which is going to be about 0.3% of error rate.
% hard to use it as a threshold. 

%% DE diff

figure, 
subplot(121), hold on
for i = 1 : 12
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')

%% Error ratio comparison from shuffling
figure, 
subplot(121), hold on
for i = 1 : 12
    touchER = cellfun(@(x) mean(x(1,:)), naive(i).errorRatio);    
    wtvER = cellfun(@(x) mean(x(2,:)), naive(i).errorRatio);    
    plot(touchER, wtvER, 'k.')
end
xlim([1 2]),ylim([1 2])
xlabel('ER w/o touch')
ylabel('ER w/o WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6
    touchER = cellfun(@(x) mean(x(1,:)), expert(i).errorRatio);    
    wtvER = cellfun(@(x) mean(x(2,:)), expert(i).errorRatio); 
    plot(touchER, wtvER, 'k.')
end
xlim([1 2]),ylim([1 2])
xlabel('ER w/o touch')
ylabel('ER w/o WTV')
title('Expert (n=6)')

%% Error ratio comparison from exclusion
figure, 
subplot(121), hold on
for i = 1 : 12
    erMat = cell2mat(naive(i).exclusionER);
    plot(erMat(:,1), erMat(:,6), 'k.')
end
xlim([1 2]),ylim([1 2])
xlabel('ER w/o touch')
ylabel('ER w/o WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6
    erMat = cell2mat(expert(i).exclusionER);
    plot(erMat(:,1), erMat(:,6), 'k.')
end
xlim([1 2]),ylim([1 2])
xlabel('ER w/o touch')
ylabel('ER w/o WTV')
title('Expert (n=6)')

%% Results: All look similar to each other. It won't matter what to use.

%% Ratio between DE diff
close all
range = logspace(-1,1,100);
DEdiffRatioCdfNL = zeros(6,length(range)-1);
for nli = 1 : 6
    i = nonlearnerInd(nli);
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempRatio = DEdiff(:,2)./DEdiff(:,1);
    tempRatio(tempRatio>10) = deal(9.95);
    DEdiffRatioCdfNL(nli,:) = histcounts(tempRatio,range,'normalization','cdf');
end

DEdiffRatioCdfNaive = zeros(6,length(range)-1);
for mi = 1 : 6
    i = matchingInd(mi);
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempRatio = DEdiff(:,2)./DEdiff(:,1);
    tempRatio(tempRatio>10) = deal(9.95);
    DEdiffRatioCdfNaive(mi,:) = histcounts(tempRatio,range,'normalization','cdf');
end

DEdiffRatioCdfExpert = zeros(6,length(range)-1);
for i = 1 : 6
    fullDEdiff = cell2mat(expert(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempRatio = DEdiff(:,2)./DEdiff(:,1);
    tempRatio(tempRatio>10) = deal(9.95);
    DEdiffRatioCdfExpert(i,:) = histcounts(tempRatio,range,'normalization','cdf');
end
rangePlot = linspace(0,10,100);
a = rangePlot(2:end);
figure, hold on
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)
tempMean = mean(DEdiffRatioCdfNL);
tempSEM = std(DEdiffRatioCdfNL)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'c-')
tempMean = mean(DEdiffRatioCdfNaive);
tempSEM = std(DEdiffRatioCdfNaive)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'b-')
tempMean = mean(DEdiffRatioCdfExpert);
tempSEM = std(DEdiffRatioCdfExpert)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'r-')
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)

xlabel('DE diff WTV/touch')
xticks([0:10])
xticklabels({'10^-^1','10^-^0^.^8','10^-^0^.^6','10^-^0^.^4','10^-^0^.^2','10^0','10^0^.^2','10^0^.^4','10^0^.^6','10^0^.^8','10^1^.^0'})
ylabel('Cumulative proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
ylim([0 1])

%% ratio or difference? 
%% 3d plot 

figure, 
subplot(121), hold on
for i = 1 : 12
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot3(DEdiffFromFull(:,1), DEdiffFromFull(:,2), DEdiffFromFull(:,2) - DEdiffFromFull(:,1), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot3(DEdiffFromFull(:,1), DEdiffFromFull(:,2), DEdiffFromFull(:,2) - DEdiffFromFull(:,1), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')

%%
figure, 
subplot(121), hold on
for i = 1 : 12
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot3(DEdiffFromFull(:,1), DEdiffFromFull(:,2), DEdiffFromFull(:,2)./DEdiffFromFull(:,1), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot3(DEdiffFromFull(:,1), DEdiffFromFull(:,2), DEdiffFromFull(:,2)./DEdiffFromFull(:,1), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')


%% Diff between DE diff
close all
range = linspace(-0.3,0.3,100);
DEdiffRatioCdfNL = zeros(6,length(range)-1);
for nli = 1 : 6
    i = nonlearnerInd(nli);
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempDiff = DEdiff(:,2)-DEdiff(:,1);
    tempDiff(tempDiff>0.3) = deal(0.299);
    DEdiffRatioCdfNL(nli,:) = histcounts(tempDiff,range,'normalization','cdf');
end

DEdiffRatioCdfNaive = zeros(6,length(range)-1);
for mi = 1 : 6
    i = matchingInd(mi);
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempDiff = DEdiff(:,2)-DEdiff(:,1);
    tempDiff(tempDiff>0.3) = deal(0.299);
    DEdiffRatioCdfNaive(mi,:) = histcounts(tempDiff,range,'normalization','cdf');
end

DEdiffRatioCdfExpert = zeros(6,length(range)-1);
for i = 1 : 6
    fullDEdiff = cell2mat(expert(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempDiff = DEdiff(:,2)-DEdiff(:,1);
    tempDiff(tempDiff>0.3) = deal(0.299);
    DEdiffRatioCdfExpert(i,:) = histcounts(tempDiff,range,'normalization','cdf');
end
rangePlot = linspace(-0.3,0.3,100);
a = rangePlot(2:end);
figure, hold on
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)
tempMean = mean(DEdiffRatioCdfNL);
tempSEM = std(DEdiffRatioCdfNL)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'c-')
tempMean = mean(DEdiffRatioCdfNaive);
tempSEM = std(DEdiffRatioCdfNaive)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'b-')
tempMean = mean(DEdiffRatioCdfExpert);
tempSEM = std(DEdiffRatioCdfExpert)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'r-')
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)

xlabel('DE diff WTV - touch')
ylabel('Cumulative proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
ylim([0 1])

%% What are the top WTV features?
wtvNaive = zeros(12,13);
for i = 1 : 12
    wtvNaive(i,:) = mean(naive(i).whiskerVariableDEdiff);
end
figure, boundedline(1:13, mean(wtvNaive), std(wtvNaive)/sqrt(12))
xticks([1:13])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
                        
%% Among the ones depending more on WTV, what are the top features?

wtvNaive = zeros(12,13);
for i = 1 : 12
    fullDEdiff = cell2mat(naive(i).DEdiff);
    DEdiff = fullDEdiff(:,[1,6]);
    tempDiff = DEdiff(:,2)-DEdiff(:,1);
    tempInd = find(tempDiff>=0);
    wtvNaive(i,:) = mean(naive(i).whiskerVariableDEdiff(tempInd,:));
end
figure, boundedline(1:13, mean(wtvNaive), std(wtvNaive)/sqrt(12))
xticks([1:13])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
%% Results: seems like most of them are equally contributing. touch count is the most impactful feature, followed by absDKappaV, absDPhi, DKappaH, DkappaV, Dtheta, slide distance, kappaHatTouch, and kappaVatTouch


%% What about wtv only fitting model?

cd(baseDir)
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');


wtvNaive = zeros(12,13);
for i = 1 : 12
    wtvNaive(i,:) = mean(wtvModel.naive(i).whiskerVariableDEdiff);
end
figure, boundedline(1:13, mean(wtvNaive), std(wtvNaive)/sqrt(12))
xticks([1:13])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
%% Results: still similar, but DkappaV stands out (2nd place), followed by phi at touch and kappaV at touch

%% Look at each individual cells
i = 1; % index of naive
j = 1; % index of cell
% figure, plot(wtvModel


%% sharpness vs DE of each
tuning = load('angle_tuning_summary');
figure, plot(tuning.naive(i).sharpness, wtvModel.naive(i).whiskerVariableDEdiff(:,13), 'k.')

%% Diff between DE diff in partial fitting
close all
range = linspace(-0.2,0.2,100);
DEdiffRatioCdfNL = zeros(6,length(range)-1);
for nli = 1 : 6
    i = nonlearnerInd(nli);
    cID = wtvModel.naive(i).cID; % touch cells    
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    tempDiff = DEdiffFromPartial(:,2)-DEdiffFromPartial(:,1);
    tempDiff(tempDiff>max(range)) = deal(max(range) - mean(diff(range))/2);
    tempDiff(tempDiff<min(range)) = deal(min(range) + mean(diff(range))/2);
    DEdiffRatioCdfNL(nli,:) = histcounts(tempDiff,range,'normalization','cdf');
end

DEdiffRatioCdfNaive = zeros(6,length(range)-1);
for mi = 1 : 6
    i = matchingInd(mi);
    cID = wtvModel.naive(i).cID; % touch cells    
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    tempDiff = DEdiffFromPartial(:,2)-DEdiffFromPartial(:,1);
    tempDiff(tempDiff>max(range)) = deal(max(range) - mean(diff(range))/2);
    tempDiff(tempDiff<min(range)) = deal(min(range) + mean(diff(range))/2);
    DEdiffRatioCdfNaive(mi,:) = histcounts(tempDiff,range,'normalization','cdf');
end

DEdiffRatioCdfExpert = zeros(6,length(range)-1);
for i = 1 : 6
    cID = wtvModel.expert(i).cID; % touch cells    
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    tempDiff = DEdiffFromPartial(:,2)-DEdiffFromPartial(:,1);
    tempDiff(tempDiff>max(range)) = deal(max(range) - mean(diff(range))/2);
    tempDiff(tempDiff<min(range)) = deal(min(range) + mean(diff(range))/2);
    DEdiffRatioCdfExpert(i,:) = histcounts(tempDiff,range,'normalization','cdf');
end
a = range(2:end);
figure, hold on
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)
tempMean = mean(DEdiffRatioCdfNL);
tempSEM = std(DEdiffRatioCdfNL)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'c-')
tempMean = mean(DEdiffRatioCdfNaive);
tempSEM = std(DEdiffRatioCdfNaive)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'b-')
tempMean = mean(DEdiffRatioCdfExpert);
tempSEM = std(DEdiffRatioCdfExpert)/sqrt(6);
boundedline(a, tempMean, tempSEM, 'r-')
plot(a, mean(DEdiffRatioCdfNL), 'c-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfNaive), 'b-', 'linewidth',2)
plot(a, mean(DEdiffRatioCdfExpert), 'r-', 'linewidth',2)

xlabel('DE diff WTV - touch')
ylabel('Cumulative proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
ylim([0 1])




%% Correlation between each wtv and angles
wtvAngleCorr = zeros(12,13);
for mi = 1:12
    mouse = mice(mi);
    session = sessions{mi}(1);
    cd(sprintf('%s%03d',baseDir,mouse))
    load(sprintf('glmWithWhiskerTouchVariables_JK%03dS%02d_R01',mouse,session), 'allPredictors', 'indPartial')
    % load(sprintf('UberJK%03dS%02d',mouse,session))
    angles = 45:15:135;
    angleInds = cell(1,length(angles));
    for i = 1 : length(angles)
        angleInds{i} = find(allPredictors{1}(:,indPartial{1}((i-1)*3+1)) > nanmean(allPredictors{1}(:,indPartial{1}((i-1)*3+1))));
    end
    wtvColumns = zeros(sum(cellfun(@length, angleInds)),13);
    angleColumn = zeros(sum(cellfun(@length, angleInds)),1);
    tempInds1 = [1, cumsum(cellfun(@length, angleInds(1:end-1)))+1];
    tempInds2 = cumsum(cellfun(@length, angleInds));
    for i = 1 : length(angles)
        angleColumn(tempInds1(i):tempInds2(i)) = angles(i);
    end
    for i = 1 : size(wtvColumns,2)
        for j = 1 : length(angles)
            wtvColumns(tempInds1(j):tempInds2(j),i) = allPredictors{1}(angleInds{j}, indPartial{6}((i-1)*3+1));
        end
    end
    
    for i = 1 : 13
        noNaNi = intersect(find(isfinite(angleColumn)), find(isfinite(wtvColumns(:,i))));
        wtvAngleCorr(mi,i) = corr(angleColumn(noNaNi), wtvColumns(noNaNi,i));
    end
end

figure,
% boundedline(1:13, mean(wtvAngleCorr), std(wtvAngleCorr)/sqrt(12))
% hold on
boundedline(1:13, mean(abs(wtvAngleCorr)), std(abs(wtvAngleCorr))/sqrt(12))
xticks([1:13])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
ylabel('Abs corr with angle')


%%
