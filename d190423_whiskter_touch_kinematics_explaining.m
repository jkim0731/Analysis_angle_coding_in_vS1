%% So far, from d190420,
%% Using either full model (partial prediction) or partial model works, since the difference is correlated
% To enhance the effect of wtv (or to remove any interaction between wtv and touch) it is going to be better
% to analyze in partial model
% It probably depends on the question. I want to know if there is "abstract" angle tuned cells, regardless of the within-angle variance.
% Angle variable has no within-angle variance, while having between-angle variance. WTV has both. Therefore, when put together, between-angle variance will be distributed in two groups of variables.
% It is going to be better to have separate fittings to see how much within-angle variance is required for spike prediction.
%% Divide into during touch and instantaneous touch features, and which ones explain more of the tuned cells
% (what about in touch but not tuned cells?)
% In each features, see which type of tuning can be explained the best.
%% Divide into whisker variable-related or touch identity-related groups
%% see how they are distributed: depth, learning, in & out of C2, 45 or 135, single, broad, complex, etc.


%% 
clear
baseDir = 'C:\JK\';
cd(baseDir)
fullModel = load('glm_cell_function_error_ratio_withWTV_shuffling', 'naive', 'expert');
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');
touchModel = load('glm_results_responseType', 'naive', 'expert');

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
matchingInd = find(cellfun(@length, sessions)>1);
nonlearnerInd = setdiff(1:length(mice), matchingInd);
%% Diff between DE diff
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

xlabel('DE(WTV) - DE(touch)')
ylabel('Cumulative proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
ylim([0 1])


%% What are the most important whisker features in predicting spikes?
%% in different groups
cd(baseDir)
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');

wtvNL = zeros(6,13);
for nli = 1 : 6
    i = nonlearnerInd(nli);
    wtvNL(nli,:) = mean(wtvModel.naive(i).whiskerVariableDEdiff);    
end
wtvNaive = zeros(6,13);
for ni = 1 : 6
    i = matchingInd(ni);
    wtvNaive(ni,:) = mean(wtvModel.naive(i).whiskerVariableDEdiff);    
end
wtvExpert = zeros(6,13);
for i = 1 : 6    
    wtvExpert(i,:) = mean(wtvModel.expert(i).whiskerVariableDEdiff);    
end

figure, 
bar([1:8]-0.2, mean(wtvNL(:,1:8)), 0.2, 'c', 'edgecolor', 'c'), hold on
bar(1:8, mean(wtvNaive(:,1:8)), 0.2, 'b', 'edgecolor', 'b'), hold on
bar([1:8]+0.2, mean(wtvExpert(:,1:8)), 0.2, 'r', 'edgecolor', 'r'), hold on
errorbar([1:8]-0.2, mean(wtvNL(:,1:8)), std(wtvNL(:,1:8))/sqrt(6), 'c.', 'linewidth', 2)
errorbar(1:8, mean(wtvNaive(:,1:8)), std(wtvNaive(:,1:8))/sqrt(6), 'b.', 'linewidth', 2)
errorbar([1:8]+0.2, mean(wtvExpert(:,1:8)), std(wtvExpert(:,1:8))/sqrt(6), 'r.', 'linewidth', 2)

bar([10:14]-0.2, mean(wtvNL(:,9:13)), 0.2, 'c', 'edgecolor', 'c'), hold on
bar(10:14, mean(wtvNaive(:,9:13)), 0.2, 'b', 'edgecolor', 'b'), hold on
bar([10:14]+0.2, mean(wtvExpert(:,9:13)), 0.2, 'r', 'edgecolor', 'r'), hold on
errorbar([10:14]-0.2, mean(wtvNL(:,9:13)), std(wtvNL(:,9:13))/sqrt(6), 'c.', 'linewidth', 2)
errorbar(10:14, mean(wtvNaive(:,9:13)), std(wtvNaive(:,9:13))/sqrt(6), 'b.', 'linewidth', 2)
errorbar([10:14]+0.2, mean(wtvExpert(:,9:13)), std(wtvExpert(:,9:13))/sqrt(6), 'r.', 'linewidth', 2)

xticks([1:8,10:14])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
ylabel('DE diff')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
set(gca, 'box', 'off')


%% all naive
cd(baseDir)
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');

wtvNaive = zeros(12,13);
for i = 1 : 12
    wtvNaive(i,:) = mean(wtvModel.naive(i).whiskerVariableDEdiff);
end
figure, 
bar(1:8, mean(wtvNaive(:,1:8))), hold on
errorbar(1:8, mean(wtvNaive(:,1:8)), std(wtvNaive(:,1:8))/sqrt(12), 'k.')
bar(10:14, mean(wtvNaive(:,9:13)))
errorbar(10:14, mean(wtvNaive(:,9:13)), std(wtvNaive(:,9:13))/sqrt(12), 'k.')
xticks([1:8,10:14])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
ylabel('DE diff')
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
set(gca, 'box', 'off')

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
bar(1:8, mean(abs(wtvAngleCorr(:,1:8)))), hold on
errorbar(1:8, mean(abs(wtvAngleCorr(:,1:8))), std(abs(wtvAngleCorr(:,1:8)))/sqrt(12), 'k.')
bar(10:14, mean(abs(wtvAngleCorr(:,9:13))))
errorbar(10:14, mean(abs(wtvAngleCorr(:,9:13))), std(abs(wtvAngleCorr(:,9:13)))/sqrt(12), 'k.')
xticks([1:8,10:14])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
ylabel('Abs correlation with object angle')
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
set(gca, 'box', 'off')






%% Correlation between DE diff of different wtv

ddCorr = zeros(13,13,12);
for mi = 1 : 12
    ddCorr(:,:,mi) = corr(wtvModel.naive(mi).whiskerVariableDEdiff);
end
figure, imagesc(mean(ddCorr,3)), axis off

%% expert
ddCorr = zeros(13,13,6);
for mi = 1 : 6
    ddCorr(:,:,mi) = corr(wtvModel.expert(mi).whiskerVariableDEdiff);
end
figure, imagesc(mean(ddCorr,3)), axis off

%% Results: There are two different groups of cells. Ones that fit to during touch features, and others that fit to at touch features.

%% How are these variables correlated
wtvBetweenCorr = zeros(13,13,12);
for mi = 1:12
    mouse = mice(mi);
    session = sessions{mi}(1);
    cd(sprintf('%s%03d',baseDir,mouse))
    load(sprintf('glmWithWhiskerTouchVariables_JK%03dS%02d_R01',mouse,session), 'allPredictors', 'indPartial')
    noNaNi = find(isfinite(sum(allPredictors{1}(:,indPartial{6}([1:3:37])),2)));
    wtvBetweenCorr(:,:,mi) = corr(allPredictors{1}(noNaNi, indPartial{6}([1:3:37])));
end
%
figure, imagesc(mean(wtvBetweenCorr,3))
axis off
%% How are these different groups of fitting distributed?
% 1. layer
% 2. C2 and non-C2
% 3. tuned and not-tuned
% 4. tuning type (single, broad, complex)
% 5. tuning sharpness
%% first of all, how can you define cells listening to during and at touch?
% mean DE diff?
range = -0.04:0.005:0.03;
ddDist = zeros(12,length(range)-1);
for mi = 1 : 12
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2);
    ddDist(mi,:) = histcounts(ddVal, range, 'normalization', 'cdf');
end
figure, boundedline(range(2:end), mean(ddDist), std(ddDist)/sqrt(12))

%%
range = -0.04:0.005:0.03;
ddDistNL = zeros(6,length(range)-1);
ddDistNaive = zeros(6,length(range)-1);
ddDistExpert = zeros(6,length(range)-1);
for i = 1 : 6
    mi = nonlearnerInd(i);
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2);
    ddDistNL(i,:) = histcounts(ddVal, range, 'normalization', 'cdf');
    mi = matchingInd(i);
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2);
    ddDistNaive(i,:) = histcounts(ddVal, range, 'normalization', 'cdf');
    mi = i;
    ddVal = mean(wtvModel.expert(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.expert(mi).whiskerVariableDEdiff(:,9:13),2);
    ddDistExpert(i,:) = histcounts(ddVal, range, 'normalization', 'cdf');    
end
figure, 
plot(range(2:end), mean(ddDistNL), 'c-'), hold on
plot(range(2:end), mean(ddDistNaive), 'b-')
plot(range(2:end), mean(ddDistExpert), 'r-')
boundedline(range(2:end), mean(ddDistNL), std(ddDistNL)/sqrt(6), 'c-')
boundedline(range(2:end), mean(ddDistNaive), std(ddDistNaive)/sqrt(6), 'b-')
boundedline(range(2:end), mean(ddDistExpert), std(ddDistExpert)/sqrt(6), 'r-')
plot(range(2:end), mean(ddDistNL), 'c-')
plot(range(2:end), mean(ddDistNaive), 'b-')
plot(range(2:end), mean(ddDistExpert), 'r-')
xlabel('meanDE(during touch WKV) - meanDE(at touch WKV)')
ylabel('Cumulative proportion')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
set(gca, 'box', 'off')

%% confirm with distribution 
corrVal = zeros(12,1);
figure, hold on
for i = 1 : 12
    plot(mean(wtvModel.naive(i).whiskerVariableDEdiff(:,1:8),2), mean(wtvModel.naive(i).whiskerVariableDEdiff(:,9:13),2), 'k.')
    corrVal(i) = corr(mean(wtvModel.naive(i).whiskerVariableDEdiff(:,1:8),2), mean(wtvModel.naive(i).whiskerVariableDEdiff(:,9:13),2));
end
xlabel('DE diff during touch WKV')
ylabel('DE diff at touch WKV')

mean(corrVal)
%% Results: About 65 % are listening more to at touch features, while 35 % are to during features.
%% Define cells with mean DE diff < 1 std as at touch cells, and > 1 std as during touch cells.

mi = 1;
ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2); 
threshold = std(ddVal);
atInd = find(ddVal < -threshold);
dtInd = find(ddVal > threshold);

length(atInd)
length(dtInd)

%% Result: too small # of cells survive through this thresholding

%% Look at distribution in different classes (mean +/- sem)

%% L2/3 vs L4

cd(baseDir)
info = load('cellFunctionRidgeDE010', 'naive', 'expert');
L23val = zeros(12,1);
L4val = zeros(12,1);
for mi = 1:12
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2); 
    L23ind = find(info.naive(mi).cellDepths(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))) < 350);
    L4ind = find(info.naive(mi).cellDepths(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))) >= 350);
    L23val(mi) = mean(ddVal(L23ind));
    L4val(mi) = mean(ddVal(L4ind));
end
figure, 
bar(1:2, [mean(L23val), mean(L4val)]), hold on
errorbar(1:2, [mean(L23val), mean(L4val)], [std(L23val)/sqrt(length(L23val)), std(L4val)/sqrt(length(L4val))], 'k.')
ttest(L23val, L4val)

%% slight tendency?
%% look at distribution (DE diff vs depth)
figure, hold on
for mi = 1 : 12
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2); 
    depth = info.naive(mi).cellDepths(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID)));
    plot(ddVal,depth, 'k.')
end
set(gca, 'YDir', 'reverse')
xlabel('DE (during touch) - DE (at touch)')
ylabel('Depth')

%% C2 vs non-C2, in L2/3 or L4
cd(baseDir)
info = load('cellFunctionRidgeDE010', 'naive', 'expert');
L23C2val = zeros(12,1);
L23nonC2val = zeros(12,1);
L4C2val = zeros(12,1);
L4nonC2val = zeros(12,1);
for mi = 1:12
    ddVal = mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,1:8),2) - mean(wtvModel.naive(mi).whiskerVariableDEdiff(:,9:13),2); 
    L23ind = find(info.naive(mi).cellDepths(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))) < 350);
    L4ind = find(info.naive(mi).cellDepths(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))) >= 350);
    C2ind = find(info.naive(mi).isC2(find(ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))));
    nonC2ind = find(info.naive(mi).isC2(find(1-ismember(info.naive(mi).cellNums, wtvModel.naive(mi).cID))));
    L23C2val(mi) = mean(ddVal(intersect(L23ind, C2ind)));
    L23nonC2val(mi) = mean(ddVal(intersect(L23ind, nonC2ind)));
    L4C2val(mi) = mean(ddVal(intersect(L4ind, C2ind)));
    L4nonC2val(mi) = mean(ddVal(intersect(L4ind, nonC2ind)));    
end
figure,
subplot(221), bar(mean(L23C2val)), hold on, errorbar(mean(L23C2val), std(L23C2val)/sqrt(12)), title('L23 C2')
subplot(222), bar(mean(L23nonC2val)), hold on, errorbar(mean(L23nonC2val), std(L23nonC2val)/sqrt(12)), title('L23 non-C2')
subplot(223), bar(mean(L4C2val)), hold on, errorbar(mean(L4C2val), std(L4C2val)/sqrt(12)), title('L4 C2')
subplot(224), bar(nanmean(L4nonC2val)), hold on, errorbar(nanmean(L4nonC2val), nanstd(L4nonC2val)/sqrt(12)), title('L4 non-C2')


%% how can we define wkv listening cells?
% using std?
range = linspace(-0.2,0.2,100);
DEdiffRatioCdf = zeros(12,length(range)-1);
propWKV = zeros(12,1);
propAngle = zeros(12,1);
thresholds = zeros(12,1);
for i = 1 : 12    
    cID = wtvModel.naive(i).cID; % touch cells    
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    tempDiff = DEdiffFromPartial(:,2)-DEdiffFromPartial(:,1);
    threshold = std(tempDiff);
    propWKV(i) = length(find(tempDiff < -threshold))/length(tempDiff);
    propAngle(i) = length(find(tempDiff > threshold))/length(tempDiff);
    thresholds(i) = threshold;
end

mean(propWKV)
std(propWKV)/sqrt(12)

mean(propAngle)
std(propAngle)/sqrt(12)

mean(thresholds)
std(thresholds)/sqrt(12)

%% What about all naive among wkv cells? (listening more to wkv instead of angle identity)
cd(baseDir)
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');

wtvNaive = zeros(12,13);
for i = 1 : 12
    wtvNaive(i,:) = mean(wtvModel.naive(i).whiskerVariableDEdiff);
end
figure, 
bar(1:8, mean(wtvNaive(:,1:8))), hold on
errorbar(1:8, mean(wtvNaive(:,1:8)), std(wtvNaive(:,1:8))/sqrt(12), 'k.')
bar(10:14, mean(wtvNaive(:,9:13)))
errorbar(10:14, mean(wtvNaive(:,9:13)), std(wtvNaive(:,9:13))/sqrt(12), 'k.')
xticks([1:8,10:14])
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
ylabel('DE diff')
xticklabels({'max\Delta\kappa_H', 'max\Delta\kappa_V', 'max\Delta\theta', 'max\Delta\phi', 'maxSlideDistance', 'maxDuration', 'max|\Delta\kappa_V|', 'max|\Delta\phi|', '\thetaAtTouch', '\phiAtTouch', '\kappa_HAtTouch', '\kappa_VAtTouch', 'touchCount'})
xtickangle(45)
set(gca, 'box', 'off')

