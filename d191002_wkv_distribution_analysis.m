%% whisker kinematic variables
%% How do their variable importance distributed?
%% Can we say whether a cell is dKv responsive, vs dPhi responsive?
%% How are they co-distributed?


%% Basic settings
clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';
loadfn = 'glmResults_WKV_touchCell_exclusion_NC';
whisker = load(sprintf('%s%s',baseDir,loadfn));
touch = load(sprintf('%sglmResults_devExp_touch_NC', baseDir));
tuning = load(sprintf('%sangle_tuning_summary_preAnswer_perTouch_NC', baseDir));

featureNames = {'maxDtheta', 'maxDphi', 'maxDkappaH', 'maxDkappaV', 'maxSlideDistance', 'maxDuration', ...    
                            'thetaAtTouch', 'phiAtTouch', 'kappaHAtTouch', 'kappaVAtTouch', 'arcLengthAtTouch', 'touchCount'};

% load features
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
features = cell(length(mice),1);
for i = 1 : length(mice)
    load(sprintf('%s%03d\\glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R01', baseDir, mice(i), mice(i), sessions{i}(1)), 'allPredictors')
    selectInds = 1:3:34;
    features{i} = allPredictors{1}(:,selectInds);
end

%% First of all, is whisker variable as important as touch variable??
touchVi = cell(length(touch.naive),1);
whiskerVi = cell(length(whisker.naive),1);
for i = 1 : length(touch.naive)
    touchInd = find(ismember(touch.naive(i).cellID, whisker.naive(i).cID));
    touchVi{i} = cellfun(@(x) x(1), touch.naive(i).partialSub(touchInd));
    whiskerVi{i} = cellfun(@(x) x(1), whisker.naive(i).DEdiff);
end

%%
figure, hold on
scatter(cell2mat(touchVi), cell2mat(whiskerVi), 'k.')
plot([0 0.8], [0 0.8], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2) 
plot([0 0.8], [0.1 0.1], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2) 
axis equal
xlim([0.1 0.8])
ylim([0 0.7])
xlabel('Touch variable importance')
ylabel('Whisker variable importance')

corr(cell2mat(touchVi), cell2mat(whiskerVi))

%% proportion of failed touch in whisker glm 
propFail = zeros(length(whiskerVi),1);
for i = 1 : length(whiskerVi)
    propFail(i) = length(find(whiskerVi{i} < 0.1)) / length(whiskerVi{i});
end

mean(propFail)
std(propFail)/sqrt(length(propFail))

%% Result: 4.21 +- 0.47 % failed to be touch cell in whisker glm

%% Distribution of each features variable importance
vi = cell(12,12); % rows: each mouse, columns: each features. e.g., dKv (4th feature) of mouse #2: vi{2,4}
for i = 1 : length(whisker.naive)
    for j = 1 : 12
        vi{i,j} = whisker.naive(i).whiskerVariableDEdiff(:,j);
    end
end

%%
binSize = 0.005;
histRange = [0:binSize:0.15, 10];
spInd = [1, 2, 5, 6, 9, 10, 3, 4, 7, 8, 11, 12]; % subplot index
figure, 
for i = 1 : 12
    subplot(3,4,spInd(i))
    temp = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), vi(:,i), 'un', 0));
    boundedline(histRange(1:end-1)+binSize/2, mean(temp), std(temp)/sqrt(length(whisker.naive)), 'cmap', [0 0 0])
    title(featureNames{i})
    ylim([0 0.1])
end

%% Calculate infliction points in all features
binSize = 0.005;
histRange = [0:binSize:0.15, 10];
inflictPoint = zeros(12,1);
for i = 1
    temp = mean(cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), vi(:,i), 'un', 0)));
%     [~, maxi] = max(diff(diff(temp)))
end
%% What if you just take max?
% how do their value look like in distribution

vi = cell(length(whisker.naive),1);
maxval = cell(length(whisker.naive),1);
maxind = cell(length(whisker.naive),1);
for i = 1 : length(whisker.naive)
    vi{i} = whisker.naive(i).whiskerVariableDEdiff;
    [maxval{i}, maxind{i}] = max(vi{i},[],2);
    
end
%%
histBinSize = 0.01;
histRange = 0:histBinSize:0.4;
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxval, 'un', 0));
figure,
boundedline(histRange(1:end-1)+histBinSize/2, mean(tempmat), std(tempmat)/sqrt(size(tempmat,1)))
xlabel('Max VI each cell')
ylabel('Proportion')
ylim([0 0.14])
mean(cellfun(@mean, maxval))
std(cellfun(@mean, maxval))/sqrt(size(maxval,1))

%%
histRange = 1:size(whisker.naive(1).whiskerVariableDEdiff,2)+1;
tempmat = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), maxind, 'un', 0));
xposition = [1:6, 8:13];
figure, hold on
bar(xposition, mean(tempmat), 'k')
errorbar(xposition, mean(tempmat), std(tempmat)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
xticks(xposition)
xtickangle(45)
xticklabels(featureNames)
ylabel('Proportion')


%% How about in not-tuned cells?

viNT = cell(length(whisker.naive),1);
maxvalNT = cell(length(whisker.naive),1);
maxindNT = cell(length(whisker.naive),1);
for i = 1 : length(whisker.naive)
    viNT{i} = whisker.naive(i).whiskerVariableDEdiff(find(tuning.naive(i).tuned==0),:);
    [maxvalNT{i}, maxindNT{i}] = max(viNT{i},[],2);
end

%%
histBinSize = 0.01;
histRange = 0:histBinSize:0.4;
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalNT, 'un', 0));
figure,
boundedline(histRange(1:end-1)+histBinSize/2, mean(tempmat), std(tempmat)/sqrt(size(tempmat,1)), 'k')
xlabel('Max VI each cell')
ylabel('Proportion')
ylim([0 0.14])
title(sprintf('Not-tuned, %.3f \\pm %.3f', mean(cellfun(@mean, maxvalNT)), std(cellfun(@mean, maxvalNT))/sqrt(size(maxvalNT,1))))

histRange = 1:size(whisker.naive(1).whiskerVariableDEdiff,2)+1;
tempmat = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), maxindNT, 'un', 0));
xposition = [1:6, 8:13];
figure, hold on
bar(xposition, mean(tempmat), 'k')
errorbar(xposition, mean(tempmat), std(tempmat)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
xticks(xposition)
xtickangle(45)
xticklabels(featureNames)
ylabel('Proportion')
title('Not-tuned')

%% Compare with tuned cells

viTuned = cell(length(whisker.naive),1);
maxvalTuned = cell(length(whisker.naive),1);
maxindTuned = cell(length(whisker.naive),1);
for i = 1 : length(whisker.naive)
    viTuned{i} = whisker.naive(i).whiskerVariableDEdiff(find(tuning.naive(i).tuned==1),:);
    [maxvalTuned{i}, maxindTuned{i}] = max(viTuned{i},[],2);
end

%
histBinSize = 0.01;
histRange = 0:histBinSize:0.4;
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalTuned, 'un', 0));
figure,
boundedline(histRange(1:end-1)+histBinSize/2, mean(tempmat), std(tempmat)/sqrt(size(tempmat,1)), 'k')
xlabel('Max VI each cell')
ylabel('Proportion')
ylim([0 0.14])
mean(cellfun(@mean, maxvalTuned))
std(cellfun(@mean, maxvalTuned))/sqrt(size(maxvalTuned,1))
title(sprintf('Tuned, %.3f \\pm %.3f', mean(cellfun(@mean, maxvalTuned)), std(cellfun(@mean, maxvalTuned))/sqrt(size(maxvalTuned,1))))

%
histRange = 1:size(whisker.naive(1).whiskerVariableDEdiff,2)+1;
tempmat = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), maxindTuned, 'un', 0));
xposition = [1:6, 8:13];
figure, hold on
bar(xposition, mean(tempmat), 'k')
errorbar(xposition, mean(tempmat), std(tempmat)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
xticks(xposition)
xtickangle(45)
xticklabels(featureNames)
ylabel('Proportion')
title('Tuned')

%% Show them together.

histBinSize = 0.01;
histRange = 0:histBinSize:0.4;
figure, hold on
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalNT, 'un', 0));
plot(histRange(1:end-1)+histBinSize/2, mean(tempmat), 'r')
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalTuned, 'un', 0));
plot(histRange(1:end-1)+histBinSize/2, mean(tempmat), 'k')
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalNT, 'un', 0));
boundedline(histRange(1:end-1)+histBinSize/2, mean(tempmat), std(tempmat)/sqrt(size(tempmat,1)), 'r')
tempmat = cell2mat(cellfun(@(x) histcounts(x,histRange,'normalization','probability'), maxvalTuned, 'un', 0));
boundedline(histRange(1:end-1)+histBinSize/2, mean(tempmat), std(tempmat)/sqrt(size(tempmat,1)), 'k')
ylim([0 0.14])
ylabel('Proportion')
xlabel('Variable importance (max)')
legend({'Not-tuned', 'Tuned'})

%%
histRange = 1:size(whisker.naive(1).whiskerVariableDEdiff,2)+1;
xposition = [1:6, 8:13];


pvalues = zeros(length(featureNames),1);
for i = 1 : length(featureNames)
    [~, pvalues(i)] = ttest(featurematTuned(:,i), featurematNT(:,i));
end
nsInds = find(pvalues > 0.05);

binwidth = 0.4;
positionAdjust = 0.2;
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.5 0.45]), hold on
featurematTuned = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), maxindTuned, 'un', 0));
featurematNT = cell2mat(cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), maxindNT, 'un', 0));
bar(xposition-positionAdjust, mean(featurematTuned), binwidth, 'k')
bar(xposition+positionAdjust, mean(featurematNT), binwidth, 'w')
errorbar(xposition-positionAdjust, mean(featurematTuned), zeros(1,length(xposition)), std(featurematTuned)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
errorbar(xposition+positionAdjust, mean(featurematNT), zeros(1,length(xposition)), std(featurematNT)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')

% nsTuned = nan(1,size(featurematTuned,2)); nsTuned(nsInds) = mean(featurematTuned(:,nsInds));
% nsNT = nan(1,size(featurematNT,2)); nsNT(nsInds) = mean(featurematNT(:,nsInds));
% bar(xposition-positionAdjust, nsTuned, binwidth, 'facecolor', ones(1,3)*0.3)
% bar(xposition+positionAdjust, nsNT, binwidth, 'facecolor', ones(1,3)*0.9)

sInds = setdiff(1:length(featureNames), nsInds);
sTuned = nan(1,size(featurematTuned,2)); sTuned(sInds) = mean(featurematTuned(:,sInds));
sNT = nan(1,size(featurematNT,2)); sNT(sInds) = mean(featurematNT(:,sInds));

bar(xposition-positionAdjust, sTuned, binwidth, 'k', 'edgecolor', 'r')
bar(xposition+positionAdjust, sNT, binwidth, 'w', 'edgecolor', 'r')


xticks(xposition)
xtickangle(45)
xticklabels(featureNames)
ylabel('Proportion')

set(gca, 'fontsize', 12, 'fontname', 'Arial')
legend({'Tuned','Not-tuned'}, 'location', 'northwest')

%%



%% There are some differences, including increase in dKv and slide distance and
%% decrease in dKh and touch counts. The increased ones (3, including KvAT) can be the candidate for angle tuning.

%% What if I compare average vi between not-tuned and tuned cells?

meanTuned = cell2mat(cellfun(@mean, viTuned, 'un', 0));
meanNT = cell2mat(cellfun(@mean, viNT, 'un', 0));

pvalues = zeros(size(meanTuned,2),1);
for i = 1 : size(meanTuned,2)
    [~, pvalues(i)] = ttest(meanTuned(:,i), meanNT(:,i));
end

nsInds = find(pvalues > 0.05);

xposition = [1:6, 8:13];
binwidth = 0.4;
positionAdjust = 0.2;
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.5 0.45]), hold on
bar(xposition-positionAdjust, mean(meanTuned), binwidth, 'k')
bar(xposition+positionAdjust, mean(meanNT), binwidth, 'w')
errorbar(xposition-positionAdjust, mean(meanTuned), zeros(1,length(xposition)), std(meanTuned)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
errorbar(xposition+positionAdjust, mean(meanNT), zeros(1,length(xposition)), std(meanNT)/sqrt(length(whisker.naive)), 'k', 'linestyle', 'none')
% nsTuned = nan(1,size(meanTuned,2)); nsTuned(nsInds) = mean(meanTuned(:,nsInds));
% nsNT = nan(1,size(meanTuned,2)); nsNT(nsInds) = mean(meanNT(:,nsInds));
% bar(xposition-positionAdjust, nsTuned, binwidth, 'facecolor', ones(1,3)*0.2)
% bar(xposition+positionAdjust, nsNT, binwidth, 'facecolor', ones(1,3)*0.8)

sInds = setdiff(1:length(featureNames), nsInds);
sTuned = nan(1,size(meanTuned,2)); sTuned(sInds) = mean(meanTuned(:,sInds));
sNT = nan(1,size(meanNT,2)); sNT(sInds) = mean(meanNT(:,sInds));

bar(xposition-positionAdjust, sTuned, binwidth, 'k', 'edgecolor', 'r')
bar(xposition+positionAdjust, sNT, binwidth, 'w', 'edgecolor', 'r')

xticks(xposition)
xtickangle(45)
xticklabels(featureNames)
ylabel('Averaged VI values')
set(gca, 'fontsize', 12, 'fontname', 'Arial')
legend({'Tuned','Not-tuned'}, 'location', 'northwest')




%% co-activation?
corrWhiskerVI = zeros(length(featureNames), length(featureNames),3);
corrWhiskerFeature = zeros(length(featureNames), length(featureNames),3);
for i = 1 : length(whisker.naive)
    corrWhiskerVI(:,:,i) = corr(whisker.naive(i).whiskerVariableDEdiff);
    corrWhiskerFeature(:,:,i) = corr(features{i}(isfinite(sum(features{i},2)), :));
end
meanCorrVI = mean(corrWhiskerVI,3);
meanCorrFeature = mean(corrWhiskerFeature,3);
%%
figure, 
subplot(121), imagesc(meanCorrVI), colorbar, axis square
title('Variable importance correlation between features')
subplot(122), imagesc(meanCorrFeature), colorbar, axis square
title('Correlation between features')



%% take the ones that are within 1 sd from the peak DE

i = 1;
    onesd = 2 * std(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    threshold = max(whisker.naive(i).whiskerVariableDEdiff, [], 2) - onesd;
    [~, maxInd] = max(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    coactive = cell(length(threshold),1);
    caMat = zeros(length(featureNames), length(featureNames));
    for ci = 1 : length(threshold)
        coactive{ci} = setdiff(find(whisker.naive(i).whiskerVariableDEdiff(ci,:) > threshold(ci)), maxInd(ci));
        caMat(maxInd(ci), coactive{ci}) = caMat(maxInd(ci), coactive{ci}) + 1;
    end
    
    figure, 
    imagesc(caMat)
    xlabel('Sub features'), ylabel('Max feature')
    
    
%% Comparison between 1st vi value and the 2nd vi value
secondIndRatio = zeros(length(mice),1);
figure, hold on
for i = 1 : length(mice)
    [valFirst, maxInd] = max(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    valSecond = zeros(size(valFirst));
    sdThreshold = 1 * std(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    threshold = max(whisker.naive(i).whiskerVariableDEdiff, [], 2) - sdThreshold;
    coactive = cell(length(threshold),1);
    for ci = 1 : length(valFirst)
        coactive{ci} = setdiff(find(whisker.naive(i).whiskerVariableDEdiff(ci,:) > threshold(ci)), maxInd(ci));
        valSecond(ci) = max(whisker.naive(i).whiskerVariableDEdiff(ci,setdiff(1:length(featureNames), maxInd(ci))));
    end
    scatter(valFirst, valSecond, 'k.'), axis equal
    secondInd = find(cellfun(@length, coactive));
    scatter(valFirst(secondInd), valSecond(secondInd), 'r.')
    secondIndRatio(i) = length(secondInd) / length(coactive);
end
xlim([0 0.3]), ylim([0 0.3])
legend({'Fit = 1', 'Fit > 1'}, 'autoupdate', 0, 'location', 'northwest')
plot([0 1], [0 1], '--', 'color', [0.7 0.7 0.7])
xlabel('First value'), ylabel('Second value')

mean(secondIndRatio)
std(secondIndRatio)/sqrt(length(mice))
%%


%%
caMatAll = zeros(length(featureNames), length(featureNames), length(mice));
for i = 1 : 12

    onesd = 2 * std(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    threshold = max(whisker.naive(i).whiskerVariableDEdiff, [], 2) - onesd;
    [~, maxInd] = max(whisker.naive(i).whiskerVariableDEdiff, [], 2);
    coactive = cell(length(threshold),1);
    caMat = zeros(length(featureNames), length(featureNames));
    for ci = 1 : length(threshold)
        coactive{ci} = setdiff(find(whisker.naive(i).whiskerVariableDEdiff(ci,:) > threshold(ci)), maxInd(ci));
        caMat(maxInd(ci), coactive{ci}) = caMat(maxInd(ci), coactive{ci}) + 1;
    end    
    caMatAll(:,:,i) = caMat./length(maxInd);
end

figure,
imagesc(nanmean(caMatAll,3)), axis square, colorbar
xlabel('Sub features'), ylabel('Max feature')



%% coactive features, tuned vs not-tuned

caMatTuned = zeros(length(featureNames), length(featureNames), length(mice));
caMatNT = zeros(length(featureNames), length(featureNames), length(mice));
for i = 1 : 12
    dd = whisker.naive(i).whiskerVariableDEdiff;
    ddTuned = dd(find(tuning.naive(i).tuned), :);
    
    onesd = 2 * std(ddTuned, [], 2);
    threshold = max(ddTuned, [], 2) - onesd;
    [~, maxInd] = max(ddTuned, [], 2);
    coactive = cell(length(threshold),1);
    caMat = zeros(length(featureNames), length(featureNames));
    for ci = 1 : length(threshold)
        coactive{ci} = setdiff(find(ddTuned(ci,:) > threshold(ci)), maxInd(ci));
        caMat(maxInd(ci), coactive{ci}) = caMat(maxInd(ci), coactive{ci}) + 1;
    end
%     occurrence = zeros(length(featureNames),1);
%     for fi = 1 : length(featureNames)
%         occurrence(fi) = length(find(maxInd == fi));
%     end
%     caMatTuned(:,:,i) = caMat./occurrence;
    caMatTuned(:,:,i) = caMat./length(maxInd);
    
    ddNT = dd(find(tuning.naive(i).tuned == 0), :);
    
    onesd = 2 * std(ddNT, [], 2);
    threshold = max(ddNT, [], 2) - onesd;
    [~, maxInd] = max(ddNT, [], 2);
    coactive = cell(length(threshold),1);
    caMat = zeros(length(featureNames), length(featureNames));
    for ci = 1 : length(threshold)
        coactive{ci} = setdiff(find(ddNT(ci,:) > threshold(ci)), maxInd(ci));
        caMat(maxInd(ci), coactive{ci}) = caMat(maxInd(ci), coactive{ci}) + 1;
    end
%     occurrence = zeros(length(featureNames),1);
%     for fi = 1 : length(featureNames)
%         occurrence(fi) = length(find(maxInd == fi));
%     end
%     caMatNT(:,:,i) = caMat./occurrence;
    caMatNT(:,:,i) = caMat./length(maxInd);
end

figure,
subplot(121),
imagesc(nanmean(caMatTuned,3)), axis square, colorbar
xlabel('Sub features'), ylabel('Max feature'), title('Tuned')

subplot(122),
imagesc(nanmean(caMatNT,3)), axis square, colorbar
xlabel('Sub features'), ylabel('Max feature'), title('Not-tuned')
