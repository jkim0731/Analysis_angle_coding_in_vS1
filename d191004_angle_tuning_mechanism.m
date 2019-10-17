%% For angle tuning mechanism
% 1. See how well the angle tuning is preserved in the full models
%     - tuned vs not
%     - tuned angle
%     - tuning modulation
% 2. Make figures explaining the method 
%     - Raw tuning, model tuning, remove-one tuning, tuning with single feature
% 3. Make figures about the results



%% Basic settings
baseDir = 'Y:\Whiskernas\JK\suite2p\';
loadFn = 'modelAngleTuning_NC';
tuning = load(sprintf('%s%s',baseDir, loadFn), 'naive'); % only with naive, for now
numMice = length(tuning.naive);

%% 1. See how well the angle tuning is preserved in the full models
% inds: 1 raw, 2 touch model, 3 whisker model
%     - tuned vs not
% simiarity is defined by sum(innerproduct)/length
similarity = zeros(numMice,2); % 1: raw vs touch model , 2: raw vs whisker model
for mi = 1 : numMice
    tempMat = tuning.naive(mi).tunedAllCell;
    similarity(mi,1) = length(find(tempMat(:,1) - tempMat(:,2) == 0)) / size(tempMat,1);
    similarity(mi,2) = length(find(tempMat(:,1) - tempMat(:,3) == 0)) / size(tempMat,1);
end

figure, hold on
bar(mean(similarity), 'k')
errorbar(mean(similarity), std(similarity)/sqrt(numMice), 'k', 'linestyle', 'none')
xticks([1 2]), xticklabels({'VS touch model', 'VS whisker model'}), xtickangle(45)
ylabel('Similarity')

mean(similarity)

%% Result: 80 % of touch model, 79% of whisker model. Which ones are wrong?
%(1) are there any "new tuning" from not-tuned cells?
NT2T = zeros(numMice,2); % 1: raw vs touch model , 2: raw vs whisker model
T2NT = zeros(numMice,2);
for mi = 1 : numMice
    tempMat = tuning.naive(mi).tunedAllCell;
    NT2T(mi,1) = sum(tempMat(:,1) < tempMat(:,2)) / size(tempMat,1);
    NT2T(mi,2) = sum(tempMat(:,1) < tempMat(:,3)) / size(tempMat,1);
    
    T2NT(mi,1) = sum(tempMat(:,1) > tempMat(:,2)) / size(tempMat,1);
    T2NT(mi,2) = sum(tempMat(:,1) > tempMat(:,3)) / size(tempMat,1);
end

figure, 
subplot(121), hold on
bar(mean(NT2T), 'k')
errorbar(mean(NT2T), std(NT2T)/sqrt(numMice), 'k' ,'linestyle', 'none')
xticks([1 2]), xticklabels({'VS touch model', 'VS whisker model'}), xtickangle(45)
ylabel('Not-tuned -> Tuned'), ylim([0 0.25])
subplot(122), hold on
bar(mean(T2NT), 'k')
errorbar(mean(T2NT), std(T2NT)/sqrt(numMice), 'k' ,'linestyle', 'none')
xticks([1 2]), xticklabels({'VS touch model', 'VS whisker model'}), xtickangle(45)
ylabel('Tuned -> Not-tuned'), ylim([0 0.25])

%% Result: Most of mismatch from touch model is increase tuning, since the model imposes angle tuning already
%% Mismatch of whisker model is divided roughly in half. 
%% Among tuned, 1% of them are gone from touch model, 9% are gone in whisker model.

%% Which ones are gone in whisker model? What are their types, tuned angle, and modulation level?

mmType = zeros(numMice,3);
mmAngle = zeros(numMice,7);
mmModulation = zeros(numMice,1);
allType = zeros(numMice,3); % 1 sharp, 2 broad, 3 complex
allAngle = zeros(numMice,7);
allModulation = zeros(numMice,1);

for i = 1 : numMice
    numCell = size(tuning.naive(i).anovaPAllCell,1);
    numTuned = sum(tuning.naive(i).tunedAllCell(:,1));
    tunedInd = find(tuning.naive(i).tunedAllCell(:,1));
    allType(i,1) = sum(tuning.naive(i).unimodalSingleAllCell(tunedInd,1)) / numTuned;
    allType(i,2) = sum(tuning.naive(i).unimodalBroadAllCell(tunedInd,1)) / numTuned;
    allType(i,3) = sum(tuning.naive(i).multimodalAllCell(tunedInd,1)) / numTuned;
    angles = 45:15:135;
    for ai = 1 : length(angles)
        allAngle(i,ai) = length(find(tuning.naive(i).tuneAngleAllCell(:,1) == angles(ai))) / numTuned;
    end
    allModulation(i) = mean(tuning.naive(i).tuneModulationAllCell(tunedInd,1));
    
    mmIndTemp = find(tuning.naive(i).tunedAllCell(tunedInd,3) == 0);
    mmInd = tunedInd(mmIndTemp);
    mmType(i,1) = sum(tuning.naive(i).unimodalSingleAllCell(mmInd,1)) / length(mmInd);
    mmType(i,2) = sum(tuning.naive(i).unimodalBroadAllCell(mmInd,1)) / length(mmInd);
    mmType(i,3) = sum(tuning.naive(i).multimodalAllCell(mmInd,1)) / length(mmInd);
    for ai = 1 : length(angles)
        mmAngle(i,ai) = length(find(tuning.naive(i).tuneAngleAllCell(mmInd, 1) == angles(ai))) / length(mmInd);
    end
    mmModulation(i) = mean(tuning.naive(i).tuneModulationAllCell(mmInd,1));
end


figure,
subplot(131), hold on
xposition = 1:3;
positionAdjustment = 0.2;
barWidth = 0.4;
bar(xposition-positionAdjustment, mean(allType), barWidth, 'k')
bar(xposition+positionAdjustment, mean(mmType), barWidth, 'r')
errorbar(xposition-positionAdjustment, mean(allType), std(allType)/sqrt(numMice), 'k', 'linestyle', 'none')
errorbar(xposition+positionAdjustment, mean(mmType), std(mmType)/sqrt(numMice), 'r', 'linestyle', 'none')
xticks(xposition)
xticklabels({'Sharp', 'Broad', 'Complex'})
ylabel('Proportion')
legend({'Spikes', 'Gone in whisker model'}, 'location', 'northwest')

subplot(132), hold on
xposition = 1:7;
positionAdjustment = 0.2;
barWidth = 0.4;
bar(xposition-positionAdjustment, mean(allAngle), barWidth, 'k')
bar(xposition+positionAdjustment, mean(mmAngle), barWidth, 'r')
errorbar(xposition-positionAdjustment, mean(allAngle), std(allAngle)/sqrt(numMice), 'k', 'linestyle', 'none')
errorbar(xposition+positionAdjustment, mean(mmAngle), std(mmAngle)/sqrt(numMice), 'r', 'linestyle', 'none')
xticks(xposition)
xticklabels(angles)
ylabel('Proportion')
legend({'Spikes', 'Gone in whisker model'}, 'location', 'northwest')

subplot(133), hold on
bar(1, mean(allModulation), 'k')
bar(2,mean(mmModulation), 'r')
errorbar(1, mean(allModulation), std(allModulation)/sqrt(numMice), 'k', 'linestyle', 'none')
errorbar(2, mean(mmModulation), std(mmModulation)/sqrt(numMice), 'r', 'linestyle', 'none')
xticks([1 2])
ylabel('Modulation')

%% Results: broad, tuned in the middle angles are more prone to disappear. 

%% How much of each tuned angles are gone? (e.g., are ALL of 90d tuned cells are gone? How many of them are left?)
% (1) proportion per raw tuning, and (2) overall angle distribution from whisker model
angles = 45:15:135;
propPerRaw = zeros(numMice, length(angles));
propWhisker = zeros(numMice, length(angles));
propLeft = zeros(numMice, length(angles));
propRaw = zeros(numMice, length(angles));
for i = 1 : numMice    
    tunedInd = find(tuning.naive(i).tunedAllCell(:,1));
    mmIndTemp = find(tuning.naive(i).tunedAllCell(tunedInd,3) == 0);
    mmInd = tunedInd(mmIndTemp);
    leftIndTemp = find(tuning.naive(i).tunedAllCell(tunedInd,3));
    leftInd = tunedInd(leftIndTemp);
    for ai = 1 : length(angles)        
        propPerRaw(i,ai) = sum(tuning.naive(i).tuneAngleAllCell(mmInd,1) == angles(ai)) / ...
            sum(tuning.naive(i).tuneAngleAllCell(tunedInd,1) == angles(ai));
        propWhisker(i,ai) = sum(tuning.naive(i).tuneAngleAllCell(mmInd,1) == angles(ai)) / ...
            length(mmInd);
        propLeft(i,ai) = sum(tuning.naive(i).tuneAngleAllCell(leftInd,1) == angles(ai)) / ...
            length(leftInd);
        propRaw(i,ai) = sum(tuning.naive(i).tuneAngleAllCell(tunedInd,1) == angles(ai)) / ...
            length(tunedInd);
    end    
end

figure
subplot(131), hold on
bar(angles, nanmean(propPerRaw), 'k')
errorbar(angles, nanmean(propPerRaw), nanstd(propPerRaw)/sqrt(numMice), 'k', 'linestyle', 'none')
xticks(angles)
ylabel('Proportion (gone / raw tuning) per angle')
subplot(132),
errorbar(angles, mean(propWhisker), std(propWhisker)/sqrt(numMice), 'ko-')
xticks(angles)
ylabel('Distribution in gone tuning')
subplot(133), hold on
errorbar(angles, mean(propRaw), std(propRaw)/sqrt(numMice), 'ko-')
errorbar(angles, mean(propLeft), std(propLeft)/sqrt(numMice), 'ro-')
xticks(angles)
ylabel('Distribution')
legend({'Spikes', 'Gone in whisker model'}, 'location', 'northwest')

%% Result: tuned angle distribution remains similar in the remaining cells (still tuned in the whisker model)
% It's good to use what's remained in the whisker model as tuned.

%%     - tuned angle
tunedAngleDiffOut = zeros(numMice, 2); % 1 for touch model, 2 for whisker model
% only from the ones that are remained tuned in each model
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    tunedAngleDiffOut(i,1) = mean(abs(tuning.naive(i).tuneAngleAllCell(indTouch,1) - tuning.naive(i).tuneAngleAllCell(indTouch,2)));
    tunedAngleDiffOut(i,2) = mean(abs(tuning.naive(i).tuneAngleAllCell(indWhisker,1) - tuning.naive(i).tuneAngleAllCell(indWhisker,3)));
end

figure, hold on
bar([1 2], mean(tunedAngleDiffOut), 'k')
errorbar([1 2], mean(tunedAngleDiffOut), std(tunedAngleDiffOut)/sqrt(numMice), 'k', 'linestyle', 'none')
xticks([1 2]), xticklabels({'Touch-angle model', 'Whisker model'}), xtickangle(45)
ylabel('|\DeltaAngle\circ|')
ylim([0 90]), plot([0 3], [15 15], '--', 'color', [0.7 0.7 0.7]), xlim([0 3])

mean(tunedAngleDiffOut)
std(tunedAngleDiffOut)/sqrt(numMice)

%% Result: very minimal change in tuned angle

%%    - tuning modulation
tunedModulation = zeros(numMice, 4); % 1 for raw in touch index, 2 for touch model, 3 for raw in whisker index, 4 for whisker model
% only from the ones that are remained tuned in each model
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    tunedModulation(i,1) = mean(tuning.naive(i).tuneModulationAllCell(indTouch,1));
    tunedModulation(i,2) = mean(tuning.naive(i).tuneModulationAllCell(indTouch,2));
    tunedModulation(i,3) = mean(tuning.naive(i).tuneModulationAllCell(indWhisker,1));
    tunedModulation(i,4) = mean(tuning.naive(i).tuneModulationAllCell(indWhisker,3));
end

figure, hold on
xposition = 1:4;
b = bar(xposition, mean(tunedModulation), 'k');
b.FaceColor = 'flat';
b.CData(2,:) = [1 0 1];
b.CData(4,:) = [1 0 0];
errorbar(xposition, mean(tunedModulation), std(tunedModulation)/sqrt(numMice), 'k', 'linestyle', 'none')
xticks(xposition), xticklabels({'Raw (touch index)', 'Touch model', 'Raw(whisker index)', 'Whisker model'}), xtickangle(45)
ylabel('Modulation')


mean(tunedModulation)
std(tunedModulation)/sqrt(numMice)

%% Result: Huge reduction of modulation in both models, more so in whisker model
%% about half in touch model, about quarter in whisker model

%% How about their distribution?
binWidth = 0.1;
histRange = 0:binWidth:1.6;
modDist = zeros(numMice, length(histRange)-1, 4); % 1 for raw in touch index, 2 for touch model, 3 for raw in whisker index, 4 for whisker model
% only from the ones that are remained tuned in each model
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    modDist(i,:,1) = histcounts(tuning.naive(i).tuneModulationAllCell(indTouch,1), histRange, 'normalization', 'prob');
    modDist(i,:,2) = histcounts(tuning.naive(i).tuneModulationAllCell(indTouch,2), histRange, 'normalization', 'prob');
    modDist(i,:,3) = histcounts(tuning.naive(i).tuneModulationAllCell(indWhisker,1), histRange, 'normalization', 'prob');
    modDist(i,:,4) = histcounts(tuning.naive(i).tuneModulationAllCell(indWhisker,3), histRange, 'normalization', 'prob');
end

figure, 
subplot(211), hold on
temp = squeeze(modDist(:,:,1));
boundedline(histRange(1:end-1) + binWidth/2, mean(temp), sem(temp), 'k')
temp = squeeze(modDist(:,:,2));
boundedline(histRange(1:end-1) + binWidth/2, mean(temp), sem(temp), 'm')
ylabel('Proportion')
subplot(212), hold on
temp = squeeze(modDist(:,:,3));
boundedline(histRange(1:end-1) + binWidth/2, mean(temp), sem(temp), 'k')
temp = squeeze(modDist(:,:,4));
boundedline(histRange(1:end-1) + binWidth/2, mean(temp), sem(temp), 'r')
xlabel('Modulation'), ylabel('Proportion')

%% Results: reduced modulation level



%% Then... how about the correlation?
corrTouch = cell(numMice,1);
corrWhisker = cell(numMice,1);
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrTouch{i} = zeros(length(indTouch),1);
    corrWhisker{i} = zeros(length(indWhisker),1);
    for j = 1 : length(indTouch)
        corrTouch{i}(j) = corr(cellfun(@mean, tuning.naive(i).spkValAllCell{indTouch(j),1}), cellfun(@mean, tuning.naive(i).spkValAllCell{indTouch(j),2}));
    end
    for j = 1 : length(indWhisker)
        corrWhisker{i}(j) = corr(cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),3}));
    end
end

%% plotting mean
figure, hold on
tempTouch = cellfun(@mean, corrTouch);
tempWhisker = cellfun(@mean, corrWhisker);

bar(1, mean(tempTouch), 'm')
bar(2, mean(tempWhisker), 'r')
errorbar(1, mean(tempTouch), sem(tempTouch), 'k')
errorbar(2, mean(tempWhisker), sem(tempWhisker), 'k')
xticks([1 2]), xticklabels({'Touch model', 'Whisker model'}), xtickangle(45)
ylabel('Mean correlation')

%% plotting distribution
histBin = 0.01;
histRange = 0.1:histBin:1;
distCorrTouch = zeros(numMice, length(histRange)-1);
distCorrWhisker = zeros(numMice, length(histRange)-1);
for i = 1 : numMice
    distCorrTouch(i,:) = histcounts(corrTouch{i}, histRange, 'norm', 'cdf');
    distCorrWhisker(i,:) = histcounts(corrWhisker{i}, histRange, 'norm', 'cdf');
end
figure, hold on
boundedline(histRange(1:end-1)+histBin/2, mean(distCorrTouch), sem(distCorrTouch), 'm')
boundedline(histRange(1:end-1)+histBin/2, mean(distCorrWhisker), sem(distCorrWhisker), 'r')
xlabel('Correlation')
ylabel('Cumulative proportion')
plot(histRange([1, end-1])+histBin/2, [0.5 0.5], '--', 'color', [0.7 0.7 0.7]);
plot(histRange([1, end-1])+histBin/2, [0.1 0.1], '--', 'color', [0.7 0.7 0.7]);

%% Compare with permuted random correlation
numPerm = 100;
corrPerm = cell(numMice,1);
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));    
    corrPerm{i} = zeros(length(indTuned),1);
    for ci = 1 : length(indTuned)
        raw = cellfun(@mean,tuning.naive(i).spkValAllCell{indTuned(ci),1});
        spikes = cell2mat(tuning.naive(i).spkValAllCell{indTuned(ci),1});
        groups = [];
        for gi = 1 : length(tuning.naive(i).spkValAllCell{indTuned(ci),1})
            groups = [groups; gi*ones(size(tuning.naive(i).spkValAllCell{indTuned(ci),1}{gi}))];
        end
        tempCorr = zeros(1,numPerm);
        for pi = 1 : numPerm
            tempGroups = groups(randperm(length(groups)));
            tempVal = zeros(7,1);
            for gi = 1 : 7
                tempVal(gi) = mean(spikes(find(tempGroups == gi)));
            end
            tempCorr(pi) = corr(raw, tempVal);
        end
        corrPerm{i}(ci) = mean(tempCorr);
    end
end

%%
figure, hold on
tempTouch = cellfun(@mean, corrTouch);
tempWhisker = cellfun(@mean, corrWhisker);
tempPerm = cellfun(@mean, corrPerm);
bar(1, mean(tempTouch), 'm')
bar(2, mean(tempWhisker), 'r')
bar(3, mean(tempPerm), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(tempTouch), sem(tempTouch), 'k')
errorbar(2, mean(tempWhisker), sem(tempWhisker), 'k')
errorbar(3, mean(tempPerm), sem(tempWhisker), 'k')
xticks([1 2 3]), xticklabels({'Touch model', 'Whisker model', 'Rand perm'}), xtickangle(45)
ylabel('Mean correlation')
        


%% see which one destroyes the correlation the most.
corrFeatures = cell(numMice,12);
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeatures{i,fi} = zeros(length(indWhisker),1);        
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),3+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeatures{i,fi}(j) = tempVal;
        end
    end
end

%%
figure,
tempWhisker = cellfun(@mean, corrWhisker);
tempFeature = cellfun(@mean, corrFeatures);
hold on
errorbar(1, mean(tempWhisker), sem(tempWhisker), 'ro')
errorbar(2:13, mean(tempFeature), sem(tempFeature), 'ko-')
xticks([1:13]), xticklabels({'Full whisker', '-maxDq', '-maxDf', '-maxDkH', '-maxDkV', '-max(slide distance)', '-max(protraction duration)', ...
    '-q', '-f', '-kH', '-kV', '-arc length', '-touch count'})
xtickangle(45)
ylabel('Correlation')
%% see which one recapitulates the most.

corrFeaturesAdd = cell(numMice,12);
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesAdd{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),15+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesAdd{i,fi}(j) = tempVal;
        end
    end
end

%%
figure,
tempWhisker = cellfun(@mean, corrWhisker);
tempFeatureAdd = cellfun(@mean, corrFeaturesAdd);
hold on
errorbar(1, mean(tempWhisker), sem(tempWhisker), 'ro')
errorbar(2:13, mean(tempFeatureAdd), sem(tempFeatureAdd), 'ko-')
xticks([1:13]), xticklabels({'Full whisker', 'maxDq', 'maxDf', 'maxDkH', 'maxDkV', 'max(slide distance)', 'max(protraction duration)', ...
    'q', 'f', 'kH', 'kV', 'arc length', 'touch count'})
xtickangle(45)
ylabel('Correlation')


%% How about combination of them?
loadFns = 'modelAngleTuning_NC_combinations';
data2 = load(sprintf('%s%s',baseDir, loadFns));



%% Which angle is tuned by these? which ones are not? What about their types?

corrFeaturesComb = cell(numMice,8);
for i = 1 : numMice
    indTuned = find(tuning.naive(i).tunedAllCell(:,1));
    indTemp = find(tuning.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 8
        corrFeaturesComb{i,fi} = zeros(length(indWhisker),1);        
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, tuning.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesComb{i,fi}(j) = tempVal;
        end
    end
end

%%
figure,
tempWhisker = cellfun(@mean, corrWhisker);
tempFeature = cellfun(@mean, corrFeaturesComb);
hold on
errorbar(1, mean(tempWhisker), sem(tempWhisker), 'ro')
errorbar(2:9, mean(tempFeature), sem(tempFeature), 'ko-')
xticks([1:9]), xticklabels({'Full whisker', '-(maxDf & maxDkV)', '-(maxDf & maxSD)', '-(maxDkV & maxSD)', '-(maxDf & maxDkV & maxSD)',...
    '(maxDf & maxDkV) only', '(maxDf & maxSD) only', '(maxDkV & maxSD) only', '(maxDf & maxDkV & maxSD) only'})
xtickangle(45)
ylabel('Correlation')



%% With some other combinations







baseDir = 'Y:\Whiskernas\JK\suite2p\';
loadFn1 = 'modelAngleTuning_NC';
loadFn2 = 'modelAngleTuning_NC_combinations';
data1 = load(sprintf('%s%s',baseDir, loadFn1), 'naive'); % only with naive, for now
data2 = load(sprintf('%s%s',baseDir, loadFn2), 'naive'); % 
numMice = length(data1.naive);

%Data1
%         total 27 sets for testing angle tuning
%             1: inferred spikes
%             2: touchGLM
%             3: whiskerGLM
%             4: all-maxDtheta
%             5: all-maxDphi
%             6: all-maxDkH
%             7: all-maxDkV
%             8: all-max(Slide distance)
%             9: all-max(Touch duration)
%             10: all-theta At Touch
%             11: all-phi At Touch
%             12: all-kH At Touch
%             13: all-kV At Touch
%             14: all-arc length At Touch
%             15: all-touch counts
%             16: maxDtheta only (no others included)
%             17: maxDphi only (no others included)
%             18: maxDkH only (no others included)
%             19: maxDkV only (no others included)
%             20: max(slide distance) only (no others included)
%             21: max(touch duration) only (no others included)
%             22: theta at touch only (no others included)
%             23: phi at touch only (no others included)
%             24: kH at touch only (no others included)
%             25: kV at touch only (no others included)
%             26: arc length at touch only (no others included)
%             27: touch counts only (no others included)

%Data2
%         total 25 sets for testing angle tuning
%             1: without any whisker variable
%             2: maxDtheta + other variables
%             3: maxDphi + other variables
%             4: maxDkH + other variables
%             5: maxDkV + other variables
%             6: max(Slide distance) + other variables
%             7: max(duration) + other variables
%             8: thetaAtTouch + other variables
%             9: phiAtTouch + other variables
%             10: kHAtTouch + other variables
%             11: kVAtTouch + other variables
%             12: arc length + other variables
%             13: touch count + other variables
%             14: -(maxDphi + maxDkV)
%             15: -(maxDphi + max(Slide distance))
%             16: -(maxDkV + max(Slide distance))
%             17: -(maxDphi + maxDkV + max(Slide distance))
%             18: dPhi + dKv
%             19: dPhi + slide distance
%             20: dKv + slide distance
%             21: dPhi + dKv + slide distance
%             22: (dPhi + dKv) + other
%             23: (dPhi + slide distance) + other
%             24: (dKv + slide distance) + other
%             25: (maxDphi + maxDkV + max(Slide distance)) + other

% 
% Max model lim from touch GLM
% Min lim from others only
% 
%% (1) drop-out methods
% touch, full, 12 drop-outs, others only (16 total)
corrTouch = cell(numMice,1);
corrWhisker = cell(numMice,1);
for i = 1 : numMice
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrTouch{i} = zeros(length(indTouch),1);
    corrWhisker{i} = zeros(length(indWhisker),1);
    for j = 1 : length(indTouch)
        tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indTouch(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indTouch(j),2}));
        if tempVal < 0
            tempVal = 0;
        end
        corrTouch{i}(j) = tempVal;
    end
    for j = 1 : length(indWhisker)
        tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),3}));
        if tempVal < 0 
            tempVal = 0;
        end
        corrWhisker{i}(j) = tempVal;
    end
end

corrFeaturesOut = cell(numMice,12);
for i = 1 : numMice
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesOut{i,fi} = zeros(length(indWhisker),1);        
        for j = 1 : length(indWhisker)
            if data1.naive(i).tunedAllCell(indWhisker(j),3+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),3+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesOut{i,fi}(j) = tempVal;
        end
    end
    
end

corrOther = cell(numMice,1);
for i = 1 : numMice
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for j = 1 : length(indWhisker)
        if data2.naive(i).tunedAllCell(indWhisker(j),1) == 1
            tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j), 1}));
            if isnan(tempVal) || tempVal < 0
                tempVal = 0;
            end
        else
            tempVal = 0;
        end
            
        corrOther{i}(j) = tempVal;
    end    
end

%%
figure,
tempTouch = cellfun(@mean, corrTouch);
tempWhisker = cellfun(@mean, corrWhisker);
tempFeatureOut = cellfun(@mean, corrFeaturesOut);
tempOther = cellfun(@mean, corrOther);
hold on
bar(1, mean(tempTouch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(tempTouch), sem(tempTouch), 'color', ones(1,3)*0.7)
bar(2, mean(tempWhisker), 'r')
errorbar(2, mean(tempWhisker), sem(tempWhisker), 'r')
bar(3:14, mean(tempFeatureOut), 'k')
errorbar(3:14, mean(tempFeatureOut), sem(tempFeatureOut), 'k', 'linestyle', 'none')
bar(15, mean(tempOther), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(tempOther), sem(tempOther), 'color', ones(1,3)*0.7)
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', '-maxDq', '-maxDf', '-maxDkH', '-maxDkV', '-max(slide distance)', '-max(protraction duration)', ...
    '-q', '-f', '-kH', '-kV', '-arc length', '-touch count', 'others only'})
xtickangle(45)
ylim([0 1])
ylabel('Correlation')
%%
pfull= zeros(12,1);
for i = 1 : 12
    [~, pfull(i)] = ttest(tempWhisker, tempFeatureOut(:,i));
end



%% (2) drop-in method (including all others)
% touch, full, 12 features + others, others only

corrFeaturesIn = cell(numMice,12);
for i = 1 : numMice
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesIn{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            if data2.naive(i).tunedAllCell(indWhisker(j),1+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),1+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesIn{i,fi}(j) = tempVal;
        end
    end
end
%%
figure,
tempTouch = cellfun(@mean, corrTouch);
tempWhisker = cellfun(@mean, corrWhisker);
tempFeatureIn = cellfun(@mean, corrFeaturesIn);
tempOther = cellfun(@mean, corrOther);
hold on
bar(1, mean(tempTouch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(tempTouch), sem(tempTouch), 'color', ones(1,3)*0.7)
bar(2, mean(tempWhisker), 'facecolor', ones(1,3)*0.7)
errorbar(2, mean(tempWhisker), sem(tempWhisker), 'color', ones(1,3)*0.7)
bar(3:14, mean(tempFeatureIn), 'k')
errorbar(3:14, mean(tempFeatureIn), sem(tempFeatureIn), 'k', 'linestyle', 'none')
bar(15, mean(tempOther), 'r')
errorbar(15, mean(tempOther), sem(tempOther), 'r')
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', 'maxDq + others', 'maxDf + others', 'maxDkH + others', 'maxDkV + others', 'max(slide distance) + others', ...
    'max(protraction duration) + others', ...
    'q + others', 'f + others', 'kH + others', 'kV + others', 'arc length + others', 'touch count + others', 'others only'})
xtickangle(45)
ylim([0 1])
ylabel('Correlation')

%%
pmin= zeros(12,1);
for i = 1 : 12
    [~, pmin(i)] = ttest(tempOther, tempFeatureIn(:,i));
end

%% (3) combinations
% touch, full, drop-out combinations, drop-in combinations, others only

% 14~17, 22~25

corrFeaturesCombOut = cell(numMice,4);
corrFeaturesCombIn = cell(numMice,4);
for i = 1 : numMice
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 4
        corrFeaturesCombOut{i,fi} = zeros(length(indWhisker),1);
        corrFeaturesCombIn{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            if data2.naive(i).tunedAllCell(indWhisker(j),13+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),13+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombOut{i,fi}(j) = tempVal;
            
            if data2.naive(i).tunedAllCell(indWhisker(j), 21+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),21+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombIn{i,fi}(j) = tempVal;
            
        end
    end
end
%%
figure,
tempTouch = cellfun(@mean, corrTouch);
tempWhisker = cellfun(@mean, corrWhisker);
tempFeatureOut = cellfun(@mean, corrFeaturesOut(:, [2,5,4]));
tempFeatureCombOut = cellfun(@mean, corrFeaturesCombOut(:,[2,1,3,4]));
tempFeatureIn = cellfun(@mean, corrFeaturesIn(:, [2,5,4]));
tempFeatureCombIn = cellfun(@mean, corrFeaturesCombIn(:,[2,1,3,4]));
tempOther = cellfun(@mean, corrOther);
hold on
bar(1, mean(tempTouch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(tempTouch), sem(tempTouch), 'color', ones(1,3)*0.7)
bar(2, mean(tempWhisker), 'r')
errorbar(2, mean(tempWhisker), sem(tempWhisker), 'r')
bar(3:5, mean(tempFeatureOut), 'k')
errorbar(3:5, mean(tempFeatureOut), sem(tempFeatureOut), 'k', 'linestyle', 'none')
bar(6:9, mean(tempFeatureCombOut), 'k')
errorbar(6:9, mean(tempFeatureCombOut), sem(tempFeatureCombOut), 'k', 'linestyle', 'none')

bar(10, mean(tempOther), 'r')
errorbar(10, mean(tempOther), sem(tempOther), 'r')

bar(11:13, mean(tempFeatureIn), 'k')
errorbar(11:13, mean(tempFeatureIn), sem(tempFeatureIn), 'k', 'linestyle', 'none')
bar(14:17, mean(tempFeatureCombIn), 'k')
errorbar(14:17, mean(tempFeatureCombIn), sem(tempFeatureCombIn), 'k', 'linestyle', 'none')
% xticks([1:9, 11:18]), xticklabels({'Touch model', 'Full whisker', '-maxDf', '-maxDkV', '-max(slide distance)', '-(maxDf & maxDkV)', '-(maxDf & max(SD))', '-(maxDkV & max(SD))', '-(maxDf & maxDkV & max(SD))', ...
%     'maxDf + others', 'maxDkV + others', 'max(SD) + others', '(maxDf & maxDkV) + others', '(maxDf & max(SD)) + others', '(maxDkV & max(SD)) + others', '(maxDf + maxDkV + max(SD)) + others', ...
%     'others only'})
xticks([1:17]), xticklabels({'Touch model', 'Full whisker', '-maxDf', '-max(SD)', '-maxDkV', '-(maxDf & max(SD))', '-(maxDf & maxDkV)', '-(maxDkV & max(SD))', '-(maxDf & maxDkV & max(SD))', ...
    'others only', ...
    'maxDf + others', 'max(SD) + others', 'maxDkV + others', '(maxDf & maxDkV) + others', '(maxDf & max(SD)) + others', '(maxDkV & max(SD)) + others', '(maxDf + maxDkV + max(SD)) + others', ...
    })

xtickangle(45)
ylim([0 1])
ylabel('Correlation')





%% Is it possible to point out which cells loose angle tuning in drop-out methods?
%% (1) look at the ones that maintain pvalue > 0.05 and if they match the tuned angle (or how much the tuned angle differs)

maintainPvalue = cell(numMice, 15); % 1 for touch model, 2 for full whisker model, 3:14 for each individual whisker features, 15 for no whisker features
tunedAngleDiffOut = cell(numMice, 15);
for i = 1 : 12
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    maintainPvalue{i,1} = data1.naive(i).tunedAllCell(indTuned,2);
    tunedAngleDiffOut{i,1} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data1.naive(i).tuneAngleAllCell(indTuned,2); % not tuned ones are NaN
    
    maintainPvalue{i,2} = data1.naive(i).tunedAllCell(indTuned,3);
    tunedAngleDiffOut{i,2} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data1.naive(i).tuneAngleAllCell(indTuned,3); % not tuned ones are NaN
    
    for fi = 1 : 12
        maintainPvalue{i,2+fi} = data1.naive(i).tunedAllCell(indTuned,3+fi);
        tunedAngleDiffOut{i,2+fi} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data1.naive(i).tuneAngleAllCell(indTuned, 3+fi);
    end
    
    maintainPvalue{i,15} = data2.naive(i).tunedAllCell(indTuned,1);
    tunedAngleDiffOut{i,15} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data2.naive(i).tuneAngleAllCell(indTuned,1); % not tuned ones are NaN
    
end

%% Does the number of maintaining tuned match with correlation values?
touchMaintained = cellfun(@mean, maintainPvalue(:,1)); 
fullMaintained = cellfun(@mean, maintainPvalue(:,2));
eachMaintained = cellfun(@mean, maintainPvalue(:,3:14));
othersMaintained = cellfun(@mean, maintainPvalue(:,15));
figure, hold on
bar(1, mean(fullMaintained), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(fullMaintained), sem(fullMaintained), 'color', ones(1,3)*0.7)
bar(2, mean(fullMaintained), 'r')
errorbar(2, mean(fullMaintained), sem(fullMaintained), 'r')
bar(3:14, mean(eachMaintained), 'k')
errorbar(3:14, mean(eachMaintained), sem(eachMaintained), 'k', 'lines','no')
bar(15, mean(othersMaintained), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(othersMaintained), sem(othersMaintained), 'color', ones(1,3)*0.7)

xticks([1:15]), xticklabels({'Touch model', 'Full whisker', '-maxDq', '-maxDf', '-maxDkH', '-maxDkV', '-max(slide distance)', '-max(protraction duration)', ...
    '-q', '-f', '-kH', '-kV', '-arc length', '-touch count', 'others only'})
xtickangle(45)
ylabel('Occurrence')

%% How about tuned angle difference?
fullAD = cellfun(@(x) nanmean(abs(x))/15, tunedAngleDiffOut(:,1));
eachAD = cellfun(@(x) nanmean(abs(x))/15, tunedAngleDiffOut(:,2:13));
othersAD = cellfun(@(x) nanmean(abs(x))/15, tunedAngleDiffOut(:,14));

figure, hold on
bar(1, mean(fullAD), 'r')
errorbar(1, mean(fullAD), sem(fullAD), 'r')
bar(2:13, mean(eachAD), 'k')
errorbar(2:13, mean(eachAD), sem(eachAD), 'k', 'lines','no')
bar(14, mean(othersAD), 'facecolor', ones(1,3)*0.7)
errorbar(14, mean(othersAD), sem(othersAD), 'color', ones(1,3)*0.7)
xticks([1:14]), xticklabels({'Full whisker', 'maxDq + others', 'maxDf + others', 'maxDkH + others', 'maxDkV + others', 'max(slide distance) + others', ...
    'max(protraction duration) + others', ...
    'q + others', 'f + others', 'kH + others', 'kV + others', 'arc length + others', 'touch count + others', 'others only'})
xtickangle(45)
ylabel({'Mean absolute difference'; 'of the tuned angle/ 15\circ'})

%% What if I impose both on the tuning? (maintained tuning should have both anovaP < 0.05 AND matched tuned angle)
% which is just the occurrence of angleDiff == 0
touchMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,1)); 
fullMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,2));
eachMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,3:14));
othersMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,15));

figure, hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:14, mean(eachMatch), 'k')
errorbar(3:14, mean(eachMatch), sem(eachMatch), 'k', 'lines','no')
bar(15, mean(othersMatch), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(othersMatch), sem(othersMatch), 'color', ones(1,3)*0.7)
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', '-maxDq', '-maxDf', '-maxDkH', '-maxDkV', '-max(slide distance)', '-max(protraction duration)', ...
    '-q', '-f', '-kH', '-kV', '-arc length', '-touch count', 'others only'})

xtickangle(45)
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle = 0\circ')

%% what if i include adjacent angle tuning?
touchMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,1));
fullMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,2));
eachMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,3:14));
othersMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,15));

figure, hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:14, mean(eachMatch), 'k')
errorbar(3:14, mean(eachMatch), sem(eachMatch), 'k', 'lines','no')
bar(15, mean(othersMatch), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(othersMatch), sem(othersMatch), 'color', ones(1,3)*0.7)
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', '-maxDq', '-maxDf', '-maxDkH', '-maxDkV', '-max(slide distance)', '-max(protraction duration)', ...
    '-q', '-f', '-kH', '-kV', '-arc length', '-touch count', 'others only'})
xtickangle(45)
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle \leq 15\circ')



%% How about in drop-in methods?

tunedAngleDiffIn = cell(numMice, 12); % 1:12 for each individual whisker features, + others
for i = 1 : 12
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        tunedAngleDiffIn{i,fi} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data2.naive(i).tuneAngleAllCell(indTuned, 1+fi); % not tuned ones are NaN
    end
end


%%
touchMatch = cellfun(@(x) length(find(abs(x) == 0))/length(x), tunedAngleDiffOut(:,1));
fullMatch = cellfun(@(x) length(find(abs(x) == 0))/length(x), tunedAngleDiffOut(:,2));
eachMatchIn = cellfun(@(x) length(find(abs(x) == 0))/length(x), tunedAngleDiffIn);
othersMatch = cellfun(@(x) length(find(abs(x) == 0))/length(x), tunedAngleDiffOut(:,15));

figure, hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:14, mean(eachMatchIn), 'k')
errorbar(3:14, mean(eachMatchIn), sem(eachMatchIn), 'k', 'lines','no')
bar(15, mean(othersMatch), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(othersMatch), sem(othersMatch), 'color', ones(1,3)*0.7)
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', 'maxDq + others', 'maxDf + others', 'maxDkH + others', 'maxDkV + others', 'max(slide distance) + others', ...
    'max(protraction duration) + others', ...
    'q + others', 'f + others', 'kH + others', 'kV + others', 'arc length + others', 'touch count + others', 'others only'})
xtickangle(45)
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle = 0\circ')

%%
touchMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,1));
fullMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,2));
eachMatchIn = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffIn);
othersMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,15));

figure, hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:14, mean(eachMatchIn), 'k')
errorbar(3:14, mean(eachMatchIn), sem(eachMatchIn), 'k', 'lines','no')
bar(15, mean(othersMatch), 'facecolor', ones(1,3)*0.7)
errorbar(15, mean(othersMatch), sem(othersMatch), 'color', ones(1,3)*0.7)
xticks([1:15]), xticklabels({'Touch model', 'Full whisker', 'maxDq + others', 'maxDf + others', 'maxDkH + others', 'maxDkV + others', 'max(slide distance) + others', ...
    'max(protraction duration) + others', ...
    'q + others', 'f + others', 'kH + others', 'kV + others', 'arc length + others', 'touch count + others', 'others only'})
xtickangle(45)
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle \leq 15\circ')


%% 3 top features, both combination and each, both drop-out and -in
% 14~17, 22~25
tunedAngleDiffCombOut = cell(numMice, 12); % 1:12 for each individual whisker features, + others
tunedAngleDiffCombIn = cell(numMice, 12); % 1:12 for each individual whisker features, + others
for i = 1 : 12
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 4
        tunedAngleDiffCombOut{i,fi} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data2.naive(i).tuneAngleAllCell(indTuned, 13+fi); % not tuned ones are NaN
        tunedAngleDiffCombIn{i,fi} = data1.naive(i).tuneAngleAllCell(indTuned,1) - data2.naive(i).tuneAngleAllCell(indTuned, 21+fi); % not tuned ones are NaN
    end
end

%%

figure,
touchMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,1));
fullMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,2));
eachMatchOut = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:, [4,7,6]));
combMatchOut = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffCombOut(:, [2,1,3,4]));
othersMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,15));
eachMatchIn = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffIn(:, [2,5,4]));
combMatchIn = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffCombIn(:, [2,1,3,4]));

hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:5, mean(eachMatchOut), 'k')
errorbar(3:5, mean(eachMatchOut), sem(eachMatchOut), 'k', 'linestyle', 'none')
bar(6:9, mean(combMatchOut), 'k')
errorbar(6:9, mean(combMatchOut), sem(combMatchOut), 'k', 'linestyle', 'none')

bar(10, mean(othersMatch), 'r')
errorbar(10, mean(othersMatch), sem(othersMatch), 'r')

bar(11:13, mean(eachMatchIn), 'k')
errorbar(11:13, mean(eachMatchIn), sem(eachMatchIn), 'k', 'linestyle', 'none')
bar(14:17, mean(combMatchIn), 'k')
errorbar(14:17, mean(combMatchIn), sem(combMatchIn), 'k', 'linestyle', 'none')

xticks([1:17]), xticklabels({'Touch model', 'Full whisker', '-maxDf', '-max(SD)', '-maxDkV', '-(maxDf & max(SD))', '-(maxDf & maxDkV)', '-(maxDkV & max(SD))', '-(maxDf & maxDkV & max(SD))', ...
    'others only', ...
    'maxDf + others', 'max(SD) + others', 'maxDkV + others', '(maxDf & maxDkV) + others', '(maxDf & max(SD)) + others', '(maxDkV & max(SD)) + others', '(maxDf + maxDkV + max(SD)) + others', ...
    })

xtickangle(45)
ylim([0 1])
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle == 0\circ')




%%

figure,
touchMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,1));
fullMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,2));
eachMatchOut = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:, [4,7,6]));
combMatchOut = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffCombOut(:, [2,1,3,4]));
othersMatch = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffOut(:,15));
eachMatchIn = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffIn(:, [2,5,4]));
combMatchIn = cellfun(@(x) length(find(abs(x) <= 15))/length(x), tunedAngleDiffCombIn(:, [2,1,3,4]));

hold on
bar(1, mean(touchMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(touchMatch), sem(touchMatch), 'color', ones(1,3)*0.7)
bar(2, mean(fullMatch), 'r')
errorbar(2, mean(fullMatch), sem(fullMatch), 'r')
bar(3:5, mean(eachMatchOut), 'k')
errorbar(3:5, mean(eachMatchOut), sem(eachMatchOut), 'k', 'linestyle', 'none')
bar(6:9, mean(combMatchOut), 'k')
errorbar(6:9, mean(combMatchOut), sem(combMatchOut), 'k', 'linestyle', 'none')

bar(10, mean(othersMatch), 'r')
errorbar(10, mean(othersMatch), sem(othersMatch), 'r')

bar(11:13, mean(eachMatchIn), 'k')
errorbar(11:13, mean(eachMatchIn), sem(eachMatchIn), 'k', 'linestyle', 'none')
bar(14:17, mean(combMatchIn), 'k')
errorbar(14:17, mean(combMatchIn), sem(combMatchIn), 'k', 'linestyle', 'none')

xticks([1:17]), xticklabels({'Touch model', 'Full whisker', '-maxDf', '-max(SD)', '-maxDkV', '-(maxDf & max(SD))', '-(maxDf & maxDkV)', '-(maxDkV & max(SD))', '-(maxDf & maxDkV & max(SD))', ...
    'others only', ...
    'maxDf + others', 'max(SD) + others', 'maxDkV + others', '(maxDf & maxDkV) + others', '(maxDf & max(SD)) + others', '(maxDkV & max(SD)) + others', '(maxDf + maxDkV + max(SD)) + others', ...
    })

xtickangle(45)
ylim([0 1])
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle \leq 15\circ')


%% Just top 2 features (dKv and slide distance)

%% correlation analysis
figure,
tempWhisker = cellfun(@mean, corrWhisker);
temptempFeatureOutMean = cellfun(@mean, corrFeaturesOut);
tempFeatureOutMean = mean(temptempFeatureOutMean,2);
tempFeatureOut = cellfun(@mean, corrFeaturesOut(:, [5,4]));
tempFeatureCombOut = cellfun(@mean, corrFeaturesCombOut(:,3));
temptempFeatureInMean = cellfun(@mean, corrFeaturesIn);
tempFeatureInMean = mean(temptempFeatureInMean,2);
tempFeatureIn = cellfun(@mean, corrFeaturesIn(:, [5,4]));
tempFeatureCombIn = cellfun(@mean, corrFeaturesCombIn(:,3));
tempOther = cellfun(@mean, corrOther);
hold on
bar(1, mean(tempWhisker), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(tempWhisker), sem(tempWhisker), 'color', ones(1,3)*0.7)
bar(2, mean(tempFeatureOutMean), 'facecolor', ones(1,3)*0.7)
errorbar(2, mean(tempFeatureOutMean), sem(tempFeatureOutMean), 'color', ones(1,3)*0.7)
bar(3:4, mean(tempFeatureOut), 'k')
errorbar(3:4, mean(tempFeatureOut), sem(tempFeatureOut), 'k', 'linestyle', 'none')
bar(5, mean(tempFeatureCombOut), 'k')
errorbar(5, mean(tempFeatureCombOut), sem(tempFeatureCombOut), 'k', 'linestyle', 'none')

bar(6, mean(tempOther), 'facecolor', ones(1,3)*0.7)
errorbar(6, mean(tempOther), sem(tempOther), 'color', ones(1,3)*0.7)
bar(7, mean(tempFeatureInMean), 'facecolor', ones(1,3)*0.7)
errorbar(7, mean(tempFeatureInMean), sem(tempFeatureInMean), 'color', ones(1,3)*0.7)

bar(8:9, mean(tempFeatureIn), 'k')
errorbar(8:9, mean(tempFeatureIn), sem(tempFeatureIn), 'k', 'linestyle', 'none')
bar(10, mean(tempFeatureCombIn), 'k')
errorbar(10, mean(tempFeatureCombIn), sem(tempFeatureCombIn), 'k', 'linestyle', 'none')
xticks([1:10]), xticklabels({'Full whisker', 'mean drop-out', '-max(SD)', '-maxDkV', '-(maxDkV & max(SD))', ...
    'others only', ...
    'mean drop-in + others', 'max(SD) + others', 'maxDkV + others', '(maxDkV & max(SD)) + others'})
xtickangle(45)
ylim([0 1])
ylabel('Correlation')


%% Matched angle analysis

figure,
fullMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,2));
tempEachMatchOutMean = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,3:14));
eachMatchOutMean = mean(tempEachMatchOutMean,2);
eachMatchOut = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:, [7,6]));
combMatchOut = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffCombOut(:, 3));
othersMatch = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffOut(:,15));
tempEachMatchInMean = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffIn);
eachMatchInMean = mean(tempEachMatchInMean,2);
eachMatchIn = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffIn(:, [5,4]));
combMatchIn = cellfun(@(x) length(find(x == 0))/length(x), tunedAngleDiffCombIn(:, 3));

hold on
bar(1, mean(fullMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, mean(fullMatch), sem(fullMatch), 'color', ones(1,3)*0.7)
bar(2, mean(eachMatchOutMean), 'facecolor', ones(1,3)*0.7)
errorbar(2, mean(eachMatchOutMean), sem(eachMatchOutMean), 'color', ones(1,3)*0.7)

bar(3:4, mean(eachMatchOut), 'k')
errorbar(3:4, mean(eachMatchOut), sem(eachMatchOut), 'k', 'linestyle', 'none')
bar(5, mean(combMatchOut), 'k')
errorbar(5, mean(combMatchOut), sem(combMatchOut), 'k', 'linestyle', 'none')

bar(6, mean(othersMatch), 'facecolor', ones(1,3)*0.7)
errorbar(6, mean(othersMatch), sem(othersMatch), 'color', ones(1,3)*0.7)
bar(7, mean(eachMatchInMean), 'facecolor', ones(1,3)*0.7)
errorbar(7, mean(eachMatchInMean), sem(eachMatchInMean), 'color', ones(1,3)*0.7)


bar(8:9, mean(eachMatchIn), 'k')
errorbar(8:9, mean(eachMatchIn), sem(eachMatchIn), 'k', 'linestyle', 'none')
bar(10, mean(combMatchIn), 'k')
errorbar(10, mean(combMatchIn), sem(combMatchIn), 'k', 'linestyle', 'none')

xticks([1:10]), xticklabels({'Full whisker', 'mean drop-out', '-max(SD)', '-maxDkV', '-(maxDkV & max(SD))',...
    'others only', ...
    'mean drop-in', 'max(SD) + others', 'maxDkV + others', '(maxDkV & max(SD)) + others'})

xtickangle(45)
ylim([0 1])
ylabel('Prop. tuning maintained')
title('Tuned & \Deltaangle = 0\circ')




%% 2. Make figures explaining the method 
%     - Raw tuning, model tuning, remove-one tuning, tuning with single feature

% main example: 7th mouse (JK039) cell ID 5097 (401st cell)
% load('Y:\Whiskernas\JK\suite2p\angle_tuning_summary_preAnswer_perTouch_NC.mat')
spkCell = data1.naive(7).spkValAllCell(401,[1,3,4:15]); % 1 spikes, 2 full model, 3:12 drop-out models of 12 features
figure, hold on
plot((cellfun(@mean, spkCell{1})), 'k')
yyaxis right
plot((cellfun(@mean, spkCell{2})), 'b')
plot((cellfun(@mean, spkCell{6})), 'r-')

%%
touchGLM = load('Y:\Whiskernas\JK\suite2p\glmResults_devExp_touch_NC.mat', 'naive');
whiskerGLM = load('Y:\Whiskernas\JK\suite2p\glmResults_devExp_WKV_touchCell_NC.mat', 'naive');
tuning = load('Y:\Whiskernas\JK\suite2p\angle_tuning_summary_preAnswer_perTouch_NC.mat', 'naive');
% find the ones have highest correlation between spikes and full model, while having the smallest correlation when dropping out dKv (6th of the spkCell)
% among JK039 5th plane (369~485)
% while also having high DE in both touch & whisker model...
%%
p5Ind = find(tuning.naive(7).touchID > 5000 & tuning.naive(7).touchID < 6000);
p5Tuned = find(tuning.naive(7).tuned(p5Ind));
cellInds = p5Ind(p5Tuned);
cellIDs = tuning.naive(7).touchID(cellInds);
touchDE = touchGLM.naive(7).allDE(find(ismember(touchGLM.naive(7).cellID, cellIDs)));
whiskerDE = whiskerGLM.naive(7).allDE(cellInds);
spkCell = data1.naive(7).spkValAllCell(cellInds,[1,3,4:15]);

corrVal = zeros(length(cellInds), 2); % (:,1) for spike vs full model, (:,2) for spike vs -dKv model
for i = 1 : length(cellInds)
    corrVal(i,1) = corr(cellfun(@mean, spkCell{i,1}), cellfun(@mean, spkCell{i,2}));
    corrVal(i,2) = corr(cellfun(@mean, spkCell{i,1}), cellfun(@mean, spkCell{i,6}));
end

%% candidates: 14, 18, 26, 35, 43, 75
candInds = [14,18,26,35,43,75]; 
touchDE(candInds)
whiskerDE(candInds)
corrVal(candInds,:)
cellIDs(26)
tuning.naive(7).modulation(cellInds(26),1)


%% results: decided to have cellInds(26) as the example.
% JK039 S01 cell ID 5087.
tempInd = 26;
figure, hold on
plot((cellfun(@mean, spkCell{tempInd,1})), 'k')
yyaxis right
plot((cellfun(@mean, spkCell{tempInd,2})), 'b')
plot((cellfun(@mean, spkCell{tempInd,6})), 'r-')


%% 3. Make figures about the results


