% Check image quality across sessions.
% Before and after learning.
% Use noise level.
% Particularly compare L4 
% From same cells before and after learning.


clear
mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,16],[3,21],[1,17],[1,23],[3,21]};
mi = 1;

baseDir = 'Y:\Whiskernas\JK\suite2p\';
before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));
load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

% be careful. These should NOT be sorted.
matchedIDBefore = us.sessions(1).cellID(find(us.sessions(1).matchedRefCellID));
matchedIDAfter = us.sessions(1).matchedRefCellID(find(us.sessions(1).matchedRefCellID));

%%
indMatchedBefore = zeros(length(matchedIDBefore),1);
indMatchedAfter = zeros(length(matchedIDAfter),1);
for i = 1 : length(matchedIDBefore)
    indMatchedBefore(i) = find(before.u.cellNums == matchedIDBefore(i));
    indMatchedAfter(i) = find(after.u.cellNums == matchedIDAfter(i));
end

noises = zeros(length(matchedIDBefore),2); % (:,1) for before, (:,2) for after
noises(:,1) = before.u.noise(indMatchedBefore);
noises(:,2) = after.u.noise(indMatchedAfter);

figure, histogram(noises(:,2) - noises(:,1))
%%
figure, hold on
histogram(noises(:,1))
histogram(noises(:,2))

%% From all mice
noiseRange = 0:0.02:0.7;
noiseDiffRange = -0.3:0.02:0.3;
histNoiseBefore = zeros(length(mice), length(noiseRange)-1);
histNoiseAfter = zeros(length(mice), length(noiseRange)-1);
histNoiseDiff = zeros(length(mice), length(noiseDiffRange)-1);
for mi = 1 : length(mice)
    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));
    load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    % be careful. These should NOT be sorted.
    matchedIDBefore = us.sessions(1).cellID(find(us.sessions(1).matchedRefCellID));
    matchedIDAfter = us.sessions(1).matchedRefCellID(find(us.sessions(1).matchedRefCellID));

    indMatchedBefore = zeros(length(matchedIDBefore),1);
    indMatchedAfter = zeros(length(matchedIDAfter),1);
    for i = 1 : length(matchedIDBefore)
        indMatchedBefore(i) = find(before.u.cellNums == matchedIDBefore(i));
        indMatchedAfter(i) = find(after.u.cellNums == matchedIDAfter(i));
    end

    noises = zeros(length(matchedIDBefore),2); % (:,1) for before, (:,2) for after
    noises(:,1) = before.u.noise(indMatchedBefore);
    noises(:,2) = after.u.noise(indMatchedAfter);

    histNoiseBefore(mi,:) = histcounts(noises(:,1), noiseRange, 'normalization', 'probability');
    histNoiseAfter(mi,:) = histcounts(noises(:,2), noiseRange, 'normalization', 'probability');
    histNoiseDiff(mi,:) = histcounts(noises(:,2) - noises(:,1), noiseDiffRange, 'normalization', 'probability');    
end
%%
colors = get(gca, 'colororder');
figure,
subplot(121), hold all
errorbar(noiseRange(1:end-1), mean(histNoiseBefore), std(histNoiseBefore)/sqrt(length(mice)))
errorbar(noiseRange(1:end-1), mean(histNoiseAfter), std(histNoiseAfter)/sqrt(length(mice)))
xlabel('Noise'), ylabel('Proportion'), xlim([min(noiseRange) max(noiseRange)])
legend({'Before', 'After'})

subplot(122), hold all
errorbar(noiseDiffRange(1:end-1), mean(histNoiseDiff), std(histNoiseDiff)/sqrt(length(mice)), 'color', colors(3,:))
negFlip = [histNoiseDiff(:,1:(size(histNoiseDiff,2)/2)), flip(histNoiseDiff(:,1:(size(histNoiseDiff,2)/2)),2)];
errorbar(noiseDiffRange(1:end-1), mean(negFlip), std(negFlip)/sqrt(length(mice)), 'color', colors(4,:))
xlabel('\DeltaNoise'), ylabel('Proportion'), xlim([min(noiseDiffRange) max(noiseDiffRange)])

%% Only in L4
% Defined by the depth after learning

noiseRange = 0.1:0.02:0.7;
noiseDiffRange = -0.3:0.02:0.3;
histNoiseBefore = zeros(length(mice), length(noiseRange)-1);
histNoiseAfter = zeros(length(mice), length(noiseRange)-1);
histNoiseDiff = zeros(length(mice), length(noiseDiffRange)-1);
for mi = 1 : length(mice)
    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));
    load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    % be careful. These should NOT be sorted.
    IDL23 = after.u.cellNums(find(after.u.cellDepths < 350));
    newMatchedID = us.sessions(1).matchedRefCellID;
    newMatchedID(find(ismember(newMatchedID,IDL23))) = deal(0);
    
    matchedIDBefore = us.sessions(1).cellID(find(newMatchedID));
    matchedIDAfter = us.sessions(1).matchedRefCellID(find(newMatchedID));
    
    indMatchedBefore = zeros(length(matchedIDBefore),1);
    indMatchedAfter = zeros(length(matchedIDAfter),1);
    
    for i = 1 : length(matchedIDBefore)
        indMatchedBefore(i) = find(before.u.cellNums == matchedIDBefore(i));
        indMatchedAfter(i) = find(after.u.cellNums == matchedIDAfter(i));
    end

    noises = zeros(length(matchedIDBefore),2); % (:,1) for before, (:,2) for after
    noises(:,1) = before.u.noise(indMatchedBefore);
    noises(:,2) = after.u.noise(indMatchedAfter);

    histNoiseBefore(mi,:) = histcounts(noises(:,1), noiseRange, 'normalization', 'probability');
    histNoiseAfter(mi,:) = histcounts(noises(:,2), noiseRange, 'normalization', 'probability');
    histNoiseDiff(mi,:) = histcounts(noises(:,2) - noises(:,1), noiseDiffRange, 'normalization', 'probability');    
end

colors = get(gca, 'colororder');
figure,
subplot(121), hold all
errorbar(noiseRange(1:end-1), mean(histNoiseBefore), std(histNoiseBefore)/sqrt(length(mice)))
errorbar(noiseRange(1:end-1), mean(histNoiseAfter), std(histNoiseAfter)/sqrt(length(mice)))
xlabel('Noise'), ylabel('Proportion'), xlim([min(noiseRange) max(noiseRange)])
legend({'Before', 'After'})

subplot(122), hold all
errorbar(noiseDiffRange(1:end-1), mean(histNoiseDiff), std(histNoiseDiff)/sqrt(length(mice)), 'color', colors(3,:))
negFlip = [histNoiseDiff(:,1:(size(histNoiseDiff,2)/2)), flip(histNoiseDiff(:,1:(size(histNoiseDiff,2)/2)),2)];
errorbar(noiseDiffRange(1:end-1), mean(negFlip), std(negFlip)/sqrt(length(mice)), 'color', colors(4,:))
xlabel('\DeltaNoise'), ylabel('Proportion'), xlim([min(noiseDiffRange) max(noiseDiffRange)])

 
%% Results: 
% There are a little bit of increase in noise level after learning, compared to before
% So... is the reduced proportion of touch response cells in L4 after
% learning affected by this?

% Do the noise matching, same as in d190531_sparse_L4.

clear
mice = [25,27,30,36,39,52];
indMice = [1,2,3,4,6,9]; % for cellFunctionRidgeDE010
sessions = {[4,19],[3,16],[3,21],[1,17],[1,23],[3,21]};

% mice = [25,27,30,36,52];
% indMice = [1,2,3,4,9]; % for cellFunctionRidgeDE010
% sessions = {[4,19],[3,16],[3,21],[1,17],[3,21]};

baseDir = 'Y:\Whiskernas\JK\suite2p\';
load(sprintf('%scellFunctionRidgeDE010', baseDir), 'naive', 'expert')
noiseBinNum = 10;
prctileThresh = 5;

numTouchL23 = zeros(2,length(mice),noiseBinNum);
numTouchL4 = zeros(2,length(mice),noiseBinNum);

for mi = 1 : length(mice)

    % load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));

    indL23{1} = find(before.u.cellDepths < 350);
    indL23{2} = find(after.u.cellDepths < 350);

    indL4{1} = find(before.u.cellDepths >= 350);
    indL4{2} = find(after.u.cellDepths >= 350);

    noiseL23 = [prctile(union(before.u.noise(indL23{1}), after.u.noise(indL23{2})), prctileThresh), prctile(union(before.u.noise(indL23{1}), after.u.noise(indL23{2})), 100-prctileThresh)];
    noiseL4 = [prctile(union(before.u.noise(indL4{1}), after.u.noise(indL4{2})), prctileThresh), prctile(union(before.u.noise(indL4{1}), after.u.noise(indL4{2})), 100-prctileThresh)];

    noiseRangeL23 = min(noiseL23) : (max(noiseL23) - min(noiseL23)) / noiseBinNum : max(noiseL23);
    noiseRangeL4 = min(noiseL4) : (max(noiseL4) - min(noiseL4)) / noiseBinNum : max(noiseL4);

    for bi = 1 : noiseBinNum
        indBeforeNoiseL23{bi} = indL23{1}(find(before.u.noise(indL23{1}) >= noiseRangeL23(bi) & before.u.noise(indL23{1}) < noiseRangeL23(bi+1)));
        indAfterNoiseL23{bi} = indL23{2}(find(after.u.noise(indL23{2}) >= noiseRangeL23(bi) & after.u.noise(indL23{2}) < noiseRangeL23(bi+1)));
        numCellL23{mi}(bi) = min( [length(indBeforeNoiseL23{bi}), length(indAfterNoiseL23{bi})] );

        indBeforeNoiseL4{bi} = indL4{1}(find(before.u.noise(indL4{1}) >= noiseRangeL4(bi) & before.u.noise(indL4{1}) < noiseRangeL4(bi+1)));
        indAfterNoiseL4{bi} = indL4{2}(find(after.u.noise(indL4{2}) >= noiseRangeL4(bi) & after.u.noise(indL4{2}) < noiseRangeL4(bi+1)));
        numCellL4{mi}(bi) = min( [length(indBeforeNoiseL4{bi}), length(indAfterNoiseL4{bi})] );

        tempIndBeforeNoiseL23 = indBeforeNoiseL23{bi}(randperm(length(indBeforeNoiseL23{bi}), numCellL23{mi}(bi)));
        tempIndAfterNoiseL23 = indAfterNoiseL23{bi}(randperm(length(indAfterNoiseL23{bi}), numCellL23{mi}(bi)));

        tempIndBeforeNoiseL4 = indBeforeNoiseL4{bi}(randperm(length(indBeforeNoiseL4{bi}), numCellL4{mi}(bi)));
        tempIndAfterNoiseL4 = indAfterNoiseL4{bi}(randperm(length(indAfterNoiseL4{bi}), numCellL4{mi}(bi)));

        numTouchL23(1,mi,bi) = sum(ismember(before.u.cellNums(tempIndBeforeNoiseL23), naive(indMice(mi)).touchID));
        numTouchL23(2,mi,bi) = sum(ismember(after.u.cellNums(tempIndAfterNoiseL23), expert(mi).touchID));

        numTouchL4(1,mi,bi) = sum(ismember(before.u.cellNums(tempIndBeforeNoiseL4), naive(indMice(mi)).touchID));
        numTouchL4(2,mi,bi) = sum(ismember(after.u.cellNums(tempIndAfterNoiseL4), expert(mi).touchID));
    end

    propTouchL23(mi,1) = sum(numTouchL23(1,mi,:)) / sum(numCellL23{mi});
    propTouchL23(mi,2) = sum(numTouchL23(2,mi,:)) / sum(numCellL23{mi});

    propTouchL4(mi,1) = sum(numTouchL4(1,mi,:)) / sum(numCellL4{mi});
    propTouchL4(mi,2) = sum(numTouchL4(2,mi,:)) / sum(numCellL4{mi});
end

%%
figure, hold all
bar(0.8, mean(propTouchL23(:,1)), 0.4, 'w')
bar(1.2, mean(propTouchL23(:,2)), 0.4, 'k')
errorbar(0.8, mean(propTouchL23(:,1)), std(propTouchL23(:,1))/sqrt(length(mice)), 'k.')
errorbar(1.2, mean(propTouchL23(:,2)), std(propTouchL23(:,2))/sqrt(length(mice)), 'k.')
bar(1.8, mean(propTouchL4(:,1)), 0.4, 'w')
bar(2.2, mean(propTouchL4(:,2)), 0.4, 'k')
errorbar(1.8, mean(propTouchL4(:,1)), std(propTouchL4(:,1))/sqrt(length(mice)), 'k.')
errorbar(2.2, mean(propTouchL4(:,2)), std(propTouchL4(:,2))/sqrt(length(mice)), 'k.')
legend({'Naive', 'Expert'})
xticks([1, 2]), xticklabels({'L2/3', 'L4'}), ylabel('Proportion')
title('Touch cells / active cells')

%%
[~, p23] = ttest(propTouchL23(:,1), propTouchL23(:,2))
[~, p4] = ttest(propTouchL4(:,1), propTouchL4(:,2))

%% Results: There is a tendency to decrease only in L4, but not statistically significant.
% What about matched cells?

clear
mice = [25,27,30,36,39,52];
indMice = [1,2,3,4,6,9]; % for cellFunctionRidgeDE010
sessions = {[4,19],[3,16],[3,21],[1,17],[1,23],[3,21]};
baseDir = 'Y:\Whiskernas\JK\suite2p\';
load(sprintf('%scellFunctionRidgeDE010', baseDir), 'naive', 'expert')
noiseBinNum = 10;
prctileThresh = 5;

for mi = 1 : length(mice)

    load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));

    indMatchedBefore = find(us.sessions(1).matchedRefCellID);    
    indMatchedAfter = find(ismember(us.sessions(2).cellID, us.sessions(1).matchedRefCellID));
    
    indL23{1} = intersect(find(before.u.cellDepths < 350), indMatchedBefore);
    indL23{2} = intersect(find(after.u.cellDepths < 350), indMatchedAfter);

    indL4{1} = intersect(find(before.u.cellDepths >= 350), indMatchedBefore);
    indL4{2} = intersect(find(after.u.cellDepths >= 350), indMatchedAfter);

    noiseL23 = [prctile(union(before.u.noise(indL23{1}), after.u.noise(indL23{2})), prctileThresh), prctile(union(before.u.noise(indL23{1}), after.u.noise(indL23{2})), 100-prctileThresh)];
    noiseL4 = [prctile(union(before.u.noise(indL4{1}), after.u.noise(indL4{2})), prctileThresh), prctile(union(before.u.noise(indL4{1}), after.u.noise(indL4{2})), 100-prctileThresh)];

    noiseRangeL23 = min(noiseL23) : (max(noiseL23) - min(noiseL23)) / noiseBinNum : max(noiseL23);
    noiseRangeL4 = min(noiseL4) : (max(noiseL4) - min(noiseL4)) / noiseBinNum : max(noiseL4);

    numTouchL23 = zeros(2,1,noiseBinNum);
    numTouchL4 = zeros(2,1,noiseBinNum);

    for bi = 1 : noiseBinNum
        indBeforeNoiseL23{bi} = indL23{1}(find(before.u.noise(indL23{1}) >= noiseRangeL23(bi) & before.u.noise(indL23{1}) < noiseRangeL23(bi+1)));
        indAfterNoiseL23{bi} = indL23{2}(find(after.u.noise(indL23{2}) >= noiseRangeL23(bi) & after.u.noise(indL23{2}) < noiseRangeL23(bi+1)));
        numCellL23(bi) = min( [length(indBeforeNoiseL23{bi}), length(indAfterNoiseL23{bi})] );

        indBeforeNoiseL4{bi} = indL4{1}(find(before.u.noise(indL4{1}) >= noiseRangeL4(bi) & before.u.noise(indL4{1}) < noiseRangeL4(bi+1)));
        indAfterNoiseL4{bi} = indL4{2}(find(after.u.noise(indL4{2}) >= noiseRangeL4(bi) & after.u.noise(indL4{2}) < noiseRangeL4(bi+1)));
        numCellL4(bi) = min( [length(indBeforeNoiseL4{bi}), length(indAfterNoiseL4{bi})] );

        tempIndBeforeNoiseL23 = indBeforeNoiseL23{bi}(randperm(length(indBeforeNoiseL23{bi}), numCellL23(bi)));
        tempIndAfterNoiseL23 = indAfterNoiseL23{bi}(randperm(length(indAfterNoiseL23{bi}), numCellL23(bi)));

        tempIndBeforeNoiseL4 = indBeforeNoiseL4{bi}(randperm(length(indBeforeNoiseL4{bi}), numCellL4(bi)));
        tempIndAfterNoiseL4 = indAfterNoiseL4{bi}(randperm(length(indAfterNoiseL4{bi}), numCellL4(bi)));

        numTouchL23(1,1,bi) = sum(ismember(before.u.cellNums(tempIndBeforeNoiseL23), naive(indMice(mi)).touchID));
        numTouchL23(2,1,bi) = sum(ismember(after.u.cellNums(tempIndAfterNoiseL23), expert(mi).touchID));

        numTouchL4(1,1,bi) = sum(ismember(before.u.cellNums(tempIndBeforeNoiseL4), naive(indMice(mi)).touchID));
        numTouchL4(2,1,bi) = sum(ismember(after.u.cellNums(tempIndAfterNoiseL4), expert(mi).touchID));
    end

    propTouchL23(mi,1) = sum(numTouchL23(1,1,:)) / sum(numCellL23);
    propTouchL23(mi,2) = sum(numTouchL23(2,1,:)) / sum(numCellL23);

    propTouchL4(mi,1) = sum(numTouchL4(1,1,:)) / sum(numCellL4);
    propTouchL4(mi,2) = sum(numTouchL4(2,1,:)) / sum(numCellL4);
end
%%
figure, hold all
bar(0.8, mean(propTouchL23(:,1)), 0.4, 'w')
bar(1.2, mean(propTouchL23(:,2)), 0.4, 'k')
errorbar(0.8, mean(propTouchL23(:,1)), std(propTouchL23(:,1))/sqrt(length(mice)), 'k.')
errorbar(1.2, mean(propTouchL23(:,2)), std(propTouchL23(:,2))/sqrt(length(mice)), 'k.')
bar(1.8, mean(propTouchL4(:,1)), 0.4, 'w')
bar(2.2, mean(propTouchL4(:,2)), 0.4, 'k')
errorbar(1.8, mean(propTouchL4(:,1)), std(propTouchL4(:,1))/sqrt(length(mice)), 'k.')
errorbar(2.2, mean(propTouchL4(:,2)), std(propTouchL4(:,2))/sqrt(length(mice)), 'k.')
legend({'Naive', 'Expert'})
xticks([1, 2]), xticklabels({'L2/3', 'L4'}), ylabel('Proportion')
title('Touch cells / active cells (matched cells)')
%%
[~, p23] = ttest(propTouchL23(:,1), propTouchL23(:,2))
[~, p4] = ttest(propTouchL4(:,1), propTouchL4(:,2))

%% Results: Not significantly reduced touch response cells in L4 among matched cells.
% Maybe because many of the matched cells are touch response cells.
% So... What is the proportion of being detected (active) after learning, between touch response and not (before)?
% Dividing into L2/3 and L4.

clear
mice = [25,27,30,36,39,52];
indMice = [1,2,3,4,6,9]; % for cellFunctionRidgeDE010
sessions = {[4,19],[3,16],[3,21],[1,17],[1,23],[3,21]};
baseDir = 'Y:\Whiskernas\JK\suite2p\';
load(sprintf('%scellFunctionRidgeDE010', baseDir), 'naive', 'expert')

touchMatched = zeros(length(mice),1);
nontouchMatched = zeros(length(mice),1);

for mi = 1 : length(mice)
    
    load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));
    
    indTouchBeforeL23 = intersect(find(before.u.cellDepths < 350), find(ismember(before.u.cellNums, naive(indMice(mi)).touchID)));
    indTouchBeforeL4 = intersect(find(before.u.cellDepths >= 350), find(ismember(before.u.cellNums, naive(indMice(mi)).touchID)));
    indNonTouchBeforeL23 = intersect(find(before.u.cellDepths < 350), find(1-ismember(before.u.cellNums, naive(indMice(mi)).touchID)));
    indNonTouchBeforeL4 = intersect(find(before.u.cellDepths >= 350), find(1-ismember(before.u.cellNums, naive(indMice(mi)).touchID)));
    
    touchMatchedL23(mi) = length(find(us.sessions(1).matchedRefCellID(indTouchBeforeL23))) / length(indTouchBeforeL23);
    nonTouchMatchedL23(mi) = length(find(us.sessions(1).matchedRefCellID(indNonTouchBeforeL23))) / length(indNonTouchBeforeL23);
    touchMatchedL4(mi) = length(find(us.sessions(1).matchedRefCellID(indTouchBeforeL4))) / length(indTouchBeforeL4);
    nonTouchMatchedL4(mi) = length(find(us.sessions(1).matchedRefCellID(indNonTouchBeforeL4))) / length(indNonTouchBeforeL4);
    
end

%%
figure, hold all
bar(0.8, mean(touchMatchedL23), 0.4, 'w')
bar(1.2, mean(nonTouchMatchedL23), 0.4, 'k')
errorbar(0.8, mean(touchMatchedL23), std(touchMatchedL23)/sqrt(length(mice)), 'k.')
errorbar(1.2, mean(nonTouchMatchedL23), std(nonTouchMatchedL23)/sqrt(length(mice)), 'k.')
bar(1.8, mean(touchMatchedL4), 0.4, 'w')
bar(2.2, mean(nonTouchMatchedL4), 0.4, 'k')
errorbar(1.8, mean(touchMatchedL4), std(touchMatchedL4)/sqrt(length(mice)), 'k.')
errorbar(2.2, mean(nonTouchMatchedL4), std(nonTouchMatchedL4)/sqrt(length(mice)), 'k.')
legend({'Touch', 'Non touch'})
xticks([1, 2]), xticklabels({'L2/3', 'L4'}), xlim([0.5 2.5]), ylabel('Proportion')
title('Remained cells')

%% Results:
% A little bit different trend between L2/3 and L4, but none of them are
% significant (paired t-test)

%% Conclusion:
% There are a little bit of increase in noise level after learning, compared to before




%% P.S.
% What if I don't consider noise at all?


clear
mice = [25,27,30,36,39,52];
indMice = [1,2,3,4,6,9]; % for cellFunctionRidgeDE010
sessions = {[4,19],[3,16],[3,21],[1,17],[1,23],[3,21]};

% mice = [25,27,30,36,52];
% indMice = [1,2,3,4,9]; % for cellFunctionRidgeDE010
% sessions = {[4,19],[3,16],[3,21],[1,17],[3,21]};

baseDir = 'Y:\Whiskernas\JK\suite2p\';
load(sprintf('%scellFunctionRidgeDE010', baseDir), 'naive', 'expert')

for mi = 1 : length(mice)

    % load(sprintf('%s%03d\\cellIDmatch_JK%03d', baseDir, mice(mi), mice(mi)))

    before = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
    after = load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(2)));

    indL23{1} = find(before.u.cellDepths < 350);
    indL23{2} = find(after.u.cellDepths < 350);

    indL4{1} = find(before.u.cellDepths >= 350);
    indL4{2} = find(after.u.cellDepths >= 350);

    
    numTouchL23(1,mi) = sum(ismember(before.u.cellNums(indL23{1}), naive(indMice(mi)).touchID));
    numTouchL23(2,mi) = sum(ismember(after.u.cellNums(indL23{2}), expert(mi).touchID));

    numTouchL4(1,mi) = sum(ismember(before.u.cellNums(indL4{1}), naive(indMice(mi)).touchID));
    numTouchL4(2,mi) = sum(ismember(after.u.cellNums(indL4{2}), expert(mi).touchID));

    
    propTouchL23(mi,1) = numTouchL23(1,mi) / length(indL23{1});
    propTouchL23(mi,2) = numTouchL23(2,mi) / length(indL23{2});

    propTouchL4(mi,1) = numTouchL4(1,mi) / length(indL4{1});
    propTouchL4(mi,2) = numTouchL4(2,mi) / length(indL4{2});
end

%%
figure, hold all
bar(0.8, mean(propTouchL23(:,1)), 0.4, 'w')
bar(1.2, mean(propTouchL23(:,2)), 0.4, 'k')
errorbar(0.8, mean(propTouchL23(:,1)), std(propTouchL23(:,1))/sqrt(length(mice)), 'k.')
errorbar(1.2, mean(propTouchL23(:,2)), std(propTouchL23(:,2))/sqrt(length(mice)), 'k.')
bar(1.8, mean(propTouchL4(:,1)), 0.4, 'w')
bar(2.2, mean(propTouchL4(:,2)), 0.4, 'k')
errorbar(1.8, mean(propTouchL4(:,1)), std(propTouchL4(:,1))/sqrt(length(mice)), 'k.')
errorbar(2.2, mean(propTouchL4(:,2)), std(propTouchL4(:,2))/sqrt(length(mice)), 'k.')
legend({'Naive', 'Expert'})
xticks([1, 2]), xticklabels({'L2/3', 'L4'}), ylabel('Proportion')
title('Touch cells / active cells')

%%
[~, p23] = ttest(propTouchL23(:,1), propTouchL23(:,2))
[~, p4] = ttest(propTouchL4(:,1), propTouchL4(:,2))


%% Results: It was not statistically significant even before matching noise...




