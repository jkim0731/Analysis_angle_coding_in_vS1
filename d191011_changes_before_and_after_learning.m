%% What changes and what does not change, before and after learning?
%% Appearred cells, disappeared cells, matching cells
%% Touch or non-touch, tuned or not, tuned angle
%% Across depth


%% basic data loading and settings

clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
matchFn = 'cellMatching_beforeNafter.mat';
functionFn = 'cellFunctionLasso_NC';
angleTuningFn = 'angle_tuning_summary_preAnswer_perTouch_NC.mat';

load(sprintf('%s%s',baseDir, matchFn), 'match');
cellFunction = load(sprintf('%s%s',baseDir, functionFn), 'expert', 'naive');
tune = load(sprintf('%s%s',baseDir, angleTuningFn), 'expert', 'naive');

learnerInd = [1,2,3,4,7,9];
numMice = size(match,1);

%% How much does # of active cell change?
numActive = cellfun(@length, match(:,[1,3]));

propChange = numActive(:,2)./numActive(:,1);
mean(propChange)
sem(propChange)
[~,p] = ttest(propChange-1)

%% # of active cells decreases, by about 12 % (88.32 % /+- 3.88 %; compared to naive. p = 0.0297)
%% Where to they come from? L4, or L2/3?

numActiveNaive = zeros(numMice,4); % 1 L2/3 C2, 2 L2/3 non-C2, 3 L4 C2, 4 L4 non-C2
numActiveExpert = zeros(numMice,4); % 1 L2/3 C2, 2 L2/3 non-C2, 3 L4 C2, 4 L4 non-C2
for mi = 1 : numMice
    idActive = match{mi,1};
    indActive = find(ismember(cellFunction.naive(learnerInd(mi)).allCell.cellNums, idActive));
    depths = cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indActive);
    indL23 = indActive(find(depths < 350));
    indL4 = indActive(find(depths >= 350));    
    isC2 = cellFunction.naive(learnerInd(mi)).allCell.isC2(indActive);
    indC2 = indActive(find(isC2 == 1));
    indNonC2 = indActive(find(isC2 == 0));
    
    numActiveNaive(mi,1) = length(intersect(indL23, indC2));
    numActiveNaive(mi,2) = length(intersect(indL23, indNonC2));
    numActiveNaive(mi,3) = length(intersect(indL4, indC2));
    numActiveNaive(mi,4) = length(intersect(indL4, indNonC2));

    idActive = match{mi,3};
    indActive = find(ismember(cellFunction.expert(mi).allCell.cellNums, idActive));
    depths = cellFunction.expert(mi).allCell.cellDepths(indActive);
    indL23 = indActive(find(depths < 350));
    indL4 = indActive(find(depths >= 350));    
    isC2 = cellFunction.naive(learnerInd(mi)).allCell.isC2(indActive);
    indC2 = indActive(find(isC2 == 1));
    indNonC2 = indActive(find(isC2 == 0));
    
    numActiveExpert(mi,1) = length(intersect(indL23, indC2));
    numActiveExpert(mi,2) = length(intersect(indL23, indNonC2));
    numActiveExpert(mi,3) = length(intersect(indL4, indC2));
    numActiveExpert(mi,4) = length(intersect(indL4, indNonC2));
end

barWidth = 0.4;
posAdj = 0.2;
figure, hold on
bar([1:4]-posAdj, mean(numActiveNaive), barWidth, 'k')
bar([1:4]+posAdj, mean(numActiveExpert), barWidth, 'c')
errorbar([1:4]-posAdj, mean(numActiveNaive), sem(numActiveNaive), 'k', 'lines', 'no')
errorbar([1:4]+posAdj, mean(numActiveExpert), sem(numActiveExpert), 'k', 'lines', 'no')
legend({'Naive', 'Expert'}, 'box', 'off')
xticks(1:4)
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})

%%
tempMat = numActiveExpert./numActiveNaive;
figure, hold on
bar([1:4], mean(tempMat), 'k')
errorbar([1:4], mean(tempMat), sem(tempMat), 'k', 'lines', 'no')
xvals = xlim();
plot(xvals, [1; 1], '--', 'color', ones(1,3)*0.7)
xticks(1:4)
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
ylabel('Ratio # active cells (expert / naive)')
[~,pvals] = ttest(tempMat-1)


%% decrese only in L2/3, not in L4, but all of them shows no significance.
%% what if I don't divide into C2 and non-C2?


numActiveNaive = zeros(numMice,2); % 1: L2/3, 2: L4 
numActiveExpert = zeros(numMice,2); % 1: L2/3, 2: L4
for mi = 1 : numMice
    idActive = match{mi,1};
    indActive = find(ismember(cellFunction.naive(learnerInd(mi)).allCell.cellNums, idActive));
    depths = cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indActive);
    indL23 = indActive(find(depths < 350));
    indL4 = indActive(find(depths >= 350));    
    
    numActiveNaive(mi,1) = length(indL23);
    numActiveNaive(mi,2) = length(indL4);
    
    idActive = match{mi,3};
    indActive = find(ismember(cellFunction.expert(mi).allCell.cellNums, idActive));
    depths = cellFunction.expert(mi).allCell.cellDepths(indActive);
    indL23 = indActive(find(depths < 350));
    indL4 = indActive(find(depths >= 350));    
    
    numActiveExpert(mi,1) = length(indL23);
    numActiveExpert(mi,2) = length(indL4);
end

tempMat = numActiveExpert./numActiveNaive;
figure, hold on
bar([1:2], mean(tempMat), 'k')
errorbar([1:2], mean(tempMat), sem(tempMat), 'k', 'lines', 'no')
xvals = xlim();
plot(xvals, [1; 1], '--', 'color', ones(1,3)*0.7)
xticks(1:2)
xticklabels({'L2/3', 'L4'})
ylabel('Ratio # active cells (expert / naive)')
[~,pvals] = ttest(tempMat-1)


%% Is there a difference in touch-responsiveness in appearred, disappearred, and matching cells?
%% How about whisking cells?


propTouch = zeros(numMice,4); % 1: naive disappearred, 2: naive remained, 3: expert remained, 4: expert appearred
propWhisking = zeros(numMice,4);
propMixed = zeros(numMice,4);

for mi = 1 : numMice
    idNaiveDis = match{mi,1}(find(match{mi,2} == 0));
    idNaiveRem = match{mi,1}(find(match{mi,2}));
    naiveTouch = cellFunction.naive(learnerInd(mi)).touchID;
    naiveWhisking = cellFunction.naive(learnerInd(mi)).whiskingID;
    naiveMixed = intersect(naiveTouch, naiveWhisking);
    
    propTouch(mi,1) = length(intersect(idNaiveDis, naiveTouch)) / length(idNaiveDis);
    propTouch(mi,2) = length(intersect(idNaiveRem, naiveTouch)) / length(idNaiveRem);
    propWhisking(mi,1) = length(intersect(idNaiveDis, naiveWhisking)) / length(idNaiveDis);
    propWhisking(mi,2) = length(intersect(idNaiveRem, naiveWhisking)) / length(idNaiveRem);
    propMixed(mi,1) = length(intersect(idNaiveDis, naiveMixed)) / length(idNaiveDis);
    propMixed(mi,2) = length(intersect(idNaiveRem, naiveMixed)) / length(idNaiveRem);
    
    
    idExpertRem = match{mi,2}(find(match{mi,2}));
    idExpertApp = setdiff(match{mi,3}, idExpertRem);
    
    expertTouch = cellFunction.expert(mi).touchID;
    expertWhisking = cellFunction.expert(mi).whiskingID;
    expertMixed = intersect(expertTouch, expertWhisking);
    
    propTouch(mi,3) = length(intersect(idExpertRem, expertTouch)) / length(idExpertRem);
    propTouch(mi,4) = length(intersect(idExpertApp, expertTouch)) / length(idExpertApp);
    propWhisking(mi,3) = length(intersect(idExpertRem, expertWhisking)) / length(idExpertRem);
    propWhisking(mi,4) = length(intersect(idExpertApp, expertWhisking)) / length(idExpertApp);
    propMixed(mi,3) = length(intersect(idExpertRem, expertMixed)) / length(idExpertRem);
    propMixed(mi,4) = length(intersect(idExpertApp, expertMixed)) / length(idExpertApp);
end

%%
figure
subplot(131), hold on
bar(mean(propTouch),'k')
errorbar(mean(propTouch), sem(propTouch), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Proportion')
title('Touch')


subplot(132), hold on
bar(mean(propWhisking),'k')
errorbar(mean(propWhisking), sem(propWhisking), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Proportion')
title('Whisking')



subplot(133), hold on
bar(mean(propMixed),'k')
errorbar(mean(propMixed), sem(propMixed), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Proportion')
title('Mixed')

%% Not so much difference in proportion of functions between sub populations.

%% What if I just sum them up?
propTouch = zeros(numMice,2); % 1: naive, 2: expert
propWhisking = zeros(numMice,2);
propMixed = zeros(numMice,2);

for mi = 1 : numMice
    naiveTouch = cellFunction.naive(learnerInd(mi)).touchID;
    naiveWhisking = cellFunction.naive(learnerInd(mi)).whiskingID;
    naiveMixed = intersect(naiveTouch, naiveWhisking);
    
    propTouch(mi,1) = length(intersect(match{mi,1},naiveTouch)) / length(match{mi,1});
    propWhisking(mi,1) = length(intersect(match{mi,1},naiveWhisking)) / length(match{mi,1});
    propMixed(mi,1) = length(intersect(match{mi,1},naiveMixed)) / length(match{mi,1});

    expertTouch = cellFunction.expert(mi).touchID;
    expertWhisking = cellFunction.expert(mi).whiskingID;
    expertMixed = intersect(expertTouch, expertWhisking);
    
    propTouch(mi,2) = length(intersect(match{mi,3},expertTouch)) / length(match{mi,3});
    propWhisking(mi,2) = length(intersect(match{mi,3},expertWhisking)) / length(match{mi,3});
    propMixed(mi,2) = length(intersect(match{mi,3},expertMixed)) / length(match{mi,3});
end

figure
subplot(131), hold on
bar(mean(propTouch),'k')
errorbar(mean(propTouch), sem(propTouch), 'k', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
ylabel('Proportion')
title('Touch')


subplot(132), hold on
bar(mean(propWhisking),'k')
errorbar(mean(propWhisking), sem(propWhisking), 'k', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
ylabel('Proportion')
title('Whisking')



subplot(133), hold on
bar(mean(propMixed),'k')
errorbar(mean(propMixed), sem(propMixed), 'k', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
ylabel('Proportion')
title('Mixed')

%% In case of touch, show each mouse
figure, hold on
for i = 1 : size(propTouch,1)
    plot(propTouch(i,:), 'ko-')
end
errorbar(mean(propTouch), sem(propTouch), 'ro', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
xlim([0.5 2.5])
ylabel('Proportion (/active cells)')
title('Touch cells')


%%
[~,p] = ttest(propTouch(:,1), propTouch(:,2))
%% How about tuned response?

propTuned = zeros(numMice,4); % 1: naive disappearred, 2: naive remained, 3: expert remained, 4: expert appearred
for mi = 1 : numMice
    idNaiveDis = match{mi,1}(find(match{mi,2} == 0));
    idNaiveRem = match{mi,1}(find(match{mi,2}));
    naiveTouch = tune.naive(learnerInd(mi)).touchID;
    naiveTuned = naiveTouch(find(tune.naive(learnerInd(mi)).tuned));
    
    propTuned(mi,1) = length(intersect(idNaiveDis, naiveTuned)) / length(intersect(idNaiveDis, naiveTouch));
    propTuned(mi,2) = length(intersect(idNaiveRem, naiveTuned)) / length(intersect(idNaiveRem, naiveTouch));
    
    idExpertRem = match{mi,2}(find(match{mi,2}));
    idExpertApp = setdiff(match{mi,3}, idExpertRem);
    expertTouch = tune.expert(mi).touchID;
    expertTuned = expertTouch(find(tune.expert(mi).tuned));
    
    propTuned(mi,3) = length(intersect(idExpertRem, expertTuned)) / length(intersect(idExpertRem, expertTouch));
    propTuned(mi,4) = length(intersect(idExpertApp, expertTuned)) / length(intersect(idExpertApp, expertTouch));

end

figure, hold on
bar(mean(propTuned),'k')
errorbar(mean(propTuned), sem(propTuned), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Proportion')
title('Tuned')

%%
propTuned = zeros(numMice,2); % 1: naive 2: expert
for mi = 1 : numMice
    idNaive = match{mi,1};    
    naiveTouch = tune.naive(learnerInd(mi)).touchID;
    naiveTuned = naiveTouch(find(tune.naive(learnerInd(mi)).tuned));
    
    propTuned(mi,1) = length(intersect(idNaive, naiveTuned)) / length(intersect(idNaive, naiveTouch));
    
    idExpert = match{mi,3};
    expertTouch = tune.expert(mi).touchID;
    expertTuned = expertTouch(find(tune.expert(mi).tuned));
    
    propTuned(mi,2) = length(intersect(idExpert, expertTuned)) / length(intersect(idExpert, expertTouch));

end

%%
figure, hold on
for i = 1 : size(propTuned,1)
    plot(propTuned(i,:), 'ko-')
end
errorbar(mean(propTuned), sem(propTuned), 'ro', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
xlim([0.5 2.5])
ylabel('Proportion (/touch cells)')
title('Tuned cells')

[~,p] = ttest(propTuned(:,1) - propTuned(:,2))
%% Tuned prop increases, even in matched populations.
%% Not so much difference between remained and not.

%% How about modulation, sharpness, and type?

modulations = zeros(numMice,4);
sharpness = zeros(numMice,4);
types = zeros(numMice,4,3); % 3rd dim, 1: sharp, 2: broad, 3: complex
for mi = 1 : numMice
    idNaiveDis = match{mi,1}(find(match{mi,2} == 0));
    idNaiveRem = match{mi,1}(find(match{mi,2}));
    naiveTouch = tune.naive(learnerInd(mi)).touchID;
    naiveTuned = naiveTouch(find(tune.naive(learnerInd(mi)).tuned));
    naiveTunedDis = intersect(idNaiveDis, naiveTuned);
    naiveTunedRem = intersect(idNaiveRem, naiveTuned);
    naiveTunedDisInd = find(ismember(naiveTouch, naiveTunedDis));
    naiveTunedRemInd = find(ismember(naiveTouch, naiveTunedRem));
    
    modulations(mi,1) = mean(tune.naive(learnerInd(mi)).modulation(naiveTunedDisInd));
    modulations(mi,2) = mean(tune.naive(learnerInd(mi)).modulation(naiveTunedRemInd));
    sharpness(mi,1) = mean(tune.naive(learnerInd(mi)).sharpness(naiveTunedDisInd));
    sharpness(mi,2) = mean(tune.naive(learnerInd(mi)).sharpness(naiveTunedRemInd));
    
    types(mi,1,1) = length(find(tune.naive(learnerInd(mi)).unimodalSingle(naiveTunedDisInd))) / length(naiveTunedDisInd);
    types(mi,1,2) = length(find(tune.naive(learnerInd(mi)).unimodalBroad(naiveTunedDisInd))) / length(naiveTunedDisInd);
    types(mi,1,3) = length(find(tune.naive(learnerInd(mi)).multimodal(naiveTunedDisInd))) / length(naiveTunedDisInd);
    
    types(mi,2,1) = length(find(tune.naive(learnerInd(mi)).unimodalSingle(naiveTunedRemInd))) / length(naiveTunedRemInd);
    types(mi,2,2) = length(find(tune.naive(learnerInd(mi)).unimodalBroad(naiveTunedRemInd))) / length(naiveTunedRemInd);
    types(mi,2,3) = length(find(tune.naive(learnerInd(mi)).multimodal(naiveTunedRemInd))) / length(naiveTunedRemInd);
    
    idExpertRem = match{mi,2}(find(match{mi,2}));
    idExpertApp = setdiff(match{mi,3}, idExpertRem);
    expertTouch = tune.expert(mi).touchID;
    expertTuned = expertTouch(find(tune.expert(mi).tuned));
    expertTunedRem = intersect(idExpertRem, expertTuned);
    expertTunedApp = intersect(idExpertApp, expertTuned);    
    expertTunedRemInd = find(ismember(expertTouch, expertTunedRem));
    expertTunedAppInd = find(ismember(expertTouch, expertTunedApp));
    
    modulations(mi,3) = mean(tune.expert(mi).modulation(expertTunedRemInd));
    modulations(mi,4) = mean(tune.expert(mi).modulation(expertTunedAppInd));
    sharpness(mi,3) = mean(tune.expert(mi).sharpness(expertTunedRemInd));
    sharpness(mi,4) = mean(tune.expert(mi).sharpness(expertTunedAppInd));
    
    types(mi,3,1) = length(find(tune.expert(mi).unimodalSingle(expertTunedRemInd))) / length(expertTunedRemInd);
    types(mi,3,2) = length(find(tune.expert(mi).unimodalBroad(expertTunedRemInd))) / length(expertTunedRemInd);
    types(mi,3,3) = length(find(tune.expert(mi).multimodal(expertTunedRemInd))) / length(expertTunedRemInd);

    types(mi,4,1) = length(find(tune.expert(mi).unimodalSingle(expertTunedAppInd))) / length(expertTunedAppInd);
    types(mi,4,2) = length(find(tune.expert(mi).unimodalBroad(expertTunedAppInd))) / length(expertTunedAppInd);
    types(mi,4,3) = length(find(tune.expert(mi).multimodal(expertTunedAppInd))) / length(expertTunedAppInd);
end

%%
figure, 
subplot(131), hold on
bar(mean(modulations),'k')
errorbar(mean(modulations), sem(modulations), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Modulation')

subplot(132), hold on
bar(mean(sharpness),'k')
errorbar(mean(sharpness), sem(sharpness), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Sharpness')

subplot(133), hold on
bar(mean(sum(types, 3)), 'facecolor', ones(1,3)*1)
bar(mean(sum(types(:,:,1:2), 3)), 'facecolor', ones(1,3)*0.7)
bar(mean(squeeze(types(:,:,1))), 'facecolor', ones(1,3)*0)
legend({'Complex', 'Broad', 'Sharp'}, 'box', 'off', 'autoupdate', false)
bar(mean(sum(types, 3)), 'facecolor', ones(1,3)*1)
errorbar(mean(sum(types, 3)), sem(squeeze(types(:,:,3))), 'k', 'lines', 'no')
bar(mean(sum(types(:,:,1:2), 3)), 'facecolor', ones(1,3)*0.7)
errorbar(mean(sum(types(:,:,1:2), 3)), sem(squeeze(types(:,:,2))), 'k', 'lines', 'no')
bar(mean(squeeze(types(:,:,1))), 'facecolor', ones(1,3)*0)
errorbar(mean(squeeze(types(:,:,1))), sem(squeeze(types(:,:,1))), 'k', 'lines', 'no')
xticks([1:4])
xticklabels({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xtickangle(45)
ylabel('Types')


%% Increased modulation in remained, compared to changed populations.
%% Types distributions are similar across populations.


%% how about distribution of tuned angles?

angles = 45:15:135;
angleDist = zeros(numMice, length(angles), 4); 
for mi = 1 : numMice
    idNaiveDis = match{mi,1}(find(match{mi,2} == 0));
    idNaiveRem = match{mi,1}(find(match{mi,2}));
    naiveTouch = tune.naive(learnerInd(mi)).touchID;
    naiveTuned = naiveTouch(find(tune.naive(learnerInd(mi)).tuned));
    naiveTunedDis = intersect(idNaiveDis, naiveTuned);
    naiveTunedRem = intersect(idNaiveRem, naiveTuned);
    naiveTunedDisInd = find(ismember(naiveTouch, naiveTunedDis));
    naiveTunedRemInd = find(ismember(naiveTouch, naiveTunedRem));
    naiveTunedAngles = tune.naive(learnerInd(mi)).tunedAngle;
    
    idExpertRem = match{mi,2}(find(match{mi,2}));
    idExpertApp = setdiff(match{mi,3}, idExpertRem);
    expertTouch = tune.expert(mi).touchID;
    expertTuned = expertTouch(find(tune.expert(mi).tuned));
    expertTunedRem = intersect(idExpertRem, expertTuned);
    expertTunedApp = intersect(idExpertApp, expertTuned);    
    expertTunedRemInd = find(ismember(expertTouch, expertTunedRem));
    expertTunedAppInd = find(ismember(expertTouch, expertTunedApp));
    expertTunedAngles = tune.expert(mi).tunedAngle;
    
    for ai = 1 : length(angles)
        angleDist(mi,ai,1) = length(find(naiveTunedAngles(naiveTunedDisInd) == angles(ai))) / length(naiveTunedDisInd);
        angleDist(mi,ai,2) = length(find(naiveTunedAngles(naiveTunedRemInd) == angles(ai))) / length(naiveTunedRemInd);
        angleDist(mi,ai,3) = length(find(expertTunedAngles(expertTunedRemInd) == angles(ai))) / length(expertTunedRemInd);
        angleDist(mi,ai,4) = length(find(expertTunedAngles(expertTunedAppInd) == angles(ai))) / length(expertTunedAppInd);
    end    
end

colors = {'c', 'b', 'r', 'm'};
figure, hold on
for i = 1 : 4
    errorbar(angles, mean(squeeze(angleDist(:,:,i))), sem(squeeze(angleDist(:,:,i))), colors{i}, 'o')
end
%%
legend({'naive disapp', 'naive remained', 'expert remained', 'expert app'})
xticks(angles)
xlabel('Object angle (\circ)')
ylabel('Proportion')

%% No difference between groups.


%% How does tuned angle change, within matched cells?
%% first, how much of them "change" touch / not-touch responsiveness, and
%% how much of them "change" tuned / not-tuned?
changesTouch = zeros(numMice,4); % 1 touch -> touch, 2 touch -> non-touch, 3 non-touch -> touch, 4 non-touch -> non-touch
for mi = 1 : numMice

    idTouchNaive = tune.naive(learnerInd(mi)).touchID;
    indTouchNaive = find(ismember(match{mi,1}, idTouchNaive));
    indNTNaive = find(1-ismember(match{mi,1}, idTouchNaive));
    refIdMatchedTouchNaive = match{mi,2}(indTouchNaive);
    refIdMatchedNTNaive = match{mi,2}(indNTNaive);
    refIdMatchedTouchNaive = refIdMatchedTouchNaive(find(refIdMatchedTouchNaive));
    refIdMatchedNTNaive = refIdMatchedNTNaive(find(refIdMatchedNTNaive));

    refIdMatchedTouchExpert = intersect(match{mi,2}, tune.expert(mi).touchID);
    refIdMatchedNTExpert = setdiff(match{mi,2}, tune.expert(mi).touchID);

    changesTouch(mi,1) = length(intersect(refIdMatchedTouchExpert, refIdMatchedTouchNaive)) / length(refIdMatchedTouchNaive);
    changesTouch(mi,2) = length(intersect(refIdMatchedNTExpert, refIdMatchedTouchNaive)) / length(refIdMatchedTouchNaive);

    changesTouch(mi,3) = length(intersect(refIdMatchedTouchExpert, refIdMatchedNTNaive)) / length(refIdMatchedNTNaive);
    changesTouch(mi,4) = length(intersect(refIdMatchedNTExpert, refIdMatchedNTNaive)) / length(refIdMatchedNTNaive);
end
%%
figure, hold on
bar(1:2, mean(changesTouch(:,1:2)), 'k')
bar(3:4, mean(changesTouch(:,3:4)), 'w')
errorbar(mean(changesTouch), sem(changesTouch), 'k', 'lines', 'no')
xticks(1:4), xticklabels({'T \rightarrow T', 'T \rightarrow NT', 'NT \rightarrow T', 'NT \rightarrow NT'})
ylabel('Proportion')
legend({'Touch at naive', 'Non-touch at naive'}, 'location', 'northwest', 'box', 'off')

%% where does the reduction come from? layers and C2
% for depth, follow the value in naive
% first, compare proportion between C2 and non-C2, and L2/3 and 4, between
% matched and all

propAll = zeros(numMice,4); % 1 L2/3 C2, 2 L2/3 non-C2, 3 L4 C2, 4 L4 non-C2
propMatched = zeros(numMice,4);
for mi = 1 : numMice

    
    indL23 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths < 350);
    indL4 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths >= 350);
    indC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2 == 1);
    indNonC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2 == 0);
    numCell = length(cellFunction.naive(learnerInd(mi)).allCell.cellNums);
    
    propAll(mi,1) = length(intersect(indL23, indC2)) / numCell;
    propAll(mi,2) = length(intersect(indL23, indNonC2)) / numCell;
    propAll(mi,3) = length(intersect(indL4, indC2)) / numCell;
    propAll(mi,4) = length(intersect(indL4, indNonC2)) / numCell;
    
    indMatched = find(match{mi,2});
    
    indL23 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) < 350);
    indL4 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) >= 350);
    indC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 1);
    indNonC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 0);
    numCell = length(indMatched);
    
    propMatched(mi,1) = length(intersect(indL23, indC2)) / numCell;
    propMatched(mi,2) = length(intersect(indL23, indNonC2)) / numCell;
    propMatched(mi,3) = length(intersect(indL4, indC2)) / numCell;
    propMatched(mi,4) = length(intersect(indL4, indNonC2)) / numCell;
    
end

posAdj = 0.2;
barWidth = 0.4;
figure, hold on
bar([1:4] - posAdj, mean(propAll), barWidth, 'k')
bar([1:4] + posAdj, mean(propMatched), barWidth, 'c')
errorbar([1:4] - posAdj, mean(propAll), sem(propAll), 'k', 'lines', 'no')
errorbar([1:4] + posAdj, mean(propMatched), sem(propMatched), 'k', 'lines', 'no')
legend({'All', 'Matched'}, 'box', 'off')
xticks(1:4)
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})

ylabel('Proportion')

%% Result: there is no change in where the active cells come from.


%% Is there a preferred location of T->NT cells?
% conditional probability of T -> NT, in each location. 
% so, first, need to show prop of T cells in each location.
propTouchAll = zeros(numMice,4); % 1 L2/3 C2, 2 L2/3 non-C2, 3 L4 C2, 4 L4 non-C2
propTouchMatched = zeros(numMice,4);
for mi = 1 : numMice

    
    indL23 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths < 350);
    indL4 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths >= 350);
    indC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2 == 1);
    indNonC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2 == 0);
    indTouch = find(ismember(cellFunction.naive(learnerInd(mi)).allCell.cellNums, cellFunction.naive(learnerInd(mi)).touchID));
    
    propTouchAll(mi,1) = length(intersect(intersect(indL23, indC2), indTouch)) / length(intersect(indL23, indC2));
    propTouchAll(mi,2) = length(intersect(intersect(indL23, indNonC2), indTouch)) / length(intersect(indL23, indNonC2));
    propTouchAll(mi,3) = length(intersect(intersect(indL4, indC2), indTouch)) / length(intersect(indL4, indC2));
    propTouchAll(mi,4) = length(intersect(intersect(indL4, indNonC2), indTouch)) / length(intersect(indL4, indNonC2));
    
    indMatched = find(match{mi,2});
    
    indL23 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) < 350);
    indL4 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) >= 350);
    indC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 1);
    indNonC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 0);
    indTouch = find(ismember(cellFunction.naive(learnerInd(mi)).allCell.cellNums(indMatched), cellFunction.naive(learnerInd(mi)).touchID));
    
    propTouchMatched(mi,1) = length(intersect(intersect(indL23, indC2), indTouch)) / length(intersect(indL23, indC2));
    propTouchMatched(mi,2) = length(intersect(intersect(indL23, indNonC2), indTouch)) / length(intersect(indL23, indNonC2));
    propTouchMatched(mi,3) = length(intersect(intersect(indL4, indC2), indTouch)) / length(intersect(indL4, indC2));
    propTouchMatched(mi,4) = length(intersect(intersect(indL4, indNonC2), indTouch)) / length(intersect(indL4, indNonC2));
end
%%
posAdj = 0.2;
barWidth = 0.4;
figure, hold on
bar([1:4] - posAdj, mean(propTouchAll), barWidth, 'k')
bar([1:4] + posAdj, nanmean(propTouchMatched), barWidth, 'c')
errorbar([1:4] - posAdj, mean(propTouchAll), sem(propTouchAll), 'k', 'lines', 'no')
errorbar([1:4] + posAdj, nanmean(propTouchMatched), sem(propTouchMatched), 'k', 'lines', 'no')
legend({'All', 'Matched'}, 'box', 'off')
xticks(1:4)
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})

ylabel('Proportion of touch cell in each location')

%% Results: distributions and proportions are similar between all cells and matched cells.

%% Where are T->NT cells located?
propT2NT = zeros(numMice,4); % 1 L2/3 C2, 2 L2/3 non-C2, 3 L4 C2, 4 L4 non-C2

for mi = 1 : numMice
    
    indMatched = find(match{mi,2});
    
    indL23 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) < 350);
    indL4 = find(cellFunction.naive(learnerInd(mi)).allCell.cellDepths(indMatched) >= 350);
    indC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 1);
    indNonC2 = find(cellFunction.naive(learnerInd(mi)).allCell.isC2(indMatched) == 0);
    indTouch = find(ismember(cellFunction.naive(learnerInd(mi)).allCell.cellNums(indMatched), cellFunction.naive(learnerInd(mi)).touchID));
    
    refIdNT = intersect(match{mi,3}, tune.expert(mi).touchID);
    indFutureNT = find(ismember(match{mi,2}, refIdNT));
    
    propT2NT(mi,1) = length(intersect(intersect(indL23, indC2), intersect(indTouch, indFutureNT))) / length(intersect(intersect(indL23, indC2), indTouch));
    propT2NT(mi,2) = length(intersect(intersect(indL23, indNonC2), intersect(indTouch, indFutureNT))) / length(intersect(intersect(indL23, indNonC2), indTouch));
    propT2NT(mi,3) = length(intersect(intersect(indL4, indC2), intersect(indTouch, indFutureNT))) / length(intersect(intersect(indL4, indC2), indTouch));
    propT2NT(mi,4) = length(intersect(intersect(indL4, indNonC2), intersect(indTouch, indFutureNT))) / length(intersect(intersect(indL4, indNonC2), indTouch));
    
end

figure, hold on
bar(nanmean(propT2NT), 'k')
errorbar(nanmean(propT2NT), sem(propT2NT), 'k', 'lines', 'no')

xticks(1:4)
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
ylabel('Proportion T \rightarrow NT')


%% Do the same thing for tuned cells.
changesTuned = zeros(numMice,6); % 1 tuned -> tuned, 2 tuned -> not-tuned, 3 tuned -> non-touch, 4 not-tuned -> tuned, 5 not-tuned -> not-tuned, 6 not-tuned -> non-touch
% because very few non-touch cells become touch cells, these are excluded from the analysis
for mi = 1 : numMice
    idTouchNaive = tune.naive(learnerInd(mi)).touchID;
    idTunedNaive = idTouchNaive(find(tune.naive(learnerInd(mi)).tuned == 1));
    idNTNaive = idTouchNaive(find(tune.naive(learnerInd(mi)).tuned == 0));
    indTunedNaive = find(ismember(match{mi,1}, idTunedNaive));
    indNTNaive = find(ismember(match{mi,1}, idNTNaive));
    refIdMatchedTunedNaive = match{mi,2}(indTunedNaive);
    refIdMatchedNTNaive = match{mi,2}(indNTNaive);
    refIdMatchedTunedNaive = refIdMatchedTunedNaive(find(refIdMatchedTunedNaive));
    refIdMatchedNTNaive = refIdMatchedNTNaive(find(refIdMatchedNTNaive));

    refIdMatchedTouchExpert = intersect(match{mi,2}, tune.expert(mi).touchID);
    refIdmatchedNonTouchExpert = setdiff(match{mi,2}, tune.expert(mi).touchID);
    refIdMatchedTunedExpert = intersect(match{mi,2}, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    refIdMatchedNTExpert = setdiff(refIdMatchedTouchExpert, refIdMatchedTunedExpert);

    changesTuned(mi,1) = length(intersect(refIdMatchedTunedExpert, refIdMatchedTunedNaive)) / length(refIdMatchedTunedNaive);
    changesTuned(mi,2) = length(intersect(refIdMatchedNTExpert, refIdMatchedTunedNaive)) / length(refIdMatchedTunedNaive);
    changesTuned(mi,3) = length(intersect(refIdmatchedNonTouchExpert, refIdMatchedTunedNaive)) / length(refIdMatchedTunedNaive);

    changesTuned(mi,4) = length(intersect(refIdMatchedTunedExpert, refIdMatchedNTNaive)) / length(refIdMatchedNTNaive);
    changesTuned(mi,5) = length(intersect(refIdMatchedNTExpert, refIdMatchedNTNaive)) / length(refIdMatchedNTNaive);
    changesTuned(mi,6) = length(intersect(refIdmatchedNonTouchExpert, refIdMatchedNTNaive)) / length(refIdMatchedNTNaive);
end
%%
figure, hold on
bar(1:3, mean(changesTuned(:,1:3)), 'k')
bar(4:6, mean(changesTuned(:,4:6)), 'w')
errorbar(mean(changesTuned), sem(changesTuned), 'k', 'lines', 'no')
xticks(1:6), xticklabels({'Tuned \rightarrow Tuned', 'Tuned \rightarrow Not-tuned', 'Tuned \rightarrow Non-touch', ...
    'Not-tuned \rightarrow Tuned', 'Not-tuned \rightarrow Not-tuned', 'Not-tuned \rightarrow Non-touch'})
ylabel('Proportion')
legend({'Tuned at naive', 'Not-tuned at naive'}, 'location', 'northeast', 'box', 'off')
%%
xtickangle(45)



%% Finally, how much angle do they change on average?
% only among tuned -> tuned

taChange = cell(numMice,1); % tuned-angle Change
for mi = 1 : numMice
    indTunedNaive = find(ismember( match{mi,1}, tune.naive(learnerInd(mi)).touchID(find(tune.naive(learnerInd(mi)).tuned)) ));
    indMatchedNaive = find(match{mi,2});
    refIdTunedNaive = match{mi,2}(intersect(indTunedNaive, indMatchedNaive));
    refIdTunedExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    idExpertBothTuned = intersect(refIdTunedNaive, refIdTunedExpert); % sorted in order
    
    idNaiveBothTunedSorted = zeros(length(idExpertBothTuned),1);
    tunedAngleNaive = zeros(length(idExpertBothTuned),1);
    tunedAngleExpert = zeros(length(idExpertBothTuned),1);
    for idi = 1 : length(idExpertBothTuned)
        idNaiveBothTunedSorted(idi) = match{mi,1}(find(match{mi,2} == idExpertBothTuned(idi)));
        tunedAngleNaive(idi) = tune.naive(learnerInd(mi)).tunedAngle( find( tune.naive(learnerInd(mi)).touchID == idNaiveBothTunedSorted(idi) ) );
        tunedAngleExpert(idi) = tune.expert(mi).tunedAngle( find ( tune.expert(mi).touchID == idExpertBothTuned(idi) ) );
    end
    
    taChange{mi} = tunedAngleExpert - tunedAngleNaive;
    
end

%%
histRange = 0:15:105;
histChange = zeros(numMice, length(histRange)-1);
for i = 1 : numMice
    histChange(i,:) = histcounts(abs(taChange{i}), histRange, 'norma', 'prob');
end
figure, hold on
bar(histRange(1:end-1), mean(histChange), 'k')
errorbar(histRange(1:end-1), mean(histChange), sem(histChange), 'k', 'lines', 'no')
xticks(histRange(1:end-1))
xlabel('\DeltaAngle (\circ)')
ylabel('Proportion')


%% Results: More than half remained their tuned angle, > 80% remained within 1 angle bin

%% How about within sharp tuned at naive cells?

taChange = cell(numMice,1); % tuned-angle Change
for mi = 1 : numMice
    indSharpTunedNaive = find(ismember( match{mi,1}, tune.naive(learnerInd(mi)).touchID(find(tune.naive(learnerInd(mi)).unimodalSingle)) ));
    indMatchedNaive = find(match{mi,2});
    refIdTunedNaive = match{mi,2}(intersect(indSharpTunedNaive, indMatchedNaive));
    refIdTunedExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    idExpertBothTuned = intersect(refIdTunedNaive, refIdTunedExpert); % sorted in order
    
    idNaiveBothTunedSorted = zeros(length(idExpertBothTuned),1);
    tunedAngleNaive = zeros(length(idExpertBothTuned),1);
    tunedAngleExpert = zeros(length(idExpertBothTuned),1);
    for idi = 1 : length(idExpertBothTuned)
        idNaiveBothTunedSorted(idi) = match{mi,1}(find(match{mi,2} == idExpertBothTuned(idi)));
        tunedAngleNaive(idi) = tune.naive(learnerInd(mi)).tunedAngle( find( tune.naive(learnerInd(mi)).touchID == idNaiveBothTunedSorted(idi) ) );
        tunedAngleExpert(idi) = tune.expert(mi).tunedAngle( find ( tune.expert(mi).touchID == idExpertBothTuned(idi) ) );
    end
    
    taChange{mi} = tunedAngleExpert - tunedAngleNaive;
    
end


histRange = 0:15:105;
histChange = zeros(numMice, length(histRange)-1);
for i = 1 : numMice
    histChange(i,:) = histcounts(abs(taChange{i}), histRange, 'norma', 'prob');
end
figure, hold on
bar(histRange(1:end-1), mean(histChange), 'k')
errorbar(histRange(1:end-1), mean(histChange), sem(histChange), 'k', 'lines', 'no')
xticks(histRange(1:end-1))
xlabel('\DeltaAngle (\circ)')
ylabel('Proportion')
title('Sharp tuned cells in naive')

%% Results: More than 80% remained their tuned angle, > 95% remained within 1 angle bin


%% How about tuned features? (from whisker model)
% only the top features
% only the top features w/o touch counts
% including sub features

loadfn = 'glmResults_WKV_touchCell_exclusion_NC';
whisker = load(sprintf('%s%s',baseDir,loadfn), 'naive', 'expert');
featureNames = {'maxDq', 'maxDf', 'maxDkH', 'maxDkV', 'max(slide distance)', 'max(duration)', ...    
                            'q', 'f', 'kH', 'kV', 'arc length', 'touch count'};

fitCountsTop = zeros(length(featureNames), length(featureNames), numMice);
fitCountsTopWoTC = zeros(length(featureNames), length(featureNames), numMice);
fitCountsIncSub = zeros(length(featureNames), length(featureNames), numMice);

for mi = 1 : numMice
    idMatchNaive = match{mi,1}(find(match{mi,2}));
    idTuneNaive = tune.naive(learnerInd(mi)).touchID(find(tune.naive(learnerInd(mi)).tuned));
    idMatchTuneNaive = intersect(idMatchNaive, idTuneNaive);
    idRefMatchTuneNaive = match{mi,2}(find(ismember(match{mi,1}, idMatchTuneNaive)));

    idMatchExpert = match{mi,2};
    idTuneExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    idMatchTuneExpert = intersect(idMatchExpert, idTuneExpert);

    idRefBothTune = intersect(idRefMatchTuneNaive, idMatchTuneExpert);
    idNaiveBothTune = match{mi,1}(find(ismember(match{mi,2}, idRefBothTune)));
    idExpertBothTune = zeros(length(idNaiveBothTune),1);
    for i = 1 : length(idNaiveBothTune)
        idExpertBothTune(i) = match{mi,2}(find(match{mi,1} == idNaiveBothTune(i)));
    end

    for i = 1 : length(idNaiveBothTune)
        whiskerIndNaive = find(whisker.naive(learnerInd(mi)).cID == idNaiveBothTune(i));
        viNaive = whisker.naive(learnerInd(mi)).whiskerVariableDEdiff(whiskerIndNaive,:);
        [~, maxFeatIndNaive] = max(viNaive);

        whiskerIndExpert = find(whisker.expert(mi).cID == idExpertBothTune(i));
        viExpert = whisker.expert(mi).whiskerVariableDEdiff(whiskerIndExpert,:);
        [~, maxFeatIndExpert] = max(viExpert);

        fitCountsTop(maxFeatIndNaive, maxFeatIndExpert, mi) = fitCountsTop(maxFeatIndNaive, maxFeatIndExpert, mi) + 1;
    end
end

%%
fitRateTop = fitCountsTop ./ sum(fitCountsTop,2);
fitRateTopMouse = fitCountsTop ./ sum(sum(fitCountsTop,2));

% figure, imagesc(mean(fitCountsTop,3)), colorbar, axis square
figure, imagesc(mean(fitRateTopMouse,3)), colorbar, axis square
xticks([1:12])
xticklabels(featureNames)
xtickangle(45)
xlabel('Expert')
yticks([1:12])
yticklabels(featureNames)
ylabel('Naive')
title('Top whisker features')
% figure, imagesc(nanmean(fitRateTop,3)), colorbar, axis square
%%
trace(mean(fitRateTopMouse,3))
sum(sum(mean(fitRateTopMouse,3))) - trace(mean(fitRateTopMouse,3))
%%
trace(sum(fitCountsTop,3))
sum(sum(sum(fitCountsTop,3))) - trace(sum(fitCountsTop,3))

%% Results: Not so much of correct prediction from naive to expert

%% How does the distribution look like, compared between naive and expert?
% from tuned cells
% from both tuned cells
% from top features and including sub features

topFeatDist = zeros(numMice, length(featureNames), 2); % 1 for naive, 2 for expert
for mi = 1 : numMice
    tunedIdNaive = tune.naive(learnerInd(mi)).touchID(find(tune.naive(learnerInd(mi)).tuned));
    tunedWhiskerIndNaive = find(ismember(whisker.naive(learnerInd(mi)).cID, tunedIdNaive));
    [~, distTemp] = max(whisker.naive(learnerInd(mi)).whiskerVariableDEdiff(tunedWhiskerIndNaive,:), [], 2);
    topFeatDist(mi,:,1) = histcounts(distTemp, 1:13, 'norm', 'prob');
    tunedIdExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    tunedWhiskerIndExpert = find(ismember(whisker.expert(mi).cID, tunedIdExpert));
    [~, distTemp] = max(whisker.expert(mi).whiskerVariableDEdiff(tunedWhiskerIndExpert,:), [], 2);
    topFeatDist(mi,:,2) = histcounts(distTemp, 1:13, 'norm', 'prob');
end

%%
posAdj = 0.2;
binWidth = 0.4;
figure, hold on
naiveMat = squeeze(topFeatDist(:,:,1));
expertMat = squeeze(topFeatDist(:,:,2));
bar([1:12] - posAdj, mean(naiveMat), binWidth, 'k')
bar([1:12] + posAdj, mean(expertMat), binWidth, 'c')
errorbar([1:12] - posAdj, mean(naiveMat), zeros(1,length(featureNames)), sem(naiveMat), 'k', 'lines', 'no')
errorbar([1:12] + posAdj, mean(expertMat), zeros(1,length(featureNames)), sem(expertMat), 'k', 'lines', 'no')
xticks([1:12])
xticklabels(featureNames)
xtickangle(45)
legend({'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')

%% results: quite complicated... top features change a little.

% from both tuned cells
topFeatDistBoth = zeros(numMice, length(featureNames), 2); % 1 for naive, 2 for expert
for mi = 1 : numMice
    idMatchNaive = match{mi,1}(find(match{mi,2}));
    idTuneNaive = tune.naive(learnerInd(mi)).touchID(find(tune.naive(learnerInd(mi)).tuned));
    idMatchTuneNaive = intersect(idMatchNaive, idTuneNaive);
    idRefMatchTuneNaive = match{mi,2}(find(ismember(match{mi,1}, idMatchTuneNaive)));

    idMatchExpert = match{mi,2};
    idTuneExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    idMatchTuneExpert = intersect(idMatchExpert, idTuneExpert);

    idRefBothTune = intersect(idRefMatchTuneNaive, idMatchTuneExpert);
    idNaiveBothTune = match{mi,1}(find(ismember(match{mi,2}, idRefBothTune)));
    idExpertBothTune = zeros(length(idNaiveBothTune),1);
    for i = 1 : length(idNaiveBothTune)
        idExpertBothTune(i) = match{mi,2}(find(match{mi,1} == idNaiveBothTune(i)));
    end

    indNaiveBothTune = find(ismember(whisker.naive(learnerInd(mi)).cID, idNaiveBothTune));
    [~, featIndTemp] = max(whisker.naive(learnerInd(mi)).whiskerVariableDEdiff(indNaiveBothTune,:), [], 2);
    topFeatDistBoth(mi,:,1) = histcounts(featIndTemp, 1:13, 'norm', 'prob');
    
    indExpertBothTune = find(ismember(whisker.expert(mi).cID, idExpertBothTune));
    [~, featIndTemp] = max(whisker.expert(mi).whiskerVariableDEdiff(indExpertBothTune,:), [], 2);
    topFeatDistBoth(mi,:,2) = histcounts(featIndTemp, 1:13, 'norm', 'prob');
    
end

posAdj = 0.2;
binWidth = 0.4;
figure, hold on
naiveMat = squeeze(topFeatDistBoth(:,:,1));
expertMat = squeeze(topFeatDistBoth(:,:,2));
bar([1:12] - posAdj, mean(naiveMat), binWidth, 'k')
bar([1:12] + posAdj, mean(expertMat), binWidth, 'c')
errorbar([1:12] - posAdj, mean(naiveMat), zeros(1,length(featureNames)), sem(naiveMat), 'k', 'lines', 'no')
errorbar([1:12] + posAdj, mean(expertMat), zeros(1,length(featureNames)), sem(expertMat), 'k', 'lines', 'no')
xticks([1:12])
xticklabels(featureNames)
xtickangle(45)
legend({'Naive', 'Expert'}, 'box', 'off', 'location', 'northwest')
%%
title('From both tuned cells')
%% How about most important features in angle tuning?
% correlation methods
% (1) from all tuned cells of both naive and expert
% (2) from both tuned cells of naive and expert
%%
loadFn1 = 'modelAngleTuning_NC';
loadFn2 = 'modelAngleTuning_NC_combinations';
data1 = load(sprintf('%s%s',baseDir, loadFn1), 'naive', 'expert');
data2 = load(sprintf('%s%s',baseDir, loadFn2), 'naive', 'expert');

%% (1)
% touch, full, 12 drop-outs, others only (16 total)
corrWhisker = cell(numMice,2); % 1 for naive, 2 for expert
corrFeaturesOut = cell(numMice,12,2);
corrOther = cell(numMice,2);
corrFeaturesCombOut = cell(numMice,4,2);
corrFeaturesCombIn = cell(numMice,4,2);
corrFeaturesIn = cell(numMice,12,2);
for mi = 1 : numMice
    % naive
    indTuned = find(data1.naive(learnerInd(mi)).tunedAllCell(:,1));
    indTemp = find(data1.naive(learnerInd(mi)).tunedAllCell(indTuned,3));
    indWhiskerNaive = indTuned(indTemp);
    corrWhisker{mi,1} = zeros(length(indWhiskerNaive),1);
    
    %expert
    indTuned = find(data1.expert(mi).tunedAllCell(:,1));
    indTemp = find(data1.expert(mi).tunedAllCell(indTuned,3));
    indWhiskerExpert = indTuned(indTemp);
    corrWhisker{mi,2} = zeros(length(indWhiskerExpert),1);
    
    % naive
    for j = 1 : length(indWhiskerNaive)
        tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),3}));
        if tempVal < 0 
            tempVal = 0;
        end
        corrWhisker{mi,1}(j) = tempVal;
        if data2.naive(learnerInd(mi)).tunedAllCell(indWhiskerNaive(j),1) == 1
            tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data2.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j), 1}));
            if isnan(tempVal) || tempVal < 0
                tempVal = 0;
            end
        else
            tempVal = 0;
        end
        corrOther{mi,1}(j) = tempVal;
    end
    % expert
    for j = 1 : length(indWhiskerExpert)
        tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),3}));
        if tempVal < 0 
            tempVal = 0;
        end
        corrWhisker{mi,2}(j) = tempVal;
        if data2.expert(mi).tunedAllCell(indWhiskerExpert(j),1) == 1
            tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data2.expert(mi).spkValAllCell{indWhiskerExpert(j), 1}));
            if isnan(tempVal) || tempVal < 0
                tempVal = 0;
            end
        else
            tempVal = 0;
        end
        corrOther{mi,2}(j) = tempVal;
    end
    
    % features out
    for fi = 1 : 12
        % naive
        corrFeaturesOut{mi,fi,1} = zeros(length(indWhiskerNaive),1);        
        for j = 1 : length(indWhiskerNaive)
            if data1.naive(learnerInd(mi)).tunedAllCell(indWhiskerNaive(j),3+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),3+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesOut{mi,fi,1}(j) = tempVal;
        end
        % expert
        corrFeaturesOut{mi,fi,2} = zeros(length(indWhiskerExpert),1);        
        for j = 1 : length(indWhiskerExpert)
            if data1.expert(mi).tunedAllCell(indWhiskerExpert(j),3+fi) == 1
                tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),3+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesOut{mi,fi,2}(j) = tempVal;
        end
    end
    
    % combination of features
    for fi = 1 : 4
        % naive
        corrFeaturesCombOut{mi,fi,1} = zeros(length(indWhiskerNaive),1);
        corrFeaturesCombIn{mi,fi,1} = zeros(length(indWhiskerNaive),1);
        for j = 1 : length(indWhiskerNaive)
            if data2.naive(learnerInd(mi)).tunedAllCell(indWhiskerNaive(j),13+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data2.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),13+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombOut{mi,fi,1}(j) = tempVal;
            
            if data2.naive(learnerInd(mi)).tunedAllCell(indWhiskerNaive(j), 21+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data2.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),21+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombIn{mi,fi,1}(j) = tempVal;
        end
        % expert
        corrFeaturesCombOut{mi,fi,2} = zeros(length(indWhiskerExpert),1);
        corrFeaturesCombIn{mi,fi,2} = zeros(length(indWhiskerExpert),1);
        for j = 1 : length(indWhiskerExpert)
            if data2.expert(mi).tunedAllCell(indWhiskerExpert(j),13+fi) == 1
                tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data2.expert(mi).spkValAllCell{indWhiskerExpert(j),13+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombOut{mi,fi,2}(j) = tempVal;
            
            if data2.expert(mi).tunedAllCell(indWhiskerExpert(j), 21+fi) == 1
                tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data2.expert(mi).spkValAllCell{indWhiskerExpert(j),21+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesCombIn{mi,fi,2}(j) = tempVal;
        end
    end
    
    % features in    
    for fi = 1 : 12
        % naive
        corrFeaturesIn{mi,fi,1} = zeros(length(indWhiskerNaive),1);
        for j = 1 : length(indWhiskerNaive)
            if data2.naive(learnerInd(mi)).tunedAllCell(indWhiskerNaive(j),1+fi) == 1
                tempVal = corr(cellfun(@mean, data1.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1}), cellfun(@mean, data2.naive(learnerInd(mi)).spkValAllCell{indWhiskerNaive(j),1+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesIn{mi,fi,1}(j) = tempVal;
        end
        % expert
        corrFeaturesIn{mi,fi,2} = zeros(length(indWhiskerExpert),1);
        for j = 1 : length(indWhiskerExpert)
            if data2.expert(mi).tunedAllCell(indWhiskerExpert(j),1+fi) == 1
                tempVal = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhiskerExpert(j),1}), cellfun(@mean, data2.expert(mi).spkValAllCell{indWhiskerExpert(j),1+fi}));
                if isnan(tempVal) || tempVal < 0
                    tempVal = 0;
                end
            else
                tempVal = 0;
            end
            corrFeaturesIn{mi,fi,2}(j) = tempVal;
        end
    end
end

%%




%%
% corrWhisker = cell(numMice,2); % 1 for naive, 2 for expert
% corrFeaturesOut = cell(numMice,12,2);
% corrOther = cell(numMice,2);
% corrFeaturesCombOut = cell(numMice,4,2);
% corrFeaturesCombIn = cell(numMice,4,2);
% corrFeaturesIn = cell(numMice,12,2);
posAdj = 0.2;
barWidth = 0.4;
figure, hold on
tempWhisker = cellfun(@mean, corrWhisker(:,1));
temptempFeatureOutMean = cellfun(@mean, corrFeaturesOut(:,:,1));
tempFeatureOutMean = mean(temptempFeatureOutMean,2);
tempFeatureOut = cellfun(@mean, corrFeaturesOut(:, [5,4],1));
tempFeatureCombOut = cellfun(@mean, corrFeaturesCombOut(:,3,1));
temptempFeatureInMean = cellfun(@mean, corrFeaturesIn(:,:,1));
tempFeatureInMean = mean(temptempFeatureInMean,2);
tempFeatureIn = cellfun(@mean, corrFeaturesIn(:, [5,4],1));
tempFeatureCombIn = cellfun(@mean, corrFeaturesCombIn(:,3,1));
tempOther = cellfun(@mean, corrOther(:,1));

bar(1-posAdj, mean(tempWhisker), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(1-posAdj, mean(tempWhisker), sem(tempWhisker), 'color', ones(1,3)*0.7)
bar(2-posAdj, mean(tempFeatureOutMean), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(2-posAdj, mean(tempFeatureOutMean), sem(tempFeatureOutMean), 'color', ones(1,3)*0.7)
bar([3:4]-posAdj, mean(tempFeatureOut), barWidth, 'k')
errorbar([3:4]-posAdj, mean(tempFeatureOut), sem(tempFeatureOut), 'k', 'linestyle', 'none')
bar(5-posAdj, mean(tempFeatureCombOut), barWidth, 'k')
errorbar(5-posAdj, mean(tempFeatureCombOut), sem(tempFeatureCombOut), 'k', 'linestyle', 'none')

bar(6-posAdj, mean(tempOther), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(6-posAdj, mean(tempOther), sem(tempOther), 'color', ones(1,3)*0.7)
bar(7-posAdj, mean(tempFeatureInMean), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(7-posAdj, mean(tempFeatureInMean), sem(tempFeatureInMean), 'color', ones(1,3)*0.7)

bar([8:9]-posAdj, mean(tempFeatureIn), barWidth, 'k')
errorbar([8:9]-posAdj, mean(tempFeatureIn), sem(tempFeatureIn), 'k', 'linestyle', 'none')
bar(10-posAdj, mean(tempFeatureCombIn), barWidth, 'k')
errorbar(10-posAdj, mean(tempFeatureCombIn), sem(tempFeatureCombIn), 'k', 'linestyle', 'none')

tempWhisker = cellfun(@mean, corrWhisker(:,2));
temptempFeatureOutMean = cellfun(@mean, corrFeaturesOut(:,:,2));
tempFeatureOutMean = mean(temptempFeatureOutMean,2);
tempFeatureOut = cellfun(@mean, corrFeaturesOut(:, [5,4],2));
tempFeatureCombOut = cellfun(@mean, corrFeaturesCombOut(:,3,2));
temptempFeatureInMean = cellfun(@mean, corrFeaturesIn(:,:,2));
tempFeatureInMean = mean(temptempFeatureInMean,2);
tempFeatureIn = cellfun(@mean, corrFeaturesIn(:, [5,4],2));
tempFeatureCombIn = cellfun(@mean, corrFeaturesCombIn(:,3,2));
tempOther = cellfun(@mean, corrOther(:,2));

bar(1+posAdj, mean(tempWhisker), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(1+posAdj, mean(tempWhisker), zeros(size(mean(tempWhisker))), sem(tempWhisker), 'color', ones(1,3)*0.7)
bar(2+posAdj, mean(tempFeatureOutMean), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(2+posAdj, mean(tempFeatureOutMean), zeros(size(mean(tempFeatureOutMean))), sem(tempFeatureOutMean), 'color', ones(1,3)*0.7)
bar([3:4]+posAdj, mean(tempFeatureOut), barWidth, 'c')
errorbar([3:4]+posAdj, mean(tempFeatureOut), zeros(size(mean(tempFeatureOut))), sem(tempFeatureOut), 'k', 'linestyle', 'none')
bar(5+posAdj, mean(tempFeatureCombOut), barWidth, 'c')
errorbar(5+posAdj, mean(tempFeatureCombOut), zeros(size(mean(tempFeatureCombOut))), sem(tempFeatureCombOut), 'k', 'linestyle', 'none')

bar(6+posAdj, mean(tempOther), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(6+posAdj, mean(tempOther), zeros(size(mean(tempOther))), sem(tempOther), 'color', ones(1,3)*0.7)
bar(7+posAdj, mean(tempFeatureInMean), barWidth, 'facecolor', ones(1,3)*0.7)
errorbar(7+posAdj, mean(tempFeatureInMean), zeros(size(mean(tempFeatureInMean))), sem(tempFeatureInMean), 'color', ones(1,3)*0.7)

bar([8:9]+posAdj, mean(tempFeatureIn), barWidth, 'c')
errorbar([8:9]+posAdj, mean(tempFeatureIn), zeros(size(mean(tempFeatureIn))), sem(tempFeatureIn), 'k', 'linestyle', 'none')
bar(10+posAdj, mean(tempFeatureCombIn), barWidth, 'c')
errorbar(10+posAdj, mean(tempFeatureCombIn), zeros(size(mean(tempFeatureCombIn))), sem(tempFeatureCombIn), 'k', 'linestyle', 'none')



xticks([1:10]), xticklabels({'Full whisker', 'mean drop-out', '-max(SD)', '-maxDkV', '-(maxDkV & max(SD))', ...
    'others only', ...
    'mean drop-in + others', 'max(SD) + others', 'maxDkV + others', '(maxDkV & max(SD)) + others'})
xtickangle(45)
ylim([0 1])
ylabel('Correlation')
