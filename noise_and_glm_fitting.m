%% if noise level affects L4 signals and fittings 

%% plot DE vs noise
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]}; 

cd(baseDir)
load('glm_results_responseType')
naiveInds = 1:12;
figure, hold on
for ni = naiveInds
    mouse = mice(ni);
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{ni}(1);
    ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
    load(ufn)
    plot(u.noise, naive(ni).allDE, 'k.')
end
xlabel('\DeltaF/F_0 noise')
ylabel('Deviance explained')

%%
%% Result seems likely glm fitting is not related to dF noise level
%%


%% Comparing within similar ?F/F0 noise levels in L2/3 and L4
%% first, figure out how to choose the "similar" level of ?F/F0 noise 

figure, 
plot(u.cellDepths, u.noise, 'k.')
upperInd = find(u.cellNums < 5000);
lowerInd = find(u.cellNums > 5000);

upper98th = prctile(u.noise(upperInd),98);
lower2th = prctile(u.noise(lowerInd),2);
withinInd = find(u.noise > lower2th & u.noise < upper98th);
hold on,
plot(u.cellDepths(withinInd), u.noise(withinInd), 'r.')


%% collated distribution of noise level

upperNoiseDist = cell(1,length(naive));
lowerNoiseDist = cell(1,length(naive));
range = 0.27:0.01:0.42;
histUpper = zeros(length(naive),length(range)-1); 
histLower = zeros(length(naive),length(range)-1);
histDiff = zeros(length(naive),length(range)-1);
for i = 1 : length(naive)
    
    mouse = mice(i);
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{i}(1);
    ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
    load(ufn)
    
    upperInd = find(u.cellNums < 5000);
    lowerInd = find(u.cellNums > 5000);

    upper98th = prctile(u.noise(upperInd),98);
    lower2th = prctile(u.noise(lowerInd),2);
    withinInd = find(u.noise > lower2th & u.noise < upper98th);
    
    indL23 = find(u.cellDepths < 350);
    indL4 = find(u.cellDepths >= 350);
    
    upperNoiseDist{i} = u.noise(intersect(indL23, withinInd));
    lowerNoiseDist{i} = u.noise(intersect(indL4, withinInd));
    
    histUpper(i,:) = histcounts(upperNoiseDist{i},range, 'normalization', 'probability');
    histLower(i,:) = histcounts(lowerNoiseDist{i},range, 'normalization', 'probability');
    histDiff(i,:) = histUpper(i,:) - histLower(i,:);
end
wholeHistUpper = histcounts(cell2mat(upperNoiseDist),range, 'normalization', 'probability');
wholeHistLower = histcounts(cell2mat(lowerNoiseDist),range, 'normalization', 'probability');
%
figure, 
subplot(121), plot(range(2:end), wholeHistUpper, 'k-'), hold on, plot(range(2:end), wholeHistLower, 'r-')
subplot(122), hold on
for i = 1 : length(naive)
    plot(range(2:end), histDiff(i,:), '-', 'color', [0.7 0.7 0.7])
end
plot(range(2:end), mean(histUpper - histLower), 'k-')

%% draw cell functions from noise-matched samples
cd(baseDir)
load('cellFunctionLasso_NC')
touches = zeros(length(naive), 4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskings = zeros(length(naive), 4); 
mixed = zeros(length(naive), 4);
others = zeros(length(naive), 4);
for i = 1 : length(naive)
    
    mouse = mice(i);
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{i}(1);
    ufn = sprintf('UberJK%03dS%02d_NC', mouse, session);
    load(ufn)
    
    upperInd = find(u.cellNums < 5000);
    lowerInd = find(u.cellNums > 5000);

    upper98th = prctile(u.noise(upperInd),98);
    lower2th = prctile(u.noise(lowerInd),2);
    withinInd = find(u.noise > lower2th & u.noise < upper98th);

    temp = naive(i);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    
    indTouch = find(ismember(temp.cellNums, temp.touchID));
    indWhisking = find(ismember(temp.cellNums, temp.whiskingID));
    indMixed = intersect(indTouch, indWhisking);
    indOther = find(ismember(temp.cellNums, union(temp.soundID, union(temp.rewardID, temp.lickingID))));
    
    touches(i,1) = length(intersect(intersect(indTouch,withinInd), intersect(indL23, indC2))) / length(intersect(withinInd, intersect(indL23, indC2)));
    touches(i,2) = length(intersect(intersect(indTouch,withinInd), intersect(indL23, indnonC2))) / length(intersect(withinInd, intersect(indL23, indnonC2)));
    touches(i,3) = length(intersect(intersect(indTouch,withinInd), intersect(indL4, indC2))) / length(intersect(withinInd, intersect(indL4, indC2)));
    touches(i,4) = length(intersect(intersect(indTouch,withinInd), intersect(indL4, indnonC2))) / length(intersect(withinInd, intersect(indL4, indnonC2)));
    
    whiskings(i,1) = length(intersect(intersect(indWhisking,withinInd), intersect(indL23, indC2))) / length(intersect(withinInd, intersect(indL23, indC2)));
    whiskings(i,2) = length(intersect(intersect(indWhisking,withinInd), intersect(indL23, indnonC2))) / length(intersect(withinInd, intersect(indL23, indnonC2)));
    whiskings(i,3) = length(intersect(intersect(indWhisking,withinInd), intersect(indL4, indC2))) / length(intersect(withinInd, intersect(indL4, indC2)));
    whiskings(i,4) = length(intersect(intersect(indWhisking,withinInd), intersect(indL4, indnonC2))) / length(intersect(withinInd, intersect(indL4, indnonC2)));
    
    mixed(i,1) = length(intersect(intersect(indMixed,withinInd), intersect(indL23, indC2))) / length(intersect(withinInd, intersect(indL23, indC2)));
    mixed(i,2) = length(intersect(intersect(indMixed,withinInd), intersect(indL23, indnonC2))) / length(intersect(withinInd, intersect(indL23, indnonC2)));
    mixed(i,3) = length(intersect(intersect(indMixed,withinInd), intersect(indL4, indC2))) / length(intersect(withinInd, intersect(indL4, indC2)));
    mixed(i,4) = length(intersect(intersect(indMixed,withinInd), intersect(indL4, indnonC2))) / length(intersect(withinInd, intersect(indL4, indnonC2)));
    
    others(i,1) = length(intersect(intersect(indOther,withinInd), intersect(indL23, indC2))) / length(intersect(withinInd, intersect(indL23, indC2)));
    others(i,2) = length(intersect(intersect(indOther,withinInd), intersect(indL23, indnonC2))) / length(intersect(withinInd, intersect(indL23, indnonC2)));
    others(i,3) = length(intersect(intersect(indOther,withinInd), intersect(indL4, indC2))) / length(intersect(withinInd, intersect(indL4, indC2)));
    others(i,4) = length(intersect(intersect(indOther,withinInd), intersect(indL4, indnonC2))) / length(intersect(withinInd, intersect(indL4, indnonC2)));
end
%
figure, hold all
% for legend ordering
i = 1;
touchOnly = nanmean(touches(:,i)) - nanmean(mixed(:,i));
touchNwhisking = nanmean(mixed(:,i));
whiskingOnly = nanmean(whiskings(:,i)) - nanmean(mixed(:,i));
other = nanmean(others(:,i));
bar(i, touchOnly, 'b')
bar(i, touchOnly + touchNwhisking, 'y')
bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    
for i = 1 : 4
%     subplot(2,2,i) % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
    touchOnly = nanmean(touches(:,i)) - nanmean(mixed(:,i));
    touchNwhisking = nanmean(mixed(:,i));
    whiskingOnly = nanmean(whiskings(:,i)) - nanmean(mixed(:,i));
    other = nanmean(others(:,i));
%     pie([touchOnly, touchNwhisking, whiskingOnly, other, 1-touchOnly - touchNwhisking - whiskingOnly - other])
    bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
    bar(i, touchOnly + touchNwhisking, 'y')
    bar(i, touchOnly, 'b')
end
legend({'Touch', 'Touch & Whisking', 'Whisking', 'Others'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)