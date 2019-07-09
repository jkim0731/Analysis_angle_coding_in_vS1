% Analyzing angle tuning from models
% How many of them "lose" angle tuning?
% How similar are the tuning?
% Which wkv factor affects tuning?
% 
% Data prepared from 'd190626_wkv_angle_tuning.m'
%         total 17 sets for testing angle tuning
%         (1) raw dF
%         (2) full model with touch angle
%         (3) full model with wkv
%         (4) from wkv, remove dKv (remove 5:7)
%         (5) from wkv, remove slide distance (remove 14:16)
%         (6) from wkv, remove touch count (remove 35:37)
%         (7) from wkv, remove arc length (remove 32:34)
%         (8) from wkv, remove dKh (remove 2:4)
%         (9) from wkv, remove dPhi (remove 11:13)
%         (10) from wkv, remove dKv + slide distance (remove [5:7, 14:16])
%         (11) from wkv, remove touch count + arc length (remove 32:37)
%         (12) from wkv, remove dKh + dPhi (remove [2:4, 11:13])
%         (13) from wkv, use only dKv (use [1, 5:7])
%         (14) from wkv, use only slide distance (use [1, 14:16])
%         (15) from wkv, use only dKv + slide distance (use [1, 5:7, 14:16])
%         (16) from wkv, use only touch count + arc length (use [1, 32:37])
%         (17) from wkv, use only dKh + dPhi (use [1:4, 11:13])

baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

%% first, check if tuning results of dF are the same from previous analysis

ratio = zeros(length(mice),1);
for mi = 1:12
si = 1;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = naive(mi);

ratio(mi) = length(find(spk.tuned))/length(find(currTuning.tuned));

end
1-mean(ratio)

%%

ratio = zeros(length(expert),1);
experti = find(cellfun(@(x) length(x) ==2, sessions));
for ei = 1:length(expert)
    mi = experti(ei);
si = 2;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = expert(ei);

ratio(ei) = length(find(spk.tuned))/length(find(currTuning.tuned));

end
1-mean(ratio)


%% About 10% (in naive) is lost when calculating angle tuning only from predecision touches
%% About 6% in expert

%% Are those from predecision touches all included in those from all touches?

included = zeros(length(mice),1);
for mi = 1:12
si = 1;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = naive(mi);

included(mi) = sum(ismember(find(spk.tuned), find(currTuning.tuned)))/length(find(spk.tuned));

end
mean(included)

included = zeros(length(expert),1);
experti = find(cellfun(@(x) length(x) ==2, sessions));
for ei = 1:length(expert)
    mi = experti(ei);
si = 2;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = expert(ei);

included(ei) = sum(ismember(find(spk.tuned), find(currTuning.tuned)))/length(find(spk.tuned));

end
mean(included)

%% 97.3% of naive and 99.5% of expert angle-tuned cells from predecision touch 
%% are included in all touch angle-tuned cells.

%% calculate match. (intersect)/(union)

match = zeros(length(mice),1);
for mi = 1:12
si = 1;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = naive(mi);

match(mi) = length(intersect(find(spk.tuned), find(currTuning.tuned)))/length(union(find(spk.tuned), find(currTuning.tuned)));

end
mean(match)

match = zeros(length(expert),1);
experti = find(cellfun(@(x) length(x) ==2, sessions));
for ei = 1:length(expert)
    mi = experti(ei);
si = 2;
mouse = mice(mi);
session = sessions{mi}(si);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_predecision', baseDir, mouse, mouse,session), 'spk')
currTuning = expert(ei);

match(ei) = length(intersect(find(spk.tuned), find(currTuning.tuned)))/length(union(find(spk.tuned), find(currTuning.tuned)));

end
mean(match)

%% Result: 85% & 93% match in naive and expert, respectively.



%% Gather all data in a single data file
clear
baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
experti = find(cellfun(@(x) length(x) == 2, sessions));

for i = 1 : 12
    mouse = mice(i);
    session = sessions{i}(1);
    naive(i) = load(sprintf('%s%03d\\angle_tuning_model_JK%03dS%02d', baseDir, mouse, mouse,session), 'tuneAngleAllCell', 'tunedAllCell');
end

for i = 1 : 6
    mouse = mice(experti(i));
    session = sessions{experti(i)}(2);
    expert(i) = load(sprintf('%s%03d\\angle_tuning_model_JK%03dS%02d', baseDir, mouse, mouse,session), 'tuneAngleAllCell', 'tunedAllCell');    
end
%% Comparing angle tuning between models

angleTuningNumNaive = zeros(12,17);
angleTuningNumExpert = zeros(6,17);

for i = 1 : 12
    for j = 1 : 17
        angleTuningNumNaive(i,j) = sum(naive(i).tunedAllCell(:,j)) / sum(naive(i).tunedAllCell(:,1));
    end
end

for i = 1 : 6
    for j = 1 : 17
        angleTuningNumExpert(i,j) = sum(expert(i).tunedAllCell(:,j)) / sum(expert(i).tunedAllCell(:,1));
    end
end

figure, hold on
shadedErrorBar(1:17, mean(angleTuningNumNaive), std(angleTuningNumNaive)/sqrt(12), 'lineprop', 'b')
shadedErrorBar(1:17, mean(angleTuningNumExpert), std(angleTuningNumExpert)/sqrt(6), 'lineprop', 'r')

%% Results:
% Kind of matching to what would have been expected, dKv (#3) and slide distance (#4) having the largest impact, even larger when combined together (#10)
% The impact is larger in naive, compared to expert.
% With only dKv, slide distance, or combined of these, the number of angle tuned cells are similar (even larger in naive) to dF, while touch count + arc length does not reproduce angle tuning (dKh + dPhi, a little bit)

%% Better calculation of angle tuning match.
% tuned angle similarity, using correlation values from angle tuned cells
clear
baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
experti = find(cellfun(@(x) length(x) == 2, sessions));

for i = 1 : 12
    mouse = mice(i);
    session = sessions{i}(1);
    naive(i) = load(sprintf('%s%03d\\angle_tuning_model_JK%03dS%02d', baseDir, mouse, mouse,session), 'tuneAngleAllCell', 'tunedAllCell', 'spkValAllCell');
end

for i = 1 : 6
    mouse = mice(experti(i));
    session = sessions{experti(i)}(2);
    expert(i) = load(sprintf('%s%03d\\angle_tuning_model_JK%03dS%02d', baseDir, mouse, mouse,session), 'tuneAngleAllCell', 'tunedAllCell', 'spkValAllCell');
end

%%
tuningSimilarityNaive = zeros(12,17);

for i = 1 : 12
    for j = 1 : 17
        similarityEachCell = zeros(size(naive(i).spkValAllCell,1),1);
        for k = 1 : length(similarityEachCell)
            similarityEachCell(k) = corr(cellfun(@(x) nanmean(x), naive(i).spkValAllCell{k,j}), cellfun(@(x) nanmean(x), naive(i).spkValAllCell{k,1}));
        end
        tuningSimilarityNaive(i,j) = nanmean(similarityEachCell);
    end
end

%
tuningSimilarityExpert = zeros(6,17);

for i = 1 : 6
    for j = 1 : 17
        similarityEachCell = zeros(size(expert(i).spkValAllCell,1),1);
        for k = 1 : length(similarityEachCell)
            similarityEachCell(k) = corr(cellfun(@(x) nanmean(x), expert(i).spkValAllCell{k,j}), cellfun(@(x) nanmean(x), expert(i).spkValAllCell{k,1}));
        end
        tuningSimilarityExpert(i,j) = nanmean(similarityEachCell);
    end
end

%
figure, hold on
shadedErrorBar(1:17, mean(tuningSimilarityNaive), std(tuningSimilarityNaive)/sqrt(12), 'lineprop', 'b')
shadedErrorBar(1:17, mean(tuningSimilarityExpert), std(tuningSimilarityExpert)/sqrt(6), 'lineprop', 'r')


%% What about only from angle-tuned cells?
tuningSimilarityNaive = zeros(12,17);
for i = 1 : 12
    for j = 1 : 17
        tunedCellInd = find(naive(i).tunedAllCell(:,1));
        similarityEachCell = zeros(length(tunedCellInd),1);
        for k = 1 : length(similarityEachCell)            
            similarityEachCell(k) = corr(cellfun(@(x) nanmean(x), naive(i).spkValAllCell{tunedCellInd(k),j}), cellfun(@(x) nanmean(x), naive(i).spkValAllCell{tunedCellInd(k),1}));
        end
        tuningSimilarityNaive(i,j) = nanmean(similarityEachCell);
    end
end

tuningSimilarityExpert = zeros(6,17);
for i = 1 : 6
    for j = 1 : 17
        tunedCellInd = find(expert(i).tunedAllCell(:,1));
        similarityEachCell = zeros(length(tunedCellInd),1);
        for k = 1 : length(similarityEachCell)
            similarityEachCell(k) = corr(cellfun(@(x) nanmean(x), expert(i).spkValAllCell{tunedCellInd(k),j}), cellfun(@(x) nanmean(x), expert(i).spkValAllCell{tunedCellInd(k),1}));
        end
        tuningSimilarityExpert(i,j) = nanmean(similarityEachCell);
    end
end

figure, hold on
shadedErrorBar(1:17, mean(tuningSimilarityNaive), std(tuningSimilarityNaive)/sqrt(12), 'lineprop', 'b')
shadedErrorBar(1:17, mean(tuningSimilarityExpert), std(tuningSimilarityExpert)/sqrt(6), 'lineprop', 'r')
xlabel('Model #'), ylabel('Correlation between mean angle responses')

%% Making more understandable figures, by grouping

%% group 1. naive
inds = [1,3,4,5,8,9];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityNaive,1))*i, tuningSimilarityNaive(:,inds(i)), 15, 'k')
end
for i = 1 : 12
    plot(1:length(inds), tuningSimilarityNaive(i,inds), 'k-')
end
xticks([1:length(inds)])
xticklabels({'Spikes', 'GLM', '- \Delta\kappa_V', '- Slide distance', '- \Delta\kappa_H', '- \Delta\phi', '- Touch count', '- Arc length'})
xtickangle(45)
ylabel({'Correlation between'; 'mean angle responses'})
set(gca,'fontsize',14)
%% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityNaive(:,3), tuningSimilarityNaive(:,inds(i)))
end

%% experts
inds = [1,3,4,5,8,9,6,7];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityExpert,1))*i, tuningSimilarityExpert(:,inds(i)), 15, 'k')
end
for i = 1 : 6
    plot(1:length(inds), tuningSimilarityExpert(i,inds), 'k-')
end

xticklabels({'Spikes', 'GLM', '- \Delta\kappa_V', '- Slide distance', '- \Delta\kappa_H', '- \Delta\phi', '- Touch count', '- Arc length'})
xtickangle(45)
ylabel('Correlation between mean angle responses')
% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityExpert(:,3), tuningSimilarityExpert(:,inds(i)))
end


%% group 2. naive
inds = [1,3,4,5,10,8,9,12];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityNaive,1))*i, tuningSimilarityNaive(:,inds(i)), 15, 'k')
end
for i = 1 : 12
    plot(1:length(inds), tuningSimilarityNaive(i,inds), 'k-')
end
xlim([1, length(inds)])
xticks([1:length(inds)])
xticklabels({'Spikes', 'GLM', '- \Delta\kappa_V', '- Slide distance', '-(\Delta\kappa_V & slide distance)','- \Delta\kappa_H', '- \Delta\phi', '-(\Delta\kappa_H & \Delta\phi)'})
% xticklabels({'Spikes', 'GLM', '- \Delta\kappa_V', '-(\Delta\kappa_V & slide distance)', '-(\Delta\kappa_H & \Delta\phi)', '-(Touch count & arc length)'})
xtickangle(45)
ylabel({'Correlation between'; 'mean angle responses'})
set(gca,'fontsize',14)
%% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityNaive(:,3), tuningSimilarityNaive(:,inds(i)))
end

%% group 2. experts
inds = [1,3,4,10,12,11];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityExpert,1))*i, tuningSimilarityExpert(:,inds(i)), 15, 'k')
end
for i = 1 : 6
    plot(1:length(inds), tuningSimilarityExpert(i,inds), 'k-')
end
xticks([1:length(inds)])
xticklabels({'Spikes', 'GLM', '- \Delta\kappa_V', '-(\Delta\kappa_V & slide distance)', '-(\Delta\kappa_H & \Delta\phi)', '-(Touch count & arc length)'})
xtickangle(45)
ylabel({'Correlation between'; 'mean angle responses'})
set(gca,'fontsize',14)
%% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityExpert(:,3), tuningSimilarityExpert(:,inds(i)))
end


%% group 3. naive
inds = [1,3,15,13,14,17];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityNaive,1))*i, tuningSimilarityNaive(:,inds(i)), 15, 'k')
end
for i = 1 : 12
    plot(1:length(inds), tuningSimilarityNaive(i,inds), 'k-')
end
xticks([1:length(inds)])
xticklabels({'Spikes', 'GLM', '+(\Delta\kappa_V & slide distance)', '+\Delta\kappa_V', '+Slide distance', '+(\Delta\kappa_H & \Delta\phi)'})
xtickangle(45)
ylabel({'Correlation between'; 'mean angle responses'})
set(gca,'fontsize',14)
%% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityNaive(:,3), tuningSimilarityNaive(:,inds(i)))
end

%% group 2. experts
inds = [1,3,13,14,15,17,16];
figure, hold on
for i = 1 : length(inds)
    scatter(ones(1,size(tuningSimilarityExpert,1))*i, tuningSimilarityExpert(:,inds(i)), 15, 'k')
end
for i = 1 : 6
    plot(1:length(inds), tuningSimilarityExpert(i,inds), 'k-')
end
xticks([1:length(inds)])
xticklabels({'Spikes', 'GLM', '+\Delta\kappa_V', '+Slide distance', '+(\Delta\kappa_V & slide distance)', '+(\Delta\kappa_H & \Delta\phi)', '+(Touch count & arc length)'})
xtickangle(45)
ylabel('Correlation between mean angle responses')
% calculating p-values
for i = 3 : length(inds)
    [~, p] = ttest(tuningSimilarityExpert(:,3), tuningSimilarityExpert(:,inds(i)))
end


%% Example graphs of angle-tuning.

%% First, find the right example.
mi = 1;
atind = find(naive(mi).tunedAllCell(:,1));
corrVal = zeros(length(atind),2);
for i = 1 : length(atind)
    spikeMean = cellfun(@nanmean, naive(mi).spkValAllCell{atind(i),1});
    wkvMean = cellfun(@nanmean, naive(mi).spkValAllCell{atind(i),3});
    wodkvMean = cellfun(@nanmean, naive(mi).spkValAllCell{atind(i),4});
    corrVal(i,1) = corr(spikeMean, wkvMean);
    corrVal(i,2) = corr(spikeMean, wodkvMean);
end
%%
[~,sind] = sort((corrVal(:,1) - corrVal(:,2)), 'descend');
indList = atind(sind);
%%
ci = 4;
angles = 45:15:135;
% dFCell = naive(mi).spkValAllCell{ci,1};
% maxVal = max(cell2mat(dFCell));
% minVal = min(cell2mat(dFCell));
% dFCellNorm = cellfun(@(x) (x-minVal)/(maxVal-minVal), dFCell, 'uniformoutput', false);
spikeMean = cellfun(@nanmean, naive(mi).spkValAllCell{ci,1});
spikeMean = (spikeMean - min(spikeMean)) / (max(spikeMean) - min(spikeMean));
spikeSem = cellfun(@(x) nanstd(x)/sqrt(length(x)), naive(mi).spkValAllCell{ci,1});


% wkvCell = naive(mi).spkValAllCell{ci,3};
% maxVal = max(cell2mat(dFCell));
% minVal = min(cell2mat(dFCell));
% wkvCellNorm = cellfun(@(x) (x-minVal)/(maxVal-minVal), wkvCell, 'uniformoutput', false);
wkvMean = cellfun(@nanmean, naive(mi).spkValAllCell{ci,3});
wkvMean = (wkvMean - min(wkvMean)) / (max(wkvMean) - min(wkvMean));
wkvSem = cellfun(@(x) nanstd(x)/sqrt(length(x)), naive(mi).spkValAllCell{ci,3});

% nondkvCell = naive(mi).spkValAllCell{ci,4};
% maxVal = max(cell2mat(dFCell));
% minVal = min(cell2mat(dFCell));
% nondkvCellNorm = cellfun(@(x) (x-minVal)/(maxVal-minVal), nondkvCell, 'uniformoutput', false);
wodkvMean = cellfun(@nanmean, naive(mi).spkValAllCell{ci,4});
wodkvMean = (wodkvMean - min(wodkvMean)) / (max(wodkvMean) - min(wodkvMean));
nondkvSem = cellfun(@(x) nanstd(x)/sqrt(length(x)), naive(mi).spkValAllCell{ci,4});

figure, hold on
% errorbar(angles, dFmean, dFsem, 'k.')
plot(angles, spikeMean, 'k.-', 'markersize', 20)

% errorbar(angles, wkvmean, wkvsem, 'b.')
plot(angles, wkvMean, 'b.-', 'markersize', 20)

% errorbar(angles, nondkvmean, nondkvsem, 'r.')
% plot(angles, nondkvmean, 'r.-', 'markersize', 20)

xticks(angles)
xlabel('Angle (\circ)')
ylabel('Normalized response pattern')
set(gca, 'fontsize', 15)

%
figure, hold on
% errorbar(angles, dFmean, dFsem, 'k.')
plot(angles, spikeMean, 'k.-', 'markersize', 20)

% errorbar(angles, wkvmean, wkvsem, 'b.')
plot(angles, wkvMean, 'b.-', 'markersize', 20)

% errorbar(angles, nondkvmean, nondkvsem, 'r.')
plot(angles, wodkvMean, 'r.-', 'markersize', 20)

xticks(angles)
xlabel('Angle (\circ)')
ylabel('Normalized response pattern')
set(gca, 'fontsize', 15)


%%
figure, 
scatter(spikeMean, wkvMean, 'b', 'filled')
xlabel('Spikes')
ylabel('GLM')
plt = gca;
plt.YAxis(1).Color = 'b';
set(gca, 'fontsize', 15)
%
figure, hold on
scatter(spikeMean, wkvMean, 'b', 'filled')

xlabel('Spikes')
ylabel('GLM')
yyaxis right
scatter(spikeMean, wodkvMean, 'r', 'filled')
ylabel('GLM - \DeltaK_v')

plt = gca;
plt.YAxis(1).Color = 'b';
plt.YAxis(2).Color = 'r';
set(gca, 'fontsize', 15)