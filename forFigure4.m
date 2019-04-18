%% for figure 4a
%% how much of tuned cells have whisking fit cells?

%% naive
clear
baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
tune = load('angle_tuning_summary');
glm = load('cellFunctionRidgeDE010');

whiskingProp = zeros(length(tune.naive),3); % 1: in total, 2: in touch cells, 3: in tuned cells

for i = 1 : length(tune.naive)
    tempW = glm.naive(i).whiskingID;
    whiskingProp(i,1) = length(tempW) / length(glm.naive(i).cellNums);
    whiskingProp(i,2) = length(intersect(tempW, glm.naive(i).touchID)) / length(glm.naive(i).touchID);
    whiskingProp(i,3) = length(intersect(tempW, tune.naive(i).touchID(find(tune.naive(i).tuned)))) / sum(tune.naive(i).tuned);
end
%
figure, 
bar(1:3, mean(whiskingProp), 'facecolor', 'w', 'linewidth', 2), hold on
errorbar(1:3, mean(whiskingProp), std(whiskingProp)/sqrt(length(tune.naive)), 'k.', 'linewidth', 2)
xticklabels({'In active cells', 'In touch cells', 'In angle-tuned cells'})
xtickangle(45)
ylabel('Whisking cell proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10, 'box', 'off')

%% expert, compared to matching naive

clear
baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
tune = load('angle_tuning_summary');
glm = load('cellFunctionRidgeDE010');
naiveInd = [1:4,7,9];
whiskingPropNaive = zeros(length(naiveInd),3); % 1: in total, 2: in touch cells, 3: in tuned cells
whiskingPropExpert = zeros(length(tune.expert),3); % 1: in total, 2: in touch cells, 3: in tuned cells
for i = 1 : length(naiveInd)
    ni = naiveInd(i);
    tempW = glm.naive(ni).whiskingID;
    whiskingPropNaive(i,1) = length(tempW) / length(glm.naive(ni).cellNums);
    whiskingPropNaive(i,2) = length(intersect(tempW, glm.naive(ni).touchID)) / length(glm.naive(ni).touchID);
    whiskingPropNaive(i,3) = length(intersect(tempW, tune.naive(ni).touchID(find(tune.naive(ni).tuned)))) / sum(tune.naive(ni).tuned);
    
    tempW = glm.expert(i).whiskingID;
    whiskingPropExpert(i,1) = length(tempW) / length(glm.expert(i).cellNums);
    whiskingPropExpert(i,2) = length(intersect(tempW, glm.expert(i).touchID)) / length(glm.expert(i).touchID);
    whiskingPropExpert(i,3) = length(intersect(tempW, tune.expert(i).touchID(find(tune.expert(i).tuned)))) / sum(tune.expert(i).tuned);
    
end
%
figure, 
bar([1:3]-0.2, mean(whiskingPropNaive), 0.4, 'facecolor', 'w', 'linewidth', 2), hold on
bar([1:3]+0.2, mean(whiskingPropExpert), 0.4, 'facecolor', 'k', 'linewidth', 2), hold on
errorbar([1:3]-0.2, mean(whiskingPropNaive), std(whiskingPropNaive)/sqrt(length(naiveInd)), 'k.', 'linewidth', 2)
errorbar([1:3]+0.2, mean(whiskingPropExpert), std(whiskingPropExpert)/sqrt(length(naiveInd)), 'k.', 'linewidth', 2)
xticks([1:3])
xticklabels({'In active cells', 'In touch cells', 'In angle-tuned cells'})
xtickangle(45)
ylabel('Whisking cell proportion')
legend({'Naive', 'Expert'}, 'box', 'off')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10, 'box', 'off')


%% DE comparison between touch and WTV fittings
clear
baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
decomp = load('glm_results_DEcomparison');
ercomp = load('glm_cell_function_error_ratio_withWTV');

%% first, see if DE comparison agrees with ER comparison

i = 1;
errExclusion = cellfun(@(x) x(1)/x(6), ercomp.naive(i).exclusionER); % err: error rate ratio
deDiffRatio = decomp.naive(i).touchWhiskerRatio;
figure, plot(errExclusion, deDiffRatio, 'k.')
hold on, 
inds = find(decomp.naive(i).touchDependence);
plot(errExclusion(inds), deDiffRatio(inds), 'r.')

%%

