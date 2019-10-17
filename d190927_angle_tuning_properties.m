%% 190927 Angle tuning properties

%% Purpose: Detailed description of angle tuning.
%     Modulation distribution
%     - How big is the modulation, if it’s angle-tuned?
%     - Is there a difference across tuned-angle?

%     Sharpness distribution
%     - How sharp is the response difference?
%     - How is it related to modulation, and tuned-angle?

%     % of response 
%     - How sparse is the response? 
%     - In both tuned-angle touches and all touches
 
%     Depth profile
%     - How are all these features distributed across depth?


% First, look at all naive. And then add expert and matching naive. 

%% basic settings and data loading
clear
baseDir = 'D:\TPM\JK\suite2p\';
loadFn = 'angle_tuning_summary_preAnswer_perTouch_NC';
load([baseDir, loadFn], 'naive', 'expert')


%% Modulation distribution

a = arrayfun(@(x) x.modulation, naive, 'un', 0)';
a = cell2mat(a);
figure, histogram(a)

b = sort(a, 'descend');
%% Result: max 2.7, top 16 out of 4567 > 2. Hist range 0:0.1:2;

histRange = [0:0.02:2,3];
distModulation = zeros(length(naive),length(histRange)-1);
for i = 1 : length(naive)
    distModulation(i,:) = histcounts(naive(i).modulation(find(naive(i).tuned)), histRange, 'normalization', 'probability');
end
figure,
errorbar(histRange(1:end-1), mean(distModulation), std(distModulation)/sqrt(length(naive)))
xlabel('Modulation'), ylabel('Proportion')

%% Sharpness distribution
a = arrayfun(@(x) x.sharpness, naive, 'un', 0)';
a = cell2mat(a);
figure, histogram(a)

b = sort(a, 'descend');

%%
histRange = [0:0.02:2,3];
distSharpness = zeros(length(naive),length(histRange)-1);
for i = 1 : length(naive)
    distSharpness(i,:) = histcounts(naive(i).sharpness(find(naive(i).tuned)), histRange, 'normalization', 'probability');
end
figure,
errorbar(histRange(1:end-1), mean(distSharpness), std(distSharpness)/sqrt(length(naive)))
xlabel('Sharpness'), ylabel('Proportion')

%%
figure, hold on
errorbar(histRange(1:end-1), mean(distModulation), std(distModulation)/sqrt(length(naive)))
errorbar(histRange(1:end-1), mean(distSharpness), std(distSharpness)/sqrt(length(naive)))
legend({'Modulation', 'Sharpness'}), ylabel('Proportion')


%% 
colors = parula(5);
single = [];
broad = [];
multi = [];
for i = 1 : length(naive)
    indSingle = find(naive(i).unimodalSingle);
    indBroad = find(naive(i).unimodalBroad);
    indMulti = setdiff(find(naive(i).multimodal), find(naive(i).unimodalBroad));
    single = [single; naive(i).modulation(indSingle), naive(i).sharpness(indSingle)];
    broad = [broad; naive(i).modulation(indBroad), naive(i).sharpness(indBroad)];
    multi = [multi; naive(i).modulation(indMulti), naive(i).sharpness(indMulti)];
end
single(find(single(:,2)<0),1) = -single(find(single(:,2)<0),1);
broad(find(broad(:,2)<0),1) = -broad(find(broad(:,2)<0),1);
multi(find(multi(:,2)<0),1) = -multi(find(multi(:,2)<0),1);
figure, hold on
scatter(single(:,1), single(:,2), [], colors(1,:));
scatter(broad(:,1), broad(:,2), [], colors(3,:));
scatter(multi(:,1), multi(:,2), [], colors(5,:));
legend({'Single', 'Broad', 'Complex'})
xlabel('Modultaion'), ylabel('Sharpness')

%% 
colors = parula(10);
props = zeros(length(naive), 3); % 1: single, 2: broad, 3: complex
single = [];
broad = [];
multi = [];
for i = 1 : length(naive)
    indSingle = find(naive(i).unimodalSingle);
    indBroad = find(naive(i).unimodalBroad);
    indMulti = setdiff(find(naive(i).multimodal), find(naive(i).unimodalBroad));
    single = [single; naive(i).modulation(indSingle), naive(i).sharpness(indSingle)];
    broad = [broad; naive(i).modulation(indBroad), naive(i).sharpness(indBroad)];
    multi = [multi; naive(i).modulation(indMulti), naive(i).sharpness(indMulti)];    
    props(i,1) = length(indSingle) / sum(naive(i).tuned);
    props(i,2) = length(indBroad) / sum(naive(i).tuned);
    props(i,3) = length(indMulti) / sum(naive(i).tuned);
end
single(find(single(:,2)<0),1) = -single(find(single(:,2)<0),1);
broad(find(broad(:,2)<0),1) = -broad(find(broad(:,2)<0),1);
multi(find(multi(:,2)<0),1) = -multi(find(multi(:,2)<0),1);
figure, hold on
scatter(single(:,1), single(:,2), 5, colors(1,:), 'filled');
scatter(broad(:,1), broad(:,2), 5, colors(6,:), 'filled');
scatter(multi(:,1), multi(:,2), 5, colors(9,:), 'filled');
legend({'Single', 'Broad', 'Complex'}, 'autoupdate', false, 'location', 'northwest')
xlabel('Modultaion'), ylabel('Sharpness')
axis equal
xlim([0 3]), ylim([0 3]),
plot([0 3], [0 3], '--', 'color', [0.7 0.7 0.7])


%%
bin = 0.1;
histRange = 0:bin:3.1;
figure, hold on
denominator = size(single,1) + size(broad,1) + size(multi,1);
factor = size(single,1) / denominator;
temp = histcounts(single(:,1), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(1,:))

factor = size(broad,1) / denominator;
temp = histcounts(broad(:,1), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(6,:))

factor = size(multi,1) / denominator;
temp = histcounts(multi(:,1), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(9,:))

ylim([0 0.1]), xlim([0 3])
xlabel('Modulation'), ylabel('Proportion')

%%
figure, hold on
denominator = size(single,1) + size(broad,1) + size(multi,1);
factor = size(single,1) / denominator;
temp = histcounts(single(:,2), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(1,:))

factor = size(broad,1) / denominator;
temp = histcounts(broad(:,2), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(6,:))

factor = size(multi,1) / denominator;
temp = histcounts(multi(:,2), histRange, 'norm', 'prob') * factor;
plot(histRange(1:end-1)+bin/2, temp, 'color', colors(9,:))

ylim([0 0.12]), xlim([0 3])
xlabel('Sharpness'), ylabel('Proportion')


%%
tempmat = flip(props,2);
allmean = mean(tempmat);
allsem = std(tempmat)/sqrt(length(naive));
labels = cell(3,1);
for i = 1 : 3
    labels{i} = sprintf('%.1f \\pm %.1f %%', allmean(i)*100, allsem(i)*100);
end
figure, pieH = pie(mean(tempmat), labels);
colormap(colors([9,6,1],:))
%% Depth profile of modulation 


%% No, before that, depth profile of tuning types

depthRange = [0:50:450,700];
depthType = zeros(length(naive), length(depthRange)-1, 3); % 1 single, 2 broad, 3 complex

for i = 1 : length(naive)
    for j = 1 : length(depthRange)-1
        depthType(i,j,1) = length(intersect( find(naive(i).unimodalSingle), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) )) / ...
            length(intersect( find(naive(i).tuned), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) ));
        depthType(i,j,2) = length(intersect( find(naive(i).unimodalBroad), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) )) / ...
            length(intersect( find(naive(i).tuned), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) ));
        depthType(i,j,3) = length(intersect( find(naive(i).multimodal), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) )) / ...
            length(intersect( find(naive(i).tuned), intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1))) ));
    end
end

figure, hold on

for i = 1 : 3
    temp = squeeze(depthType(:,:,i));    
    
    errorbar(depthRange(1:end-1)+25, nanmean(temp), nanstd(temp)/sqrt(length(naive)), 'color', colors((i-1)*2+1,:))
end
xlim([0 500])

legend({'Sharp', 'Broad', 'Complex'})

xlabel('Depth (mm)')
ylabel('Proportion (/tuned cell)')
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%% only from C2
colors = parula(5);
depthRange = [0:50:450,700];
depthType = zeros(length(naive), length(depthRange)-1, 3); % 1 single, 2 broad, 3 complex

for i = 1 : length(naive)
    for j = 1 : length(depthRange)-1
        indC2 = find(naive(i).isC2);
        indTuned = find(naive(i).tuned);
        indSharp = find(naive(i).unimodalSingle);
        indBroad = find(naive(i).unimodalBroad);
        indComplex = find(naive(i).multimodal);
        indDepth = intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1)));
        
        lengthTotal = length(intersect(intersect(indTuned, indC2), indDepth));
        
        depthType(i,j,1) = length( intersect( intersect(indSharp, indC2), indDepth ) ) / lengthTotal;
        depthType(i,j,2) = length( intersect( intersect(indBroad, indC2), indDepth ) ) / lengthTotal;
        depthType(i,j,3) = length( intersect( intersect(indComplex, indC2), indDepth ) ) / lengthTotal;
    end
end

figure, hold on

for i = 1 : 3
    temp = squeeze(depthType(:,:,i));    
    
    errorbar(depthRange(1:end-1)+25, nanmean(temp), nanstd(temp)/sqrt(length(naive)), 'color', colors((i-1)*2+1,:))
end
xlim([0 500])

legend({'Sharp', 'Broad', 'Complex'})

xlabel('Depth (mm)')
ylabel('Proportion (/tuned cell)')
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%% Result: When considering C2, there is no tendency in change of tuned type across depth.
%% All depth profiling analyses should be done within C2 only.
%% Do not need to consider properties between tuning type.
%% Since modulation is affected less from tuning type, just look at this, not sharpness.

depthRange = [0:50:450,700];
modulation = zeros(length(naive), length(depthRange)-1);
for i = 1 : length(naive)
    for j = 1 : length(depthRange)-1
        indC2Tuned = intersect(find(naive(i).isC2), find(naive(i).tuned));
        indDepth = intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1)));
        tempInd = intersect(indC2Tuned, indDepth);
        modulation(i,j) = mean(naive(i).modulation(tempInd));
    end        
end
figure,
errorbar(depthRange(1:end-1)+25, nanmean(modulation), nanstd(modulation)/sqrt(length(naive)))
xlabel('Depth'), ylabel('Modulation')


%% Just compare between L2/3 and L4. 
depthRange = [0, 350, 700];
modulation = zeros(length(naive), 2); % 1 L2/3, 2 L4
for i = 1 : length(naive)
    for j = 1 : 2
        indC2Tuned = intersect(find(naive(i).isC2), find(naive(i).tuned));
        indDepth = intersect(find(naive(i).depth >= depthRange(j)), find(naive(i).depth < depthRange(j+1)));
        tempInd = intersect(indC2Tuned, indDepth);
        modulation(i,j) = mean(naive(i).modulation(tempInd));
    end
end

figure, hold on
bar([1 2], mean(modulation), 'k')
errorbar([1 2], mean(modulation), std(modulation)/sqrt(length(naive)), 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'L2/3', 'L4'})
ylabel('Modulation')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%% Now, do this in noise-matched samples

clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
baseDir = 'D:\TPM\JK\suite2p\';
modulation = cell(2, length(mice));
noiseDist = cell(2, length(mice));
allNoiseDist = cell(2, length(mice));
numSamples = zeros(1,length(mice));
indNoiseMatched = cell(2, length(mice));
load(sprintf('%sangle_tuning_summary_preAnswer_perTouch_NC.mat',baseDir),'naive')

depthThresh = 350;
numNoiseBins = 20;
noiseThresh = [0.2, 0.6];

for mi = 1:length(naive)   
    fprintf('Processing %d/%d\n', mi, length(naive))
    
    indAllCellLayers = cell(2,1);
    indAllCellLayers{1} = intersect(intersect(find(naive(mi).depth < depthThresh), find(naive(mi).isC2)), find(naive(mi).tuned));
    indAllCellLayers{2} = intersect(intersect(find(naive(mi).depth >= depthThresh), find(naive(mi).isC2)), find(naive(mi).tuned));

    allNoiseDist{1,mi} = naive(mi).noise(indAllCellLayers{1});
    allNoiseDist{2,mi} = naive(mi).noise(indAllCellLayers{2});
    
    indCellNoise = intersect(find(naive(mi).noise > noiseThresh(1)), find(naive(mi).noise <= noiseThresh(2)));

    indCellLayers = cell(2,1);
    indCellLayers{1} = intersect(indCellNoise, indAllCellLayers{1});
    indCellLayers{2} = intersect(indCellNoise, indAllCellLayers{2});

    noiseRange = noiseThresh(1) : (noiseThresh(2) - noiseThresh(1)) / numNoiseBins : noiseThresh(2);
    numSample = zeros(length(noiseRange)-1,1);
    indNoiseLayers = cell(2,length(numSample));

    for ni = 1 : length(numSample)
        indNoiseLayers{1,ni} = indCellLayers{1}(intersect(find(naive(mi).noise(indCellLayers{1}) > noiseRange(ni)), find(naive(mi).noise(indCellLayers{1}) <= noiseRange(ni+1))));
        indNoiseLayers{2,ni} = indCellLayers{2}(intersect(find(naive(mi).noise(indCellLayers{2}) > noiseRange(ni)), find(naive(mi).noise(indCellLayers{2}) <= noiseRange(ni+1))));
        numSample(ni) = min(length(indNoiseLayers{1,ni}), length(indNoiseLayers{2,ni}));
    end
    numSampleInds = [0; cumsum(numSample)];
    numSamples(mi) = sum(numSample);
    for i = 1 : 2
        noiseDist{i,mi} = zeros(1,sum(numSample));
        modulation{i,mi} = zeros(1, sum(numSample));
        indGroup = [];
        for ni = 1 : length(numSample)
            tempIndGroup = indNoiseLayers{i,ni}(randperm(length(indNoiseLayers{i,ni}),numSample(ni)));
            noiseDist{i,mi}(numSampleInds(ni)+1:numSampleInds(ni+1)) = naive(mi).noise(tempIndGroup);
            modulation{i,mi}(numSampleInds(ni)+1:numSampleInds(ni+1)) = naive(mi).modulation(tempIndGroup);
            indGroup = [indGroup, tempIndGroup'];
        end
        indNoiseMatched{i,mi} = indGroup;
    end

end

%%
compInd = find(cellfun(@length, modulation(1,:))>=10);
modL23 = cellfun(@mean, modulation(1,compInd));
modL4 = cellfun(@mean, modulation(2,compInd));

figure, hold on
bar([1 2], [mean(modL23), mean(modL4)], 'k')
errorbar([1 2], [mean(modL23), mean(modL4)], [std(modL23)/sqrt(length(compInd)), std(modL4)/sqrt(length(compInd))], 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'L2/3', 'L4'})
ylabel('Modulation')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%%
histRange = 0:0.01:0.71;
figure, 
histogram(cell2mat(allNoiseDist(1,:)), histRange, 'normalization', 'cdf'), hold on
histogram(cell2mat(allNoiseDist(2,:)), histRange, 'normalization', 'cdf')

%%
histRange = 0.2:0.01:0.61;
figure, 
histogram(cell2mat(noiseDist(1,compInd)'), histRange, 'normalization', 'cdf'), hold on
histogram(cell2mat(noiseDist(2,compInd)'), histRange, 'normalization', 'cdf')

%% How about the signal? (95th percentile value?)

baseDir = 'D:\TPM\JK\suite2p\';
prctileThreshold = 95;

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

signal = cell(2,length(compInd));

for j = 1 : length(compInd)
    mouse = mice(compInd(j));
    session = sessions{compInd(j)}(1);
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session), 'u')
    for i = 1 : 2
        tempSignal = zeros(1,length(indNoiseMatched{i,compInd(j)}));
        for k = 1 : length(tempSignal)
            id = naive(compInd(j)).touchID(indNoiseMatched{i,compInd(j)}(k));
            tind = find(cellfun(@(x) ismember(id, x.neuindSession), u.trials));
            cind = find(u.trials{tind(1)}.neuindSession == id);
            tempData = cell2mat(cellfun(@(x) x.dF(cind,:), u.trials(tind)', 'un', 0));
            tempSignal(k) = prctile(tempData, prctileThreshold);
        end
        signal{i,j} = tempSignal;
    end
end

%%
histRange = [0.5:0.1:5.1, 10];
figure, 
plot(histRange(1:end-1), histcounts(cell2mat(signal(2,:)), histRange, 'normalization', 'cdf')), hold on
plot(histRange(1:end-1), histcounts(cell2mat(signal(1,:)), histRange, 'normalization', 'cdf'))
xlabel('Signal'), ylabel('Cumulative proportion')
legend({'L4', 'L2/3'})

%% Result: Signal is lower too. It's not sure if this lower signal is because of depth (scattering).

%% How about all L4 vs L2/3??
baseDir = 'D:\TPM\JK\suite2p\';
prctileThreshold = 95;

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

signal = cell(2,length(naive));

for j = 1 : length(naive)
    mouse = mice(j);
    session = sessions{j}(1);
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session), 'u')
    
        indL23 = intersect(intersect(find(naive(j).depth < 350), find(naive(j).isC2)), find(naive(j).tuned));
        tempSignal = zeros(1,length(indL23));
        for k = 1 : length(tempSignal)
            id = naive(j).touchID(indL23(k));
            tind = find(cellfun(@(x) ismember(id, x.neuindSession), u.trials));
            cind = find(u.trials{tind(1)}.neuindSession == id);
            tempData = cell2mat(cellfun(@(x) x.dF(cind,:), u.trials(tind)', 'un', 0));
            tempSignal(k) = prctile(tempData, prctileThreshold);
        end
        signal{1,j} = tempSignal;
    
        indL4 = intersect(intersect(find(naive(j).depth >= 350), find(naive(j).isC2)), find(naive(j).tuned));
        tempSignal = zeros(1,length(indL4));
        for k = 1 : length(tempSignal)
            id = naive(j).touchID(indL4(k));
            tind = find(cellfun(@(x) ismember(id, x.neuindSession), u.trials));
            cind = find(u.trials{tind(1)}.neuindSession == id);
            tempData = cell2mat(cellfun(@(x) x.dF(cind,:), u.trials(tind)', 'un', 0));
            tempSignal(k) = prctile(tempData, prctileThreshold);
        end
        signal{2,j} = tempSignal;
        
end

histRange = [0.5:0.1:5.1, 10];
figure, 
plot(histRange(1:end-1), histcounts(cell2mat(signal(2,:)), histRange, 'normalization', 'cdf')), hold on
plot(histRange(1:end-1), histcounts(cell2mat(signal(1,:)), histRange, 'normalization', 'cdf'))
xlabel('Signal'), ylabel('Cumulative proportion')
legend({'L4', 'L2/3'})

%%
figure, hold on
bar([1 2], mean(cellfun(@mean, signal')), 'k' )
errorbar([1 2], mean(cellfun(@mean, signal')), std(cellfun(@mean, signal'))/sqrt(length(naive)), 'k', 'linestyle', 'none');
xticks([1 2])
xticklabels({'L2/3', 'L4'})
ylabel('Signal')
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%% How does it look, in very low modulation tuned cells? (modulation < 0.1)

clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
baseDir = 'D:\TPM\JK\suite2p\';

load(sprintf('%sangle_tuning_summary_preAnswer_perTouch_NC.mat',baseDir),'naive')
%%
mi = 9;
mouse = mice(mi);
session = sessions{mi}(1);
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_predecision_NC',baseDir,mouse,mouse,session), 'ca')
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC',baseDir,mouse,mouse,session), 'spk', 'info')
load(sprintf('%s%03d\\UberJK%03dS%02d_NC',baseDir,mouse,mouse,session), 'u')


indTune = find(naive(mi).tuned);
mod = naive(mi).modulation(indTune);
cid = naive(mi).touchID(indTune);

[~, sorti] = sort(mod, 'ascend');

testi = 1;
mod(sorti(testi))
testCid = cid(sorti(testi));

example_angle_tuning_calcium(u, ca, spk, testCid)


%% Result: don't worry about it. anova p-value permutation test is enough.
% Lowest modulated cell among all naïve mice.
% JK052 S03 cell # 7113.
% Modulation 0.0277.
% Still looks tuned.
% 
% Tested all 12 lowest modulation cells from each mouse.

%% How does modulation look like between L2/3 C2 vs non-C2?

c2 = [1, 0];
modulation = zeros(length(naive), 2); % 1 L2/3 C2, 2 non-C2
for i = 1 : length(naive)
    for j = 1 : 2
        indL23 = find(naive(i).depth < 350);
        indC2Tuned = intersect(find(naive(i).isC2 == c2(j)), find(naive(i).tuned));
        
        tempInd = intersect(indC2Tuned, indL23);
        modulation(i,j) = mean(naive(i).modulation(tempInd));
    end
end

figure, hold on
bar([1 2], mean(modulation), 'k')
errorbar([1 2], mean(modulation), std(modulation)/sqrt(length(naive)), 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'L2/3 C2', 'non-C2'})
ylabel('Modulation')
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%% Per-trial sparseness
% It is not possible to get per-touch sparseness
% Calculated by the # of non-responsive trials (0 spikes / trial) at the tuned angle
binSize = 0.05;
histRange = [0:binSize:1];
sparseness = cell(length(naive),3); % 1 for tuned angle of tuned cells, 2 for not-tuned cells, 3 for min angle of tuned cells

for mi = 1 : length(naive)
    indTune = find(naive(mi).tuned);
    indTunedAngle = num2cell((naive(mi).tunedAngle(indTune) - 30) / 15);
    responses = cellfun(@(x,y) x{y}, naive(mi).val(indTune), indTunedAngle, 'un', 0);
    sparseness{mi,1} = cellfun(@(x) length(find(x == 0))/length(x), responses);
    
    indTune = find(naive(mi).tuned==0);    
    responses = cellfun(@(x) cell2mat(x), naive(mi).val(indTune), 'un', 0);
    sparseness{mi,2} = cellfun(@(x) length(find(x == 0))/length(x), responses);
    
    indLeast = cellfun(@(x) find(abs(cellfun(@mean, x)) == min(abs(cellfun(@mean, x))),1), naive(mi).val(indTune), 'un', 0);
    responses = cellfun(@(x,y) x{y}, naive(mi).val(indTune), indLeast, 'un', 0);
    sparseness{mi,3} = cellfun(@(x) length(find(x == 0))/length(x), responses);
    
end
figure, hold on
for i = 1 : 3
    temp = cellfun(@(x) histcounts(x, histRange, 'normalization', 'probability'), sparseness(:,i), 'un', 0);
    tempmat = cell2mat(temp);
    
    errorbar(histRange(1:end-1)+binSize/2, mean(tempmat), std(tempmat)/sqrt(length(naive)))
end

%%
legend({'Tuned (peak)', 'Not-tuned', 'Tuned (least)'})
xlabel('Sparsity'), ylabel('Proportion')
ylim([0 0.41])

%% comparison between L2/3 and L4 (only in C2)
meanSparsity = zeros(length(naive),2,2); % 1 for L2/3, 2 for L4 /// 1 for tuned cells (at tuned angle), 2 for not-tuned cells
for mi = 1 : length(naive)
    indTune = find(naive(mi).tuned);
    indDepth = find(naive(mi).depth < 350);
    indTest = intersect(indTune, indDepth);
    indTestAngle = num2cell((naive(mi).tunedAngle(indTest) - 30) / 15);
    responses = cellfun(@(x,y) x{y}, naive(mi).val(indTest), indTestAngle, 'un', 0);
    meanSparsity(mi,1,1) = mean(cellfun(@(x) length(find(x == 0))/length(x), responses));

    indDepth = find(naive(mi).depth >= 350);
    indTest = intersect(indTune, indDepth);
    indTestAngle = num2cell((naive(mi).tunedAngle(indTest) - 30) / 15);
    responses = cellfun(@(x,y) x{y}, naive(mi).val(indTest), indTestAngle, 'un', 0);
    meanSparsity(mi,2,1) = mean(cellfun(@(x) length(find(x == 0))/length(x), responses));

    
    indTune = find(naive(mi).tuned==0);
    indDepth = find(naive(mi).depth < 350);
    indTest = intersect(indTune, indDepth);
    responses = cellfun(@(x) cell2mat(x), naive(mi).val(indTest), 'un', 0);
    meanSparsity(mi,1,2) = mean(cellfun(@(x) length(find(x == 0))/length(x), responses));

    indTune = find(naive(mi).tuned==0);
    indDepth = find(naive(mi).depth >= 350);
    indTest = intersect(indTune, indDepth);
    responses = cellfun(@(x) cell2mat(x), naive(mi).val(indTest), 'un', 0);
    meanSparsity(mi,2,2) = mean(cellfun(@(x) length(find(x == 0))/length(x), responses));
end

%%
figure('units', 'normalized', 'outerposition', [0.2 0.2 0.2 0.5]), hold on
tempmat = squeeze(meanSparsity(:,1,:));
bar([0.7 2.3], mean(tempmat), 0.365, 'k')
tempmat = squeeze(meanSparsity(:,2,:));
bar([1.3 2.9], mean(tempmat), 0.365, 'w')

tempmat = squeeze(meanSparsity(:,1,:));
errorbar([0.7, 2.3], mean(tempmat), zeros(length(mean(tempmat)),1), std(tempmat)/sqrt(length(naive)), 'k', 'linestyle', 'none', 'linewidth', 1)
tempmat = squeeze(meanSparsity(:,2,:));
errorbar([1.3 2.9], mean(tempmat), zeros(length(mean(tempmat)),1), std(tempmat)/sqrt(length(naive)), 'k', 'linestyle', 'none', 'linewidth', 1)

xlim([0 3.5])
xticks([1 2.6])
xticklabels({'Tuned', 'Not-tuned'})
legend({'L2/3', 'L4'}, 'location', 'northwest')
ylabel('Sparsity')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

[~,p] = ttest(squeeze(meanSparsity(:,1,1)), squeeze(meanSparsity(:,2,1)))
[~,p] = ttest(squeeze(meanSparsity(:,1,2)), squeeze(meanSparsity(:,2,2)))
