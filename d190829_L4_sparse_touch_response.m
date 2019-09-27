% testing what might have led to sparse touch response in L4
% Is there something weird thing going on, or is it the real feature of L4?
% Possible reasons:
% Points to check: spike rate (it's maximum is controlled from MLspike) How
% aobut DE? 

% Results posted at 190524 Object angle coding interim summary.pptx


%% first, show the relationship with noise and neural activity
%% show if it's related to whether neurons are touch response or tuned

clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
baseDir = 'D:\TPM\JK\suite2p\';

allEventRates = cell(1, length(mice));
allNoise = cell(1, length(mice));

load(sprintf('%sangle_tuning_summary_predecision_NC.mat',baseDir),'naive')
for mi = 1:length(mice)
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(mi), mice(mi), sessions{mi}(1)))
    fprintf('Processing %d/%d\n', mi, length(mice))
    
    allEventRates = zeros(length(u.cellNums),1);
    for j = 1 : length(u.cellNums)
        tempCellId = u.cellNums(j);
        indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
        indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
        indPlane = mod(tempCellId-1,4)+1;
        poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
        allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));

        allEventRates(j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
    end
    noiseDist{i,mi} = zeros(sum(numSample),1);
    tempTouch = 0;
    tempTuned = 0;
    for ni = 1 : length(numSample)
        tempIndGroup = indNoiseLayers{i,ni}(randperm(length(indNoiseLayers{i,ni}),numSample(ni)));
        for j = 1  : length(tempIndGroup)
            tempCellId = u.cellNums(tempIndGroup(j));
            indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
            indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
            indPlane = mod(tempCellId-1,4)+1;
            poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
            allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));

            eventRatesLayers{i,mi}(numSampleInds(ni)+j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);

        end
        noiseDist{i,mi}(numSampleInds(ni)+1:numSampleInds(ni+1)) = u.noise(tempIndGroup);
        tempTouch = tempTouch + length(find(ismember(naive(mi).touchID, u.cellNums(tempIndGroup))));
        tempTuned = tempTuned + length(find(ismember(naive(mi).touchID(find(naive(mi).tuned)), u.cellNums(tempIndGroup))));
    end
    numTouch(i,mi) = tempTouch;
    numTuned(i,mi) = tempTuned;

end


%%
clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
mi = 1;

% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'D:\TPM\JK\suite2p\';
% ac = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_all_cell.mat', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
% glm = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning.mat', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
load(sprintf('%scellFunctionLasso_NC.mat', baseDir), 'naive')
load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(mi), mice(mi), sessions{mi}(1)))

%% Distribution of spike rates, L2, L3, and L4
range = [1, 136, 350, 700]; % for L2, L3, L4
indAllCellLayers = cell(length(range)-1,1);
for i = 1 : length(indAllCellLayers)
    indAllCellLayers{i} = intersect(find(u.cellDepths > range(i)), find(u.cellDepths <= range(i+1)));
end

spikeRatesLayers = cell(length(range)-1,1);
for i = 1 : length(spikeRatesLayers)
    spikeRatesLayers{i} = zeros(length(indAllCellLayers{i}),1);
    for j = 1 : length(indAllCellLayers{i})
        tempCellId = u.cellNums(indAllCellLayers{i}(j));
        indTrial = find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials));
        indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);        
        allSpikes = cell2mat(cellfun(@(x) x.spk(indCellSession,:), u.trials(indTrial)', 'uniformoutput', false));
        spikeRatesLayers{i}(j) = sum(allSpikes) / length(allSpikes);
    end
end

%%
range = 0:0.01:0.6;
figure, hold on,
for i = 1 : length(spikeRatesLayers)
    plot(range(1:end-1), histcounts(spikeRatesLayers{i}, range, 'normalization', 'probability'))
end
legend({'L2','L3','L4'}), xlabel('Spike rate (modulated)'), ylabel('Proportion'), title(sprintf('JK%03d', mice(mi)))

% L4 spike rate is slightly lower than L3, especially at > 0.1. 
% What about during pole up periods?

%% Distribution of spike rates, L2, L3, and L4 only during pole up periods
range = [1, 136, 350, 700]; % for L2, L3, L4
indAllCellLayers = cell(length(range)-1,1);
for i = 1 : length(indAllCellLayers)
    indAllCellLayers{i} = intersect(find(u.cellDepths > range(i)), find(u.cellDepths <= range(i+1)));
end

spikeRatesLayers = cell(length(range)-1,1);
for i = 1 : length(spikeRatesLayers)
    spikeRatesLayers{i} = zeros(length(indAllCellLayers{i}),1);
    for j = 1 : length(indAllCellLayers{i})
        tempCellId = u.cellNums(indAllCellLayers{i}(j));
        indTrial = find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials));
        indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
        indPlane = mod(tempCellId-1,4)+1;
        poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
        allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));
        spikeRatesLayers{i}(j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
    end
end

%%
range = 0:0.01:0.6;
figure, hold on,
for i = 1 : length(spikeRatesLayers)
    plot(range(1:end-1), histcounts(spikeRatesLayers{i}, range, 'normalization', 'probability'))
end
legend({'L2','L3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title(sprintf('JK%03d - During pole up', mice(mi)))

% Similar to all frames

%% from all mice (naive)

clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
baseDir = 'Y:\Whiskernas\JK\suite2p\';
% depthRange = [1, 136, 350, 700]; % for L2, L3, L4
depthRange = [1, 350, 700]; % for L2/3 and L4
histRange = 0:0.01:0.3;
spikeRatesLayers = cell(length(depthRange)-1, length(mice));

for mi = 1 : length(mice)
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(mi), mice(mi), sessions{mi}(1)))
    fprintf('Processing %d/%d\n', mi, length(mice))
    indAllCellLayers = cell(length(depthRange)-1,1);
    for i = 1 : length(indAllCellLayers)
        indAllCellLayers{i} = intersect(find(u.cellDepths > depthRange(i)), find(u.cellDepths <= depthRange(i+1)));
    end
    for i = 1 : size(spikeRatesLayers,1)
        spikeRatesLayers{i,mi} = zeros(length(indAllCellLayers{i}),1);
        for j = 1 : length(indAllCellLayers{i})
            tempCellId = u.cellNums(indAllCellLayers{i}(j));
            indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
            indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
            indPlane = mod(tempCellId-1,4)+1;
            poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
            allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));
            spikeRatesLayers{i,mi}(j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
        end
    end
end

histSpikeRatesLayers = zeros(length(depthRange)-1, length(mice), length(histRange)-1);

for di = 1 : length(depthRange) - 1
    for mi = 1 : length(mice)
        histSpikeRatesLayers(di,mi,:) = histcounts(spikeRatesLayers{di,mi}, histRange, 'normalization', 'probability');
    end
end

colors = get(gca,'colororder');
figure, hold all,
for i = 1 : length(depthRange) - 1
    tempMat = squeeze(histSpikeRatesLayers(i,:,:));
    errorbar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(mice)), 'color',colors(i,:))
end
for i = 1 : length(depthRange) - 1
    plot(histRange(1:end-1), mean(squeeze(histSpikeRatesLayers(i,:,:))), 'color',colors(i,:))
end
% legend({'L2','L3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title('Averaged from all mice (naive)')
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title('Averaged from all mice (naive)')
ylim([0 0.25])

%%
cumhistRange = 0:0.01:0.5;
cumhistSpikeRatesLayers = zeros(length(depthRange)-1, length(mice), length(cumhistRange)-1);

for di = 1 : length(depthRange) - 1
    for mi = 1 : length(mice)
        cumhistSpikeRatesLayers(di,mi,:) = histcounts(spikeRatesLayers{di,mi}, cumhistRange, 'normalization', 'cdf');
    end
end

colors = get(gca,'colororder');
figure, hold all,
for i = 1 : length(depthRange) - 1
    tempMat = squeeze(cumhistSpikeRatesLayers(i,:,:));
    errorbar(cumhistRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(mice)), 'color',colors(i,:))
end
for i = 1 : length(depthRange) - 1
    plot(cumhistRange(1:end-1), mean(squeeze(cumhistSpikeRatesLayers(i,:,:))), 'color',colors(i,:))
end
% legend({'L2','L3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title('Averaged from all mice (naive)')
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Cumulative proportion'), title('Averaged from all mice (naive)')


%% what if I use the similar noise level between layers?
% Just compare between L2/3 and L4
% match the noise distribution between the two layers, by dividing into 10
% bins and then randomly sampling from these bins (matched with the lower
% number of samples)

depthThresh = 350;
prctileThresh = 20; % percentile in either direction.
histRange = 0:0.01:0.3;
numNoiseBins = 10;
mi = 1;
load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(mi), mice(mi), sessions{mi}(1)))

indAllCellLayers = cell(2,1);
indAllCellLayers{1} = find(u.cellDepths < depthThresh);
indAllCellLayers{2} = find(u.cellDepths >= depthThresh);

noiseThresh = zeros(1,2);
noiseThresh(1) = prctile(u.noise(indAllCellLayers{2}), prctileThresh);
noiseThresh(2) = prctile(u.noise(indAllCellLayers{1}), 100-prctileThresh);
indCellNoise = intersect(find(u.noise > noiseThresh(1)), find(u.noise <= noiseThresh(2)));

indCellLayers = cell(2,1);
indCellLayers{1} = intersect(indCellNoise, indAllCellLayers{1});
indCellLayers{2} = intersect(indCellNoise, indAllCellLayers{2});

noiseRange = noiseThresh(1) : (noiseThresh(2) - noiseThresh(1)) / numNoiseBins : noiseThresh(2);
noiseRange = sort(noiseRange);
numSample = zeros(length(noiseRange)-1,1);
indNoiseLayers = cell(2,length(numSample));

for ni = 1 : length(numSample)
    indNoiseLayers{1,ni} = indCellLayers{1}(intersect(find(u.noise(indCellLayers{1}) > noiseRange(ni)), find(u.noise(indCellLayers{1}) <= noiseRange(ni+1))));
    indNoiseLayers{2,ni} = indCellLayers{2}(intersect(find(u.noise(indCellLayers{2}) > noiseRange(ni)), find(u.noise(indCellLayers{2}) <= noiseRange(ni+1))));
    numSample(ni) = min(length(indNoiseLayers{1,ni}), length(indNoiseLayers{2,ni}));
end
numSampleInds = [0; cumsum(numSample)];

eventRatesLayers = cell(2,1);
indGroups = cell(2,length(numSample));
for i = 1 : 2
    eventRatesLayers{i} = zeros(sum(numSample),1);
    for ni = 1 : length(numSample)
        tempIndGroup = indNoiseLayers{i,ni}(randperm(length(indNoiseLayers{i,ni}),numSample(ni)));
        indGroups{i,ni} = tempIndGroup;
        for j = 1  : length(tempIndGroup)
            tempCellId = u.cellNums(tempIndGroup(j));
            indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
            indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
            indPlane = mod(tempCellId-1,4)+1;
            poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
            allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));
            eventRatesLayers{i}(numSampleInds(ni)+j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
        end
    end
end

figure, hold on,
for i = 1 : 2
    plot(histRange(1:end-1), histcounts(eventRatesLayers{i}, histRange, 'normalization', 'probability'))
end
% noise distribution
figure, hold on,
for i = 1 : 2
    histogram(u.noise(cell2mat(indGroups(i,:))), noiseRange)
end
%% from all mice (naive)
clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
baseDir = 'D:\TPM\JK\suite2p\';
eventRatesLayers = cell(2, length(mice));
allEventRatesLayers = cell(2, length(mice));
noiseDist = cell(2, length(mice));
allNoiseDist = cell(2, length(mice));
numSamples = zeros(1,length(mice));
numTouch = zeros(2,length(mice));
numTuned = zeros(2,length(mice));
depthThresh = 350;
prctileThresh = 20; % percentile in either direction.
histRange = 0:0.01:0.3;
numNoiseBins = 20;

load(sprintf('%sangle_tuning_summary_predecision_NC.mat',baseDir),'naive')
for mi = 1:length(mice)
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(mi), mice(mi), sessions{mi}(1)))
    fprintf('Processing %d/%d\n', mi, length(mice))
    
    indAllCellLayers = cell(2,1);
    indAllCellLayers{1} = find(u.cellDepths < depthThresh);
    indAllCellLayers{2} = find(u.cellDepths >= depthThresh);

    allNoiseDist{1,mi} = u.noise(indAllCellLayers{1});
    allNoiseDist{2,mi} = u.noise(indAllCellLayers{2});
    
    noiseThresh = zeros(1,2);
    noiseThresh(1) = prctile(u.noise(indAllCellLayers{2}), prctileThresh);
    noiseThresh(2) = prctile(u.noise(indAllCellLayers{1}), 100-prctileThresh);
    indCellNoise = intersect(find(u.noise > noiseThresh(1)), find(u.noise <= noiseThresh(2)));

    indCellLayers = cell(2,1);
    indCellLayers{1} = intersect(indCellNoise, indAllCellLayers{1});
    indCellLayers{2} = intersect(indCellNoise, indAllCellLayers{2});

    noiseRange = noiseThresh(1) : (noiseThresh(2) - noiseThresh(1)) / numNoiseBins : noiseThresh(2);
    numSample = zeros(length(noiseRange)-1,1);
    indNoiseLayers = cell(2,length(numSample));

    for ni = 1 : length(numSample)
        indNoiseLayers{1,ni} = indCellLayers{1}(intersect(find(u.noise(indCellLayers{1}) > noiseRange(ni)), find(u.noise(indCellLayers{1}) <= noiseRange(ni+1))));
        indNoiseLayers{2,ni} = indCellLayers{2}(intersect(find(u.noise(indCellLayers{2}) > noiseRange(ni)), find(u.noise(indCellLayers{2}) <= noiseRange(ni+1))));
        numSample(ni) = min(length(indNoiseLayers{1,ni}), length(indNoiseLayers{2,ni}));
    end
    numSampleInds = [0; cumsum(numSample)];
    numSamples(mi) = sum(numSample);
    for i = 1 : 2
        eventRatesLayers{i,mi} = zeros(sum(numSample),1);
        allEventRatesLayers{i,mi} = zeros(length(u.cellNums),1);
        for j = 1 : length(u.cellNums)
            tempCellId = u.cellNums(j);
            indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
            indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
            indPlane = mod(tempCellId-1,4)+1;
            poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
            allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));

            allEventRatesLayers{i,mi}(j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
        end
        noiseDist{i,mi} = zeros(sum(numSample),1);
        tempTouch = 0;
        tempTuned = 0;
        for ni = 1 : length(numSample)
            tempIndGroup = indNoiseLayers{i,ni}(randperm(length(indNoiseLayers{i,ni}),numSample(ni)));
            for j = 1  : length(tempIndGroup)
                tempCellId = u.cellNums(tempIndGroup(j));
                indTrial = setdiff(find(cellfun(@(x) ~isempty(find(x.neuindSession == tempCellId)), u.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), u.trials)));
                indCellSession = find(u.trials{indTrial(1)}.neuindSession == tempCellId);
                indPlane = mod(tempCellId-1,4)+1;
                poleUpFrames = cellfun(@(x) find(x.tpmTime{indPlane} >= x.poleUpTime(1) & x.tpmTime{indPlane} <= x.poleUpTime(end)), u.trials(indTrial), 'uniformoutput', false);
                allSpikesPoleUp = cell2mat(cellfun(@(x,y) x.spk(indCellSession,y), u.trials(indTrial)', poleUpFrames', 'uniformoutput', false));
                
                eventRatesLayers{i,mi}(numSampleInds(ni)+j) = sum(allSpikesPoleUp) / length(allSpikesPoleUp);
                
            end
            noiseDist{i,mi}(numSampleInds(ni)+1:numSampleInds(ni+1)) = u.noise(tempIndGroup);
            tempTouch = tempTouch + length(find(ismember(naive(mi).touchID, u.cellNums(tempIndGroup))));
            tempTuned = tempTuned + length(find(ismember(naive(mi).touchID(find(naive(mi).tuned)), u.cellNums(tempIndGroup))));
        end
        numTouch(i,mi) = tempTouch;
        numTuned(i,mi) = tempTuned;
    end
end
%% Noise distribution comparison
figure, hold on
for mi = 1 : length(mice)
    temp = histcounts(noiseDist{1,mi}, noiseRange);
    plot(noiseRange(1:end-1), temp, 'k-')
    temp = histcounts(noiseDist{2,mi}, noiseRange);
    plot(noiseRange(1:end-1), temp, 'b-')
end

%%
sampleInd = find(numSamples>=10);
noiseRange = [0.25:0.01:0.45,1];
noiseHist = zeros(length(sampleInd),length(noiseRange)-1,2);
for mi = 1 : length(sampleInd)
    noiseHist(mi,:,1) = histcounts(noiseDist{1,sampleInd(mi)}, noiseRange, 'normalization', 'probability');
    noiseHist(mi,:,2) = histcounts(noiseDist{2,sampleInd(mi)}, noiseRange, 'normalization', 'probability');
end
figure, hold on
errorbar(noiseRange(1:end-1), mean(squeeze(noiseHist(:,:,1))), std(squeeze(noiseHist(:,:,1)))/sqrt(length(sampleInd)))
errorbar(noiseRange(1:end-1), mean(squeeze(noiseHist(:,:,2))), std(squeeze(noiseHist(:,:,2)))/sqrt(length(sampleInd)))
legend({'L2/3','L4'})
xlabel('Noise')
ylabel('Proportion')
title('Noise-matched')
set(gca,'fontsize',12,'fontname','Arial')
%% How about from all noise level?

noiseRange = [0.1:0.01:0.7,1];
allNoiseHist = zeros(length(mice),length(noiseRange)-1,2);
for mi = 1 : length(mice)
    allNoiseHist(mi,:,1) = histcounts(allNoiseDist{1,mi}, noiseRange, 'normalization', 'probability');
    allNoiseHist(mi,:,2) = histcounts(allNoiseDist{2,mi}, noiseRange, 'normalization', 'probability');
end
figure, hold on
errorbar(noiseRange(1:end-1), mean(squeeze(allNoiseHist(:,:,1))), std(squeeze(allNoiseHist(:,:,1)))/sqrt(length(mice)))
errorbar(noiseRange(1:end-1), mean(squeeze(allNoiseHist(:,:,2))), std(squeeze(allNoiseHist(:,:,2)))/sqrt(length(mice)))
legend({'L2/3','L4'})
xlabel('Noise')
ylabel('Proportion')
title('All neurons')
set(gca,'fontsize',12,'fontname','Arial')

%% event rate distribution - noise matched
histEventRatesLayers = zeros(2, length(mice), length(histRange)-1);

for di = 1 : 2
    for mi = 1 : length(mice)
        histEventRatesLayers(di,mi,:) = histcounts(eventRatesLayers{di,mi}, histRange, 'normalization', 'probability');
    end
end
%
colors = get(gca,'colororder');
figure, hold all,
for i = 1 : 2
    tempMat = squeeze(histEventRatesLayers(i,:,:));
    errorbar(histRange(1:end-1), nanmean(tempMat), nanstd(tempMat)/sqrt(sum(isfinite(sum(tempMat,2)))), 'color',colors(i,:))
end
for i = 1 : 2
    plot(histRange(1:end-1), nanmean(squeeze(histEventRatesLayers(i,:,:))), 'color',colors(i,:))
end
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title('Averaged from all mice (naive), noise level-matched')
title('Noise-matched')
set(gca,'fontsize',12,'fontname','Arial')
%% cumulative - noise matched
cumhistRange = 0 :0.01:0.5;
histEventRatesLayers = zeros(2, length(mice), length(cumhistRange)-1);

for di = 1 : 2
    for mi = 1 : length(mice)
        histEventRatesLayers(di,mi,:) = histcounts(eventRatesLayers{di,mi}, cumhistRange, 'normalization', 'cdf');
    end
end

colors = get(gca,'colororder');
figure, hold all,
for i = 1 : 2
    tempMat = squeeze(histEventRatesLayers(i,:,:));
    errorbar(cumhistRange(1:end-1), nanmean(tempMat), nanstd(tempMat)/sqrt(sum(isfinite(sum(tempMat,2)))), 'color',colors(i,:))
end
for i = 1 : 2
    plot(cumhistRange(1:end-1), mean(squeeze(histEventRatesLayers(i,:,:))), 'color',colors(i,:))
end
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Cumulative proportion'), title('Averaged from all mice (naive), noise level-matched')
title('Noise-matched')
set(gca,'fontsize',12,'fontname','Arial')

%% event rate distribution - all cells
histAllEventRatesLayers = zeros(2, length(mice), length(histRange)-1);

for di = 1 : 2
    for mi = 1 : length(mice)
        histAllEventRatesLayers(di,mi,:) = histcounts(eventRatesLayers{di,mi}, histRange, 'normalization', 'probability');
    end
end
%
colors = get(gca,'colororder');
figure, hold all,
for i = 1 : 2
    tempMat = squeeze(histAllEventRatesLayers(i,:,:));
    errorbar(histRange(1:end-1), nanmean(tempMat), nanstd(tempMat)/sqrt(sum(isfinite(sum(tempMat,2)))), 'color',colors(i,:))
end
for i = 1 : 2
    plot(histRange(1:end-1), nanmean(squeeze(histAllEventRatesLayers(i,:,:))), 'color',colors(i,:))
end
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Proportion'), title('Averaged from all mice (naive), noise level-matched')
title('All neurons')
set(gca,'fontsize',12,'fontname','Arial')

%% cumulative - all cells
% cumhistRange = 0 :0.01:0.5;
cumhistRange = histRange;
histAllEventRatesLayers = zeros(2, length(mice), length(cumhistRange)-1);

for di = 1 : 2
    for mi = 1 : length(mice)
        histAllEventRatesLayers(di,mi,:) = histcounts(eventRatesLayers{di,mi}, cumhistRange, 'normalization', 'cdf');
    end
end

colors = get(gca,'colororder');
figure, hold all,
for i = 1 : 2
    tempMat = squeeze(histAllEventRatesLayers(i,:,:));
    errorbar(cumhistRange(1:end-1), nanmean(tempMat), nanstd(tempMat)/sqrt(sum(isfinite(sum(tempMat,2)))), 'color',colors(i,:))
end
for i = 1 : 2
    plot(cumhistRange(1:end-1), mean(squeeze(histAllEventRatesLayers(i,:,:))), 'color',colors(i,:))
end
legend({'L2/3','L4'}), xlabel('Spike rate (modulated; /frame)'), ylabel('Cumulative proportion'), title('Averaged from all mice (naive), noise level-matched')
title('All neurons')
set(gca,'fontsize',12,'fontname','Arial')

%% proportion of touch responsive and tuned cells within these populations
propTouch = zeros(length(sampleInd),2);
propTuned = zeros(length(sampleInd),2);
for i = 1 : 2
    for mi = 1 : length(sampleInd)
        propTouch(mi,i) = numTouch(i,sampleInd(mi)) / numSamples(sampleInd(mi));
        propTuned(mi,i) = numTuned(i,sampleInd(mi)) / numSamples(sampleInd(mi));
    end
end
figure,
hold on
bar(1:2, mean(propTouch), 'b')
bar(1:2, mean(propTuned), 'facecolor', [0.7 0.7 0.7])
errorbar(1:2, mean(propTouch), std(propTouch)/sqrt(length(sampleInd)), 'b', 'linestyle', 'none')
errorbar(1:2, mean(propTuned), std(propTuned)/sqrt(length(sampleInd)), 'color', [0.7 0.7 0.7], 'linestyle', 'none')
bar(1:2, mean(propTuned), 'facecolor', [0.7 0.7 0.7])
legend({'Touch', 'Tuned'}, 'box', 'off')

xticks(1:2)
xticklabels({'L2/3', 'L4'})
ylabel('Proportion')
title('Noise-matched')
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%%
propTouch = zeros(length(sampleInd),2);
propTuned = zeros(length(sampleInd),2);
for mi = 1 : length(sampleInd)
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mice(sampleInd(mi)), mice(sampleInd(mi)), sessions{sampleInd(mi)}(1)))
    layerInds{1} = find(u.cellDepths < 350);
    layerInds{2} = find(u.cellDepths >= 350);
    touchID = naive(sampleInd(mi)).touchID;
    tunedID = naive(sampleInd(mi)).touchID(find(naive(sampleInd(mi).tuned)));
    for i = 1 : 2
        
        propTouch(mi,i) = numTouch(i,sampleInd(mi)) / numSamples(sampleInd(mi));
        propTuned(mi,i) = numTuned(i,sampleInd(mi)) / numSamples(sampleInd(mi));
    end
end
figure,
hold on
bar(1:2, mean(propTouch), 'b')
bar(1:2, mean(propTuned), 'facecolor', [0.7 0.7 0.7])
errorbar(1:2, mean(propTouch), std(propTouch)/sqrt(length(sampleInd)), 'b', 'linestyle', 'none')
errorbar(1:2, mean(propTuned), std(propTuned)/sqrt(length(sampleInd)), 'color', [0.7 0.7 0.7], 'linestyle', 'none')
bar(1:2, mean(propTuned), 'facecolor', [0.7 0.7 0.7])
legend({'Touch', 'Tuned'}, 'box', 'off')

xticks(1:2)
xticklabels({'L2/3', 'L4'})
ylabel('Proportion')
title('Noise-matched')
set(gca, 'fontsize', 12, 'fontname', 'Arial')



%% Conclusion
% Smaller number of touch responsive cells in L4 compared to L2/3 is NOT because of higher noise level