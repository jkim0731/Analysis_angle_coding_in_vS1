clear all; close all; clc;

load('Y:/Whiskernas/JK/suite2p/angle_tuning_summary');

%% params
dataSet = expert;
animalID = 1;
% angles = linspace(45, 135, 7);
angles = [45 45 45 90 135 135 135];
layers = [0 5000];
angleDiscrim = [45; 135];

%% reshape data
tableData = [];
properties = [];
for roiID = 1:length(dataSet(animalID).touchID)
    roiName = dataSet(animalID).touchID(roiID);
    if (roiName > layers(1) && roiName < layers(2))
        responses = dataSet(animalID).val{roiID};

        angleVals = [];
        responseVals = [];
        for a = 1:length(responses)
            angleVals = [angleVals; repmat(angles(a), length(responses{a}), 1)];
            responseVals = [responseVals; responses{a}];
        end
        properties = [properties; [roiName, dataSet(animalID).tuned(roiID),...
                                   dataSet(animalID).tunedAngle(roiID),...
                                   dataSet(animalID).tuneDirection(roiID),...
                                   dataSet(animalID).sharpness(roiID),...
                                   dataSet(animalID).modulation(roiID)]];
        tableData = [tableData, responseVals];
    end
end

%% pca
[coeff, score, latent, tsquared, explained, mu] = pca(tableData);
figure;
plot(explained);

figure;
plot(coeff(:, 1));

%% train classifier
nCells = size(tableData, 2);

accuracy = nan(reps, nCells);
accuracyShuffled = nan(reps, nCells);
allCellChoice = cell(reps, nCells);
allClassifiers = cell(reps, nCells);

for i = 1:nCells
    disp(['Running reps for ' num2str(i) ' cells.']);
    for j = 1:reps
        % choose a permutation of cells to use
        cellChoice = randperm(nCells, i);
        allCellChoice{j, i} = cellChoice;

        % data set
        classificationSet = [tableData(:, cellChoice), angleVals];
        dIdx = find(angleVals == angleDiscrim(1) | angleVals == angleDiscrim(2));
        classificationSet = (classificationSet(dIdx, :));
        % reorder rows (labels intact)
        classificationSet = classificationSet(randperm(size(classificationSet, 1)), :);

        % shuffled data set (reorder labels only)
        classificationSetShuffled = classificationSet;
        classificationSetShuffled(:, end) = classificationSetShuffled(randperm(size(classificationSetShuffled,1)), end);

        [trainedClassifier, validationAccuracy] = trainClassifier(classificationSet, angleDiscrim);
        [trainedClassifierShuffled, validationAccuracyShuffled] = trainClassifier(classificationSetShuffled, angleDiscrim);

        accuracy(j, i) = validationAccuracy;
        accuracyShuffled(j, i) = validationAccuracyShuffled;
        allClassifiers{j, i} = trainedClassifier;
    end
end

%% plot
figure; hold on;
shadedErrorBar(1:nCells, mean(accuracy), std(accuracy), 'b')
shadedErrorBar(1:nCells, mean(accuracyShuffled), std(accuracyShuffled), 'r')

%% t-sne
figure; hold on;
Y = tsne(tableData,'NumDimensions',3);
scatter3(Y(angleVals==45,1), Y(angleVals==45,2), Y(angleVals==45,3), 'r')
scatter3(Y(angleVals==135,1), Y(angleVals==135,2), Y(angleVals==135,3), 'b')
scatter3(Y(angleVals==90,1), Y(angleVals==90,2), Y(angleVals==90,3), 'k')
view(-135, 35);

%% weighting
cIdx = 10;
linear = allClassifiers{cIdx, nCells}.ClassificationDiscriminant.Coeffs(1,2).Linear;
weights = abs(linear);
[sorted, I] = sort(weights);
bestID = I(end-2:end);
worstID = I(1:3);
cellChoice = allCellChoice{cIdx, nCells};

figure; hold on;
scatter3(tableData(angleVals==45, cellChoice(bestID(1))), tableData(angleVals==45, cellChoice(bestID(2))),...
        tableData(angleVals==45, cellChoice(bestID(3))), 'r')
scatter3(tableData(angleVals==135, cellChoice(bestID(1))), tableData(angleVals==135, cellChoice(bestID(2))),...
        tableData(angleVals==135, cellChoice(bestID(3))), 'b')
view(-135, 35);

figure; hold on;
scatter3(tableData(angleVals==45, cellChoice(worstID(1))), tableData(angleVals==45, cellChoice(worstID(2))),...
         tableData(angleVals==45, cellChoice(worstID(3))), 'r')
scatter3(tableData(angleVals==135, cellChoice(worstID(1))), tableData(angleVals==135, cellChoice(worstID(2))),...
         tableData(angleVals==135, cellChoice(worstID(3))), 'b')
view(-135, 35);


