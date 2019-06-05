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
        dIdx = find(ismember(angleVals, angleDiscrim));
%         dIdx = find(angleVals == angleDiscrim(1) | angleVals == angleDiscrim(2));
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
shadedErrorBar(1:nCells, mean(accuracy), std(accuracy), 'lineprops','b')
shadedErrorBar(1:nCells, mean(accuracyShuffled), std(accuracyShuffled), 'lineprops','r')
ylim([0 1]);