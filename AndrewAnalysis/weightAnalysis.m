%% weight analysis
cIdx = 10;
linear = allClassifiers{cIdx, nCells}.ClassificationDiscriminant.Coeffs(1,2).Linear;
weights = abs(linear);
[sorted, I] = sort(weights);
bestID = I(end-2:end);
worstID = I(1:3);
cellChoice = allCellChoice{cIdx, nCells};

figure; hold on;
for i = 1:length(angles)
scatter3(tableData(angleVals==angles(i), cellChoice(bestID(1))), tableData(angleVals==angles(i), cellChoice(bestID(2))),...
         tableData(angleVals==angles(i), cellChoice(bestID(3))), 25, cmap{i}, 'filled')
end
view(-135, 35);

figure; hold on;
for i = 1:length(angles)
scatter3(tableData(angleVals==angles(i), cellChoice(worstID(1))), tableData(angleVals==angles(i), cellChoice(worstID(2))),...
         tableData(angleVals==angles(i), cellChoice(worstID(3))), 25, cmap{i}, 'filled')
end
view(-135, 35);