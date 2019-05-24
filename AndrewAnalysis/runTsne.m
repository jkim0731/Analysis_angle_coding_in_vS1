%% t-sne
figure; hold on;
Y = tsne(tableData,'NumDimensions',3);

cmap = {[1,0,0],...
        [0.66, 0, 0],...
        [0.33, 0, 0],...
        [0,0,0],...
        [0,0,0.33],...
        [0,0,0.66],...
        [0,0,1]};

scatter3(Y(angleVals==45,1), Y(angleVals==45,2), Y(angleVals==45,3), 25, cmap{1}, 'filled')
scatter3(Y(angleVals==60,1), Y(angleVals==60,2), Y(angleVals==60,3), 25, cmap{2}, 'filled')
scatter3(Y(angleVals==75,1), Y(angleVals==75,2), Y(angleVals==75,3), 25, cmap{3}, 'filled')

scatter3(Y(angleVals==90,1), Y(angleVals==90,2), Y(angleVals==90,3), 25, cmap{4})

scatter3(Y(angleVals==105,1), Y(angleVals==105,2), Y(angleVals==105,3), 25, cmap{5}, 'filled')
scatter3(Y(angleVals==120,1), Y(angleVals==120,2), Y(angleVals==120,3), 25, cmap{6}, 'filled')
scatter3(Y(angleVals==135,1), Y(angleVals==135,2), Y(angleVals==135,3), 25, cmap{7}, 'filled')

view(-135, 35);

% distance t-sne space (change Y to tableData to look in response space)
% distance = squareform(pdist(Y));
distance = squareform(pdist(tableData));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);
figure;
imagesc(distance);
set(gca,'YDir','normal');
colormap('hot');

% angle space distance
angleDistance = squareform(pdist(angleVals));
figure; subplot(1,2,1);
scatter(angleDistance(:) + rand(length(angleDistance(:)), 1)*10, distance(:), 'k');
ylim([0, 1]);

subplot(1,2,2);
boxplot(distance(:), angleDistance(:))
ylim([0, 1]);