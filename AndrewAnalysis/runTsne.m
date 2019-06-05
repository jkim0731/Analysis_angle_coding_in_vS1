%% t-sne

% Y = tsne(tableData,'NumDimensions',3, 'standardize', true);
Y = tsne(tableData,'NumDimensions',3);
cmap = jet(7);

figure; hold on;
for i = 1 : length(angles)    
    scatter3(Y(angleVals==angles(i),1), Y(angleVals==angles(i),2), Y(angleVals==angles(i),3), 25, cmap(i,:), 'filled')
end

view(-135, 35);
title(miceName{animalID})


% distance t-sne space (change Y to tableData to look in response space)
distance = squareform(pdist(Y));
% distance = squareform(pdist(tableData));
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);
figure;
imagesc(distance);
set(gca,'YDir','normal');
colormap('hot');
title(miceName{animalID})
xticks(find(diff(angleVals))+0.5)
xticklabels({''})
yticks(find(diff(angleVals))+0.5)
yticklabels({''})
set(gca, 'tickdir', 'out', 'ticklength', [0.03 0], 'box', 'off')

%%
% angle space distance
angleDistance = squareform(pdist(angleVals));
figure; subplot(1,2,1);
scatter(angleDistance(:) + rand(length(angleDistance(:)), 1)*10, distance(:), 'k');
ylim([0, 1]);

subplot(1,2,2);
boxplot(distance(:), angleDistance(:))
ylim([0, 1]);
title(miceName{animalID})
