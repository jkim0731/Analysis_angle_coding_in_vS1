% shuffling error ratio VS exclusion error ratio

load('glm_cell_function_error_ratio')
%%
figure, hold on
for i = 1 : 12
    plot(cellfun(@(x) x(1), naive(i).exclusionER), cellfun(@(x) mean(x(1,:)), naive(i).errorRatio), 'k.'), axis equal
end
xlabel('Exclusion error ratio')
ylabel('Shuffling error ratio')
xlim([1 2.2])
ylim([1 2.2])
title('Naive touch')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

figure, hold on
for i = 1 : 12
    plot(cellfun(@(x) x(4), naive(i).exclusionER), cellfun(@(x) mean(x(4,:)), naive(i).errorRatio), 'k.'), axis equal
end
xlabel('Exclusion error ratio')
ylabel('Shuffling error ratio')
xlim([1 2.2])
ylim([1 2.2])
title('Naive whisking')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)
%%
figure, hold on
for i = 1 : 6
    tempInds = find(cellfun(@(x) ~isempty(x), expert(i).exclusionER));
    plot(cellfun(@(x) x(1), expert(i).exclusionER(tempInds)), cellfun(@(x) mean(x(1,:)), expert(i).errorRatio(tempInds)), 'k.'), axis equal
end
xlabel('Exclusion error ratio')
ylabel('Shuffling error ratio')
xlim([1 2.2])
ylim([1 2.2])
title('Expert touch')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%% DE diff vs shuffling error ratio

figure, hold on
for i = 1 : 12
    tempInds = find(cellfun(@(x) ~isempty(x), naive(i).exclusionER));
    plot(cellfun(@(x) x(1), naive(i).DEdiff(tempInds)), cellfun(@(x) mean(x(1,:)), naive(i).errorRatio(tempInds)), 'k.')
end
plot([0.1 0.1], [1 2.2], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
ylabel('Shuffling error ratio')
xlabel('Exclusion DE difference')
xlim([0 0.6])
ylim([1 2.2])
title('Naive touch')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%%
figure, hold on
for i = 1 : 12
    tempInds = find(cellfun(@(x) ~isempty(x), naive(i).exclusionER));
    plot(cellfun(@(x) x(4), naive(i).DEdiff(tempInds)), cellfun(@(x) mean(x(4,:)), naive(i).errorRatio(tempInds)), 'k.')
end
plot([0.1 0.1], [1 2.2], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
ylabel('Shuffling error ratio')
xlabel('Exclusion DE difference')
xlim([0 0.6])
ylim([1 2.2])
title('Naive whisking')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%%

figure, hold on
for i = 1 : 6
    tempInds = find(cellfun(@(x) ~isempty(x), expert(i).exclusionER));
    plot(cellfun(@(x) x(1), expert(i).DEdiff(tempInds)), cellfun(@(x) mean(x(1,:)), expert(i).errorRatio(tempInds)), 'k.')
end
plot([0.1 0.1], [1 2.2], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
ylabel('Shuffling error ratio')
xlabel('Exclusion DE difference')
xlim([0 0.6])
ylim([1 2.2])
title('Expert touch')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%%
figure, hold on
for i = 1 : 6
    tempInds = find(cellfun(@(x) ~isempty(x), expert(i).exclusionER));
    plot(cellfun(@(x) x(4), expert(i).DEdiff(tempInds)), cellfun(@(x) mean(x(4,:)), expert(i).errorRatio(tempInds)), 'k.')
end
plot([0.1 0.1], [1 2.2], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
ylabel('Shuffling error ratio')
xlabel('Exclusion DE difference')
xlim([0 0.6])
ylim([1 2.2])
title('Expert whisking')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)
