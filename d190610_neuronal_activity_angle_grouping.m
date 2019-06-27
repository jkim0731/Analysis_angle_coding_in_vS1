%%
% after running AndrewAnalysis\master.m

%%
tempind = 3;
figure, hold on
for i = 1 : size(expertPerformance,1)
    plot([1, 2], [naivePerformance(i,tempind), expertPerformance(i,tempind)], 'k-')
end
scatter(ones(size(expertPerformance,1),1), naivePerformance(:,tempind))
scatter(ones(size(expertPerformance,1),1)*2, expertPerformance(:,tempind))
xlim([0.5, 2.5])
title('L4 KNN')

%% test of normality
for i = 1 : 3
    kstest(naivePerformance(:,i))
    kstest(expertPerformance(:,i))
end

%% results: 
% everything is not normally distributed

%% signed rank test
for i = 1 : 3
    signrank(naivePerformance(:,i), expertPerformance(:,i))
end

%%
for i = 1 : 3
    [~, p] = ttest(naivePerformance(:,i), expertPerformance(:,i))
end