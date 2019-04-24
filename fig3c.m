% for fig3c. Composition of touch cells and whisking cells. and mixed.
%     in L2/3 C2, L2/3 non-C2, L4 C2, and L4 non-C2. Pie graph

load('cellFunctionRidgeDE010')
%% all naive
touches = zeros(length(naive), 4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskings = zeros(length(naive), 4); 
mixed = zeros(length(naive), 4);
others = zeros(length(naive), 4);
for i = 1 : length(naive)
    temp = naive(i);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    
    indTouch = find(ismember(temp.cellNums, temp.touchID));
    indWhisking = find(ismember(temp.cellNums, temp.whiskingID));
    indMixed = intersect(indTouch, indWhisking);
    indOther = find(ismember(temp.cellNums, union(temp.soundID, union(temp.rewardID, temp.lickingID))));
    
    touches(i,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touches(i,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touches(i,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touches(i,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    whiskings(i,1) = length(intersect(indWhisking, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    whiskings(i,2) = length(intersect(indWhisking, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    whiskings(i,3) = length(intersect(indWhisking, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    whiskings(i,4) = length(intersect(indWhisking, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    mixed(i,1) = length(intersect(indMixed, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    mixed(i,2) = length(intersect(indMixed, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    mixed(i,3) = length(intersect(indMixed, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    mixed(i,4) = length(intersect(indMixed, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    others(i,1) = length(intersect(indOther, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    others(i,2) = length(intersect(indOther, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    others(i,3) = length(intersect(indOther, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    others(i,4) = length(intersect(indOther, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

%%
figure, hold all
% for legend ordering
i = 1;
touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
touchNwhisking = mean(mixed(:,i));
whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
other = mean(others(:,i));
bar(i, touchOnly, 'b')
bar(i, touchOnly + touchNwhisking, 'y')
bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    
for i = 1 : 4
%     subplot(2,2,i) % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
    touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
    touchNwhisking = mean(mixed(:,i));
    whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
    other = mean(others(:,i));
%     pie([touchOnly, touchNwhisking, whiskingOnly, other, 1-touchOnly - touchNwhisking - whiskingOnly - other])
    bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
    bar(i, touchOnly + touchNwhisking, 'y')
    bar(i, touchOnly, 'b')
end
legend({'Touch', 'Touch & Whisking', 'Whisking', 'Others'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%% expert
    
touches = zeros(length(expert), 4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskings = zeros(length(expert), 4); 
mixed = zeros(length(expert), 4);
others = zeros(length(expert), 4);
for i = 1 : length(expert)
    temp = expert(i);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    
    indTouch = find(ismember(temp.cellNums, temp.touchID));
    indWhisking = find(ismember(temp.cellNums, temp.whiskingID));
    indMixed = intersect(indTouch, indWhisking);
    indOther = find(ismember(temp.cellNums, union(temp.soundID, union(temp.rewardID, temp.lickingID))));
    
    touches(i,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touches(i,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touches(i,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touches(i,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    whiskings(i,1) = length(intersect(indWhisking, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    whiskings(i,2) = length(intersect(indWhisking, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    whiskings(i,3) = length(intersect(indWhisking, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    whiskings(i,4) = length(intersect(indWhisking, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    mixed(i,1) = length(intersect(indMixed, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    mixed(i,2) = length(intersect(indMixed, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    mixed(i,3) = length(intersect(indMixed, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    mixed(i,4) = length(intersect(indMixed, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    others(i,1) = length(intersect(indOther, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    others(i,2) = length(intersect(indOther, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    others(i,3) = length(intersect(indOther, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    others(i,4) = length(intersect(indOther, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

%%
figure, hold all
% for legend ordering
i = 1;
touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
touchNwhisking = mean(mixed(:,i));
whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
other = mean(others(:,i));
bar(i, touchOnly, 'b')
bar(i, touchOnly + touchNwhisking, 'y')
bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    
for i = 1 : 4
%     subplot(2,2,i) % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
    touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
    touchNwhisking = mean(mixed(:,i));
    whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
    other = mean(others(:,i));
%     pie([touchOnly, touchNwhisking, whiskingOnly, other, 1-touchOnly - touchNwhisking - whiskingOnly - other])
    bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
    bar(i, touchOnly + touchNwhisking, 'y')
    bar(i, touchOnly, 'b')
end
legend({'Touch', 'Touch & Whisking', 'Whisking', 'Others'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)


%% from Scnn1a

touches = zeros(2, 4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskings = zeros(2, 4); 
mixed = zeros(2, 4);
others = zeros(2, 4);
for i = 1 : 2
    temp = L4(i+2);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    
    indTouch = find(ismember(temp.cellNums, temp.touchID));
    indWhisking = find(ismember(temp.cellNums, temp.whiskingID));
    indMixed = intersect(indTouch, indWhisking);
    indOther = find(ismember(temp.cellNums, union(temp.soundID, union(temp.rewardID, temp.lickingID))));
    
    touches(i,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touches(i,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touches(i,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touches(i,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    whiskings(i,1) = length(intersect(indWhisking, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    whiskings(i,2) = length(intersect(indWhisking, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    whiskings(i,3) = length(intersect(indWhisking, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    whiskings(i,4) = length(intersect(indWhisking, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    mixed(i,1) = length(intersect(indMixed, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    mixed(i,2) = length(intersect(indMixed, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    mixed(i,3) = length(intersect(indMixed, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    mixed(i,4) = length(intersect(indMixed, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    others(i,1) = length(intersect(indOther, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    others(i,2) = length(intersect(indOther, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    others(i,3) = length(intersect(indOther, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    others(i,4) = length(intersect(indOther, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

%
figure, hold all
% for legend ordering
i = 1;
touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
touchNwhisking = mean(mixed(:,i));
whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
other = mean(others(:,i));
bar(i, touchOnly, 'b')
bar(i, touchOnly + touchNwhisking, 'y')
bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    
for i = 1 : 4
%     subplot(2,2,i) % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
    touchOnly = mean(touches(:,i)) - mean(mixed(:,i));
    touchNwhisking = mean(mixed(:,i));
    whiskingOnly = mean(whiskings(:,i)) - mean(mixed(:,i));
    other = mean(others(:,i));
%     pie([touchOnly, touchNwhisking, whiskingOnly, other, 1-touchOnly - touchNwhisking - whiskingOnly - other])
    bar(i, touchOnly + touchNwhisking + whiskingOnly + other, 'facecolor', [0.7 0.7 0.7])
    bar(i, touchOnly + touchNwhisking + whiskingOnly, 'g')
    bar(i, touchOnly + touchNwhisking, 'y')
    bar(i, touchOnly, 'b')
end
legend({'Touch', 'Touch & Whisking', 'Whisking', 'Others'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)
