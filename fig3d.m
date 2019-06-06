%% fig3d Angle tuning example
%% fig3e Angle tuning classification and distribution

%% example from JK025 S04
baseDir = 'D:\TPM\JK\suite2p\';
mouse = 25;
session = 4;

cd(sprintf('%s%03d',baseDir, mouse))
load(sprintf('JK%03dS%02dangle_tuning',mouse,session))
load(sprintf('UberJK%03dS%02d',mouse,session))

plane = 3;

touchID = spk.touchID(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);
tunedAngle = spk.tunedAngle(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);

%% pick one from tunedAngle

tai = 39;
%
cid = touchID(tai);
example_angle_tuning_calcium(u, ca, spk, cid)

%%
cind = find(u.cellNums==cid);
figure, 
imshow(mat2gray(u.mimg{plane}))
hold on
plot(u.cellx(cind), u.celly(cind), 'y.','markersize', 20)


%% proportion of tuned cells
%% in L2/3 C2, L2/3 non-C2, L4 C2, L4 non-C2
clear
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]}; 

naiveInds = 1:12;
expertInds = [1:4,7,9];
scnnInds = 13:16;

for ni = 1 : length(naiveInds)
    mouse = mice(naiveInds(ni));
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{naiveInds(ni)}(1);
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
%         if ~strcmp(fieldnames{fi}, 'val')
            naive(ni).(fieldnames{fi}) = spk.(fieldnames{fi});
%         end
    end
end

for ei = 1 : length(expertInds)
    mouse = mice(expertInds(ei));
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{expertInds(ei)}(2);
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
        
        expert(ei).(fieldnames{fi}) = spk.(fieldnames{fi});
        
    end
end

% for si = 1 : length(scnnInds)
%     mouse = mice(scnnInds(si));
%     cd(sprintf('%s%03d',baseDir, mouse))
%     session = sessions{scnnInds(si)}(1);
%     load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
%     fieldnames = fields(spk);
%     for fi = 1 : length(fieldnames)
%         
%         L4(si).(fieldnames{fi}) = spk.(fieldnames{fi});
%         
%     end
% end

cd(baseDir)
save('angle_tuning_summary.mat','naive','expert')

%%
cd(baseDir)
load('angle_tuning_summary.mat','naive','expert','L4')


%% naive
tune = load('angle_tuning_summary.mat');
info = load('cellFunctionRidgeDE010');
touch = zeros(length(naive),4); % Both tuned and not. 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
tuned = zeros(length(naive),4);
for ni = 1 : length(naive)    
    temp = info.naive(ni);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    cellIDs = temp.cellNums;
    
    temp = tune.naive(ni);
    indTouch = find(ismember(cellIDs,temp.touchID));
    indTuned = find(ismember(cellIDs,temp.touchID(find(temp.tuned))));
    
    touch(ni,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touch(ni,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touch(ni,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touch(ni,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    tuned(ni,1) = length(intersect(indTuned, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    tuned(ni,2) = length(intersect(indTuned, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    tuned(ni,3) = length(intersect(indTuned, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    tuned(ni,4) = length(intersect(indTuned, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

figure, hold on
bar(1, mean(touch(:,1)), 'facecolor', [0.7 0.7 0.7])
bar(1, mean(tuned(:,1)), 'facecolor', [0.1 0.1 0.8])
for i = 1 : 4
    bar(i, mean(touch(:,i)), 'facecolor', [0.7 0.7 0.7])
    errorbar(i, mean(touch(:,i)), std(touch(:,i) - tuned(:,i))/sqrt(size(touch,1)), 'color', [0.7 0.7 0.7], 'linewidth', 2)
end
for i = 1 : 4
    bar(i, mean(tuned(:,i)), 'facecolor', [0.1 0.1 0.8])
    errorbar(i, mean(tuned(:,i)), std(tuned(:,i))/sqrt(size(touch,1)), 'color', [0.1 0.1 0.8], 'linewidth', 2)
end
legend({'Touch', 'Tuned'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%% expert
touch = zeros(length(expert),4); % Both tuned and not. 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
tuned = zeros(length(expert),4);
for ei = 1 : length(expert)    
    temp = info.expert(ei);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    cellIDs = temp.cellNums;
    
    temp = tune.expert(ei);
    indTouch = find(ismember(cellIDs,temp.touchID));
    indTuned = find(ismember(cellIDs,temp.touchID(find(temp.tuned))));
    
    touch(ei,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touch(ei,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touch(ei,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touch(ei,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    tuned(ei,1) = length(intersect(indTuned, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    tuned(ei,2) = length(intersect(indTuned, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    tuned(ei,3) = length(intersect(indTuned, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    tuned(ei,4) = length(intersect(indTuned, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

figure, hold on
bar(1, mean(touch(:,1)), 'facecolor', [0.7 0.7 0.7])
bar(1, mean(tuned(:,1)), 'facecolor', [0.1 0.1 0.8])
for i = 1 : 4
    bar(i, mean(touch(:,i)), 'facecolor', [0.7 0.7 0.7])
    errorbar(i, mean(touch(:,i)), std(touch(:,i) - tuned(:,i))/sqrt(size(touch,1)), 'color', [0.7 0.7 0.7], 'linewidth', 2)
end
for i = 1 : 4
    bar(i, mean(tuned(:,i)), 'facecolor', [0.1 0.1 0.8])
    errorbar(i, mean(tuned(:,i)), std(tuned(:,i))/sqrt(size(touch,1)), 'color', [0.1 0.1 0.8], 'linewidth', 2)
end
legend({'Touch', 'Tuned'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%% From scnn1a. last two mice
touch = zeros(2,4); % Both tuned and not. 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
tuned = zeros(2,4);
for si = 1 : 2    
    temp = info.expert(si+2);
    indL23 = find(temp.cellDepths < 350);
    indL4 = find(temp.cellDepths >= 350);
    indC2 = find(temp.isC2);
    indnonC2 = find(temp.isC2==0);
    cellIDs = temp.cellNums;
    
    temp = tune.expert(si+2);
    indTouch = find(ismember(cellIDs,temp.touchID));
    indTuned = find(ismember(cellIDs,temp.touchID(find(temp.tuned))));
    
    touch(si,1) = length(intersect(indTouch, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    touch(si,2) = length(intersect(indTouch, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    touch(si,3) = length(intersect(indTouch, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    touch(si,4) = length(intersect(indTouch, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
    
    tuned(si,1) = length(intersect(indTuned, intersect(indL23, indC2))) / length(intersect(indL23, indC2));
    tuned(si,2) = length(intersect(indTuned, intersect(indL23, indnonC2))) / length(intersect(indL23, indnonC2));
    tuned(si,3) = length(intersect(indTuned, intersect(indL4, indC2))) / length(intersect(indL4, indC2));
    tuned(si,4) = length(intersect(indTuned, intersect(indL4, indnonC2))) / length(intersect(indL4, indnonC2));
end

figure, hold on
bar(1, mean(touch(:,1)), 'facecolor', [0.7 0.7 0.7])
bar(1, mean(tuned(:,1)), 'facecolor', [0.1 0.1 0.8])
for i = 1 : 4
    bar(i, mean(touch(:,i)), 'facecolor', [0.7 0.7 0.7])
    errorbar(i, mean(touch(:,i)), std(touch(:,i) - tuned(:,i))/sqrt(size(touch,1)), 'color', [0.7 0.7 0.7], 'linewidth', 2)
end
for i = 1 : 4
    bar(i, mean(tuned(:,i)), 'facecolor', [0.1 0.1 0.8])
    errorbar(i, mean(tuned(:,i)), std(tuned(:,i))/sqrt(size(touch,1)), 'color', [0.1 0.1 0.8], 'linewidth', 2)
end
legend({'Touch', 'Tuned'}, 'box', 'off')
xticks([1:4])
xticklabels({'L2/3 C2', 'L2/3 non-C2', 'L4 C2', 'L4 non-C2'})
xtickangle(45)
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10)

%% Tiling all angles
clear
stretchFactor = 20;
load('angle_tuning_summary.mat');
temp = naive(1);
[~, indsortSharpness] = sort(temp.sharpness, 'descend');
tempTunedAngle = temp.tunedAngle(indsortSharpness);
[~, indsortAngle] = sort(tempTunedAngle);
indsort = indsortSharpness(indsortAngle);
numNotTuned = length(find(temp.tunedAngle==0));
tilingMap = zeros(length(temp.tunedAngle) - numNotTuned, 7 * stretchFactor);
for i = 1 : length(temp.tunedAngle) - numNotTuned
    tempMap = cellfun(@mean, temp.val{indsort(i+numNotTuned)});
    normTempMap = min_max_normalization(tempMap);
    for j = 1 : 7
        tilingMap(i,(j-1)*stretchFactor+1:j*stretchFactor) = deal(normTempMap(j));
    end
end
%%
figure
imshow(tilingMap), hold on, 
plot([0 0], [0 size(tilingMap,1)+1], 'k-', 'linewidth', 2)
plot([0 size(tilingMap,2)+1], [0 0], 'k-', 'linewidth', 2)
plot([size(tilingMap,2)+1 size(tilingMap,2)+1], [0 size(tilingMap,1)+1], 'k-', 'linewidth', 2)
plot([size(tilingMap,2)+1 0], [size(tilingMap,1)+1 size(tilingMap,1)+1], 'k-', 'linewidth', 2)
xlim([0 size(tilingMap,2)+1]), ylim([0 size(tilingMap,1)+1])

%% All tuned angles
%% naive
angles = 45:15:135;
tuning = zeros(length(naive), length(angles));
for i = 1 : length(naive)
    temp = naive(i);
    for j = 1 : length(angles)
        tuning(i,j) = length(find(temp.tunedAngle == angles(j))) / sum(temp.tuned);
    end
end

figure,
errorbar(angles, mean(tuning), std(tuning)/sqrt(length(naive)), 'k-', 'linewidth', 3)
xlabel('Object angle (\circ)')
xticks([45:15:135])
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10, 'box', 'off')
    
%% expert

angles = 45:15:135;
tuning = zeros(length(expert), length(angles));
for i = 1 : length(expert)
    temp = expert(i);
    for j = 1 : length(angles)
        tuning(i,j) = length(find(temp.tunedAngle == angles(j))) / sum(temp.tuned);
    end
end

figure,
errorbar(angles, mean(tuning), std(tuning)/sqrt(length(expert)), 'k-', 'linewidth', 3)
xlabel('Object angle (\circ)')
ylabel('Proportion')
set(gca, 'linewidth', 2, 'fontweight', 'bold', 'fontsize', 10, 'box', 'off')

%% different tuning patterns
inhibited = zeros(length(naive),1);
for i = 1 : length(naive)
    inhibited(i) = length(find(naive(i).tuneDirection > 1)) / sum(naive(i).tuned);
end
mean(inhibited)

multimodal = zeros(length(naive),1);
for i = 1 : length(naive)
    multimodal(i) = sum(naive(i).multimodal) / sum(naive(i).tuned);
end
mean(multimodal)

LOO = zeros(length(naive),1);
for i = 1 : length(naive)
    LOO(i) = sum(naive(i).leaveOneOut) / sum(naive(i).tuned);
end
mean(LOO)

categorical = zeros(length(naive),1);
for i = 1 : length(naive)
    categorical(i) = sum(naive(i).categorical) / sum(naive(i).tuned);
end
mean(categorical)

ramp = zeros(length(naive),1);
for i = 1 : length(naive)
    ramp(i) = sum(naive(i).ramp) / sum(naive(i).tuned);
end
mean(ramp)

single = zeros(length(naive),1);
for i = 1 : length(naive)
    single(i) = length(setdiff(find(naive(i).unimodalSingle), find(naive(i).multimodal))) / sum(naive(i).tuned);
end
mean(single)

broad = zeros(length(naive),1);
for i = 1 : length(naive)
    broad(i) = length(setdiff(find(naive(i).unimodalBroad), find(naive(i).multimodal))) / sum(naive(i).tuned);
end
mean(broad)


%% expert
inhibited = zeros(length(expert),1);
for i = 1 : length(expert)
    inhibited(i) = length(find(expert(i).tuneDirection > 1)) / sum(expert(i).tuned);
end
mean(inhibited)

multimodal = zeros(length(expert),1);
for i = 1 : length(expert)
    multimodal(i) = sum(expert(i).multimodal) / sum(expert(i).tuned);
end
mean(multimodal)

LOO = zeros(length(expert),1);
for i = 1 : length(expert)
    LOO(i) = sum(expert(i).leaveOneOut) / sum(expert(i).tuned);
end
mean(LOO)

categorical = zeros(length(expert),1);
for i = 1 : length(expert)
    categorical(i) = sum(expert(i).categorical) / sum(expert(i).tuned);
end
mean(categorical)

ramp = zeros(length(expert),1);
for i = 1 : length(expert)
    ramp(i) = sum(expert(i).ramp) / sum(expert(i).tuned);
end
mean(ramp)

single = zeros(length(naive),1);
for i = 1 : length(naive)
    single(i) = sum(naive(i).unimodalSingle) / sum(naive(i).tuned);
end
mean(single)

broad = zeros(length(naive),1);
for i = 1 : length(naive)
    broad(i) = sum(naive(i).unimodalBroad) / sum(naive(i).tuned);
end
mean(broad)


%% Examples
tempind = find(naive(3).leaveOneOut);

example_angle_tuning_calcium(u,ca,spk,spk.touchID(tempind(3)))

