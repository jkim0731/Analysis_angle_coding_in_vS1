% looking at distribution of touch, whisking, and mixed cells in L4 of
% Scnn1a mice (only 2 mice available. JK075 and JK076)

baseDir = 'D:\TPM\JK\suite2p\';
allcellL4 = load([baseDir, 'glm_results_responseType'], 'L4');
tuningL4 = load([baseDir, 'glm_results_DEcomparison'], 'L4');

allcell = allcellL4.L4(3:4);
tuning = tuningL4.L4(3:4);

touchProp = zeros(2,2); % (:,1) for C2, (:,2) for non-C2
tunedProp = zeros(2,2); % (:,1) for C2, (:,2) for non-C2

for i = 1 : 2
    touchProp(i) = length(tuning(i).tuned) / length(allcell(i).allDE);
    tunedProp(i) = sum(tuning(i).tuned) / length(allcell(i).allDE);
end

%%


baseDir = 'D:\TPM\JK\suite2p\';
% L4 = struct;
L4(1) = glm_results_cell_function(75,4,baseDir);
L4(2) = glm_results_cell_function(76,4,baseDir);

save([baseDir, 'Scnn1a_ridgeDE010'], 'L4')

%%
touchProp = zeros(2,2); % (:,1) for C2, (:,2) for non-C2

for i = 1 : 2
    temp = intersect(L4(i).touchID, L4(i).cellNums(find(L4(i).isC2)));
    touchProp(i,1) = length(temp) / length(L4(i).cellNums);
    temp = intersect(L4(i).touchID, L4(i).cellNums(find(1-L4(i).isC2)));
    touchProp(i,2) = length(temp) / length(L4(i).cellNums);
end


%%

clear
baseDir = 'D:\TPM\JK\suite2p\';
mice = [75,76];
sessions = {[4],[4]};

for i = 1 : length(mice)
    mouse = mice(i);
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{i}(1);
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    load(sprintf('UberJK%03dS%02d',mouse,session), 'u')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
        L4tune(i).(fieldnames{fi}) = spk.(fieldnames{fi});
    end    
    L4tune(i).depth = u.cellDepths(find(ismember(u.cellNums, spk.touchID))); % spk.,touchID is sorted
end

cd(baseDir)
save('angle_tuning_Scnn1a.mat','L4tune')

%% from these 2 mice, draw L2/3 depth and L4 depth, both from C2 and non-C2
clear
baseDir = 'D:\TPM\JK\suite2p\';
load([baseDir, 'Scnn1a_ridgeDE010'], 'L4')
load('angle_tuning_Scnn1a.mat','L4tune')

touchProp = zeros(2,6); % (:,1) for C2 total, (:,2) for non-C2 total, (:,3) for C2 L2/3 depth, (:,4) for non-C2 L2/3 depth, (:,5) for C2 L4 depth, (:,6) for non-C2 L4 depth
tuneProp = zeros(2,6); % (:,1) for C2 total, (:,2) for non-C2 total, (:,3) for C2 L2/3 depth, (:,4) for non-C2 L2/3 depth, (:,5) for C2 L4 depth, (:,6) for non-C2 L4 depth
for i = 1 : 2
    temp = intersect(L4(i).touchID, L4(i).cellNums(find(L4(i).isC2)));
    touchProp(i,1) = length(temp) / length(L4(i).cellNums);
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,1) = length(temp2) / length(L4(i).cellNums);

    temp = intersect(L4(i).touchID, L4(i).cellNums(find(1-L4(i).isC2)));
    touchProp(i,2) = length(temp) / length(L4(i).cellNums);
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,2) = length(temp2) / length(L4(i).cellNums);
    
    temp = intersect(intersect(L4(i).touchID, L4(i).cellNums(find(L4(i).isC2))), L4(i).cellFitID(L4(i).cellFitIndL23));
    touchProp(i,3) = length(temp) / length(find(L4(i).cellDepths < 350));    
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,3) = length(temp2) / length(find(L4(i).cellDepths < 350));
    
    temp = intersect(intersect(L4(i).touchID, L4(i).cellNums(find(1-L4(i).isC2))), L4(i).cellFitID(L4(i).cellFitIndL23));
    touchProp(i,4) = length(temp) / length(find(L4(i).cellDepths < 350));
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,4) = length(temp2) / length(find(L4(i).cellDepths < 350));
    
    temp = setdiff(intersect(L4(i).touchID, L4(i).cellNums(find(L4(i).isC2))), L4(i).cellFitID(L4(i).cellFitIndL23));
    touchProp(i,5) = length(temp) / length(find(L4(i).cellDepths >= 350));
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,5) = length(temp2) / length(find(L4(i).cellDepths >= 350));
    
    temp = setdiff(intersect(L4(i).touchID, L4(i).cellNums(find(1-L4(i).isC2))), L4(i).cellFitID(L4(i).cellFitIndL23));
    touchProp(i,6) = length(temp) / length(find(L4(i).cellDepths >= 350));
    temp2 = intersect(temp, L4tune(i).touchID(find(L4tune(i).tuned)));
    tuneProp(i,6) = length(temp2) / length(find(L4(i).cellDepths >= 350));
end

%% draw figure

figure, hold on
for i = 1 : 3
    plot(ones(1,2) * (i-0.15), touchProp(:,(i-1)*2+1), '.', 'color', [0 0 1], 'markersize', 20)
    plot(ones(1,2) * (i+0.15), touchProp(:,i*2), '.', 'color', [0 0 1], 'markersize', 20)
    plot(ones(1,2) * (i-0.15), tuneProp(:,(i-1)*2+1), '.', 'color',[0.7 0.7 0.7], 'markersize', 20)
    plot(ones(1,2) * (i+0.15), tuneProp(:,i*2), '.', 'color',[0.7 0.7 0.7], 'markersize', 20)
end

xticks([0.85, 1.15, 1.85, 2.15, 2.85, 3.15])
xticklabels({'All C2', 'All non-C2','L3 depth C2', 'L3 depth non-C2', 'L4 depth C2', 'L4 depth non-C2'})
xtickangle(45)
ylabel('Proportion')