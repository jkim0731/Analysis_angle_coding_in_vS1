%%
% Check proportion of 90 degrees pole response touch cells and whisking
% cells, in or outside of C2, in L2/3. Compare with peron et al.

% Check consistency between glm, calcium anova, and spike anova data.

% See proportion of cell functions

% Compare between naive and expert


%%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3], [4], [4], [4], [4]}; 
% for mi = 1 : length(mice)
%     for si = 1 : length(sessions{mi})
%         Uber.buildUberArray(mice(mi), sessions{mi}(si))
%     end
% end



%%
clear
tic
baseDir = 'C:\JK\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);


for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);    
    naive(ni) = glm_results_cell_function(mouse, session, baseDir);
end

% expert = struct;
for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
    expert(ei) = glm_results_cell_function(mouse, session, baseDir);
end

% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
L4mice = [70];
L4sessions = [6];

% L4 = struct;
for mi = 1 : length(L4mice)
    mouse = L4mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    session = L4sessions(mi);    
    L4(mi) = glm_results_cell_function(mouse, session, baseDir);
end

save('Y:\Whiskernas\JK\suite2p\cellFunctionRidgeDE005.mat', 'naive', 'expert', 'L4')
toc
%% C2 L2/3 touch in 90 degrees. both naive and expert
%% All naive
props = zeros(length(naive.glm),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(naive.glm)
    ind = i;
    C2L23Cells = naive.glm(ind).cellNums(intersect( find(naive.glm(ind).isC2), find(naive.glm(ind).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(naive.glm(ind).touchID, naive.glm(ind).tunedID(naive.glm(ind).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(naive.glm(ind).cellFitID(find(cellfun(@(x) x(4), naive.glm(ind).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)

%% Matching naive
props = zeros(length(expertInd),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(expertInd)
    ind = expertInd(i);
    C2L23Cells = naive.glm(ind).cellNums(intersect( find(naive.glm(ind).isC2), find(naive.glm(ind).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(naive.glm(ind).touchID, naive.glm(ind).tunedID(naive.glm(ind).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(naive.glm(ind).cellFitID(find(cellfun(@(x) x(4), naive.glm(ind).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)

%% experts
props = zeros(length(expert.glm),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(expert.glm)
    C2L23Cells = expert.glm(i).cellNums(intersect( find(expert.glm(i).isC2), find(expert.glm(i).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(expert.glm(i).touchID, expert.glm(i).tunedID(expert.glm(i).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(expert.glm(i).cellFitID(find(cellfun(@(x) x(4), expert.glm(i).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)


%% Non-C2 L2/3 touch in 90 degrees. both naive and expert
%% All naive
props = zeros(length(naive.glm),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(naive.glm)
    ind = i;
    C2L23Cells = naive.glm(ind).cellNums(intersect( find(1-naive.glm(ind).isC2), find(naive.glm(ind).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(naive.glm(ind).touchID, naive.glm(ind).tunedID(naive.glm(ind).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(naive.glm(ind).cellFitID(find(cellfun(@(x) x(4), naive.glm(ind).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)

%% Matching naive
props = zeros(length(expertInd),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(expertInd)
    ind = expertInd(i);
    C2L23Cells = naive.glm(ind).cellNums(intersect( find(1-naive.glm(ind).isC2), find(naive.glm(ind).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(naive.glm(ind).touchID, naive.glm(ind).tunedID(naive.glm(ind).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(naive.glm(ind).cellFitID(find(cellfun(@(x) x(4), naive.glm(ind).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)
%% Experts
props = zeros(length(expert.glm),3); % 1 peronTouch, 2 whisking, 3 mixed
for i = 1 : length(expert.glm)
    C2L23Cells = expert.glm(i).cellNums(intersect( find(1-expert.glm(i).isC2), find(expert.glm(i).cellDepths<350)));
    peronC2L23TouchCells = intersect(union(expert.glm(i).touchID, expert.glm(i).tunedID(expert.glm(i).tunedAngle==90)), C2L23Cells );
    C2L23WhiskingCells = intersect(expert.glm(i).cellFitID(find(cellfun(@(x) x(4), expert.glm(i).cellFunction))), C2L23Cells);
    C2L23MixedCells = intersect(peronC2L23TouchCells, C2L23WhiskingCells);
    prop(i,1) = length(peronC2L23TouchCells)/length(C2L23Cells);
    prop(i,2) = length(C2L23WhiskingCells)/length(C2L23Cells);
    prop(i,3) = length(C2L23MixedCells)/length(C2L23Cells);
end
mean(prop)


%%
test = glm_results_cell_function(25, 4, baseDir);

%% Comparing between glm and anovas
% glm lasso DE 0.1
% 1) how many of touch cells and tuned cells from GLM are in Calcium ANOVA
% and spike ANOVA?
% 2) Are tuned angles agree with each other?

load('cellFunctionLassoDE010.mat')
%% all naive
tunedAgreement = zeros(length(naive.glm),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
touchAgreement = zeros(length(naive.glm),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
tunedPropInTuned = zeros(length(naive.glm),2); % 1 glm in ca, 2 glm in spk
tunedPropInTouch = zeros(length(naive.glm),2); % 1 glm in ca, 2 glm in spk
touchPropInTuned = zeros(length(naive.glm),2); % 1 glm in ca, 2 glm in spk
touchPropInTouch = zeros(length(naive.glm),2); % 1 glm in ca, 2 glm in spk

for i = 1 : length(naive.glm)
    tunedAgreement(i,1) = length(intersect(naive.glm(i).tunedID, naive.ca(i).tunedID)) / length(union(naive.glm(i).tunedID, naive.ca(i).tunedID));
    tunedAgreement(i,2) = length(intersect(naive.glm(i).tunedID, naive.spk(i).tunedID)) / length(union(naive.glm(i).tunedID, naive.spk(i).tunedID));
    tunedAgreement(i,3) = length(intersect(naive.ca(i).tunedID, naive.spk(i).tunedID)) / length(union(naive.ca(i).tunedID, naive.spk(i).tunedID));

    tunedPropInTuned(i,1) = length(find(ismember(naive.glm(i).tunedID, naive.ca(i).tunedID))) / length(naive.glm(i).tunedID);
    tunedPropInTuned(i,2) = length(find(ismember(naive.glm(i).tunedID, naive.spk(i).tunedID))) / length(naive.glm(i).tunedID);

    tunedPropInTouch(i,1) = length(find(ismember(naive.glm(i).tunedID, naive.ca(i).touch))) / length(naive.glm(i).tunedID);
    tunedPropInTouch(i,2) = length(find(ismember(naive.glm(i).tunedID, naive.spk(i).touch))) / length(naive.glm(i).tunedID);


    touchAgreement(i,1) = length(intersect(naive.glm(i).touchID, naive.ca(i).touch)) / length(union(naive.glm(i).touchID, naive.ca(i).touch));
    touchAgreement(i,2) = length(intersect(naive.glm(i).touchID, naive.spk(i).touch)) / length(union(naive.glm(i).touchID, naive.spk(i).touch));
    touchAgreement(i,3) = length(intersect(naive.ca(i).touch, naive.spk(i).touch)) / length(union(naive.ca(i).touch, naive.spk(i).touch));

    touchPropInTouch(i,1) = length(find(ismember(naive.glm(i).touchID, naive.ca(i).touch))) / length(naive.glm(i).touchID);
    touchPropInTouch(i,2) = length(find(ismember(naive.glm(i).touchID, naive.spk(i).touch))) / length(naive.glm(i).touchID);

    touchPropInTuned(i,1) = length(find(ismember(naive.glm(i).touchID, naive.ca(i).tunedID))) / length(naive.glm(i).touchID);
    touchPropInTuned(i,2) = length(find(ismember(naive.glm(i).touchID, naive.spk(i).tunedID))) / length(naive.glm(i).touchID);
end

mean(tunedAgreement)
mean(touchAgreement)
mean(tunedPropInTuned)
mean(tunedPropInTouch)
mean(touchPropInTouch)
mean(touchPropInTuned)

%% Experts
tunedAgreement = zeros(length(expert.glm),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
touchAgreement = zeros(length(expert.glm),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
tunedPropInTuned = zeros(length(expert.glm),2); % 1 glm in ca, 2 glm in spk
tunedPropInTouch = zeros(length(expert.glm),2); % 1 glm in ca, 2 glm in spk
touchPropInTuned = zeros(length(expert.glm),2); % 1 glm in ca, 2 glm in spk
touchPropInTouch = zeros(length(expert.glm),2); % 1 glm in ca, 2 glm in spk

for i = 1 : length(expert.glm)
    tunedAgreement(i,1) = length(intersect(expert.glm(i).tunedID, expert.ca(i).tunedID)) / length(union(expert.glm(i).tunedID, expert.ca(i).tunedID));
    tunedAgreement(i,2) = length(intersect(expert.glm(i).tunedID, expert.spk(i).tunedID)) / length(union(expert.glm(i).tunedID, expert.spk(i).tunedID));
    tunedAgreement(i,3) = length(intersect(expert.ca(i).tunedID, expert.spk(i).tunedID)) / length(union(expert.ca(i).tunedID, expert.spk(i).tunedID));

    tunedPropInTuned(i,1) = length(find(ismember(expert.glm(i).tunedID, expert.ca(i).tunedID))) / length(expert.glm(i).tunedID);
    tunedPropInTuned(i,2) = length(find(ismember(expert.glm(i).tunedID, expert.spk(i).tunedID))) / length(expert.glm(i).tunedID);

    tunedPropInTouch(i,1) = length(find(ismember(expert.glm(i).tunedID, expert.ca(i).touch))) / length(expert.glm(i).tunedID);
    tunedPropInTouch(i,2) = length(find(ismember(expert.glm(i).tunedID, expert.spk(i).touch))) / length(expert.glm(i).tunedID);


    touchAgreement(i,1) = length(intersect(expert.glm(i).touchID, expert.ca(i).touch)) / length(union(expert.glm(i).touchID, expert.ca(i).touch));
    touchAgreement(i,2) = length(intersect(expert.glm(i).touchID, expert.spk(i).touch)) / length(union(expert.glm(i).touchID, expert.spk(i).touch));
    touchAgreement(i,3) = length(intersect(expert.ca(i).touch, expert.spk(i).touch)) / length(union(expert.ca(i).touch, expert.spk(i).touch));

    touchPropInTouch(i,1) = length(find(ismember(expert.glm(i).touchID, expert.ca(i).touch))) / length(expert.glm(i).touchID);
    touchPropInTouch(i,2) = length(find(ismember(expert.glm(i).touchID, expert.spk(i).touch))) / length(expert.glm(i).touchID);

    touchPropInTuned(i,1) = length(find(ismember(expert.glm(i).touchID, expert.ca(i).tunedID))) / length(expert.glm(i).touchID);
    touchPropInTuned(i,2) = length(find(ismember(expert.glm(i).touchID, expert.spk(i).tunedID))) / length(expert.glm(i).touchID);
end

mean(tunedAgreement)
mean(touchAgreement)
mean(tunedPropInTuned)
mean(tunedPropInTouch)
mean(touchPropInTouch)
mean(touchPropInTuned)


%% matched naive
tunedAgreement = zeros(length(expertInd),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
touchAgreement = zeros(length(expertInd),3); % 1 glm-ca, 2 glm-spk, 3 ca-spk
tunedPropInTuned = zeros(length(expertInd),2); % 1 glm in ca, 2 glm in spk
tunedPropInTouch = zeros(length(expertInd),2); % 1 glm in ca, 2 glm in spk
touchPropInTuned = zeros(length(expertInd),2); % 1 glm in ca, 2 glm in spk
touchPropInTouch = zeros(length(expertInd),2); % 1 glm in ca, 2 glm in spk

for ind = 1 : length(expertInd)
    i = expertInd(ind);
    tunedAgreement(ind,1) = length(intersect(naive.glm(i).tunedID, naive.ca(i).tunedID)) / length(union(naive.glm(i).tunedID, naive.ca(i).tunedID));
    tunedAgreement(ind,2) = length(intersect(naive.glm(i).tunedID, naive.spk(i).tunedID)) / length(union(naive.glm(i).tunedID, naive.spk(i).tunedID));
    tunedAgreement(ind,3) = length(intersect(naive.ca(i).tunedID, naive.spk(i).tunedID)) / length(union(naive.ca(i).tunedID, naive.spk(i).tunedID));

    tunedPropInTuned(ind,1) = length(find(ismember(naive.glm(i).tunedID, naive.ca(i).tunedID))) / length(naive.glm(i).tunedID);
    tunedPropInTuned(ind,2) = length(find(ismember(naive.glm(i).tunedID, naive.spk(i).tunedID))) / length(naive.glm(i).tunedID);

    tunedPropInTouch(ind,1) = length(find(ismember(naive.glm(i).tunedID, naive.ca(i).touch))) / length(naive.glm(i).tunedID);
    tunedPropInTouch(ind,2) = length(find(ismember(naive.glm(i).tunedID, naive.spk(i).touch))) / length(naive.glm(i).tunedID);


    touchAgreement(ind,1) = length(intersect(naive.glm(i).touchID, naive.ca(i).touch)) / length(union(naive.glm(i).touchID, naive.ca(i).touch));
    touchAgreement(ind,2) = length(intersect(naive.glm(i).touchID, naive.spk(i).touch)) / length(union(naive.glm(i).touchID, naive.spk(i).touch));
    touchAgreement(ind,3) = length(intersect(naive.ca(i).touch, naive.spk(i).touch)) / length(union(naive.ca(i).touch, naive.spk(i).touch));

    touchPropInTouch(ind,1) = length(find(ismember(naive.glm(i).touchID, naive.ca(i).touch))) / length(naive.glm(i).touchID);
    touchPropInTouch(ind,2) = length(find(ismember(naive.glm(i).touchID, naive.spk(i).touch))) / length(naive.glm(i).touchID);

    touchPropInTuned(ind,1) = length(find(ismember(naive.glm(i).touchID, naive.ca(i).tunedID))) / length(naive.glm(i).touchID);
    touchPropInTuned(ind,2) = length(find(ismember(naive.glm(i).touchID, naive.spk(i).tunedID))) / length(naive.glm(i).touchID);
end

mean(tunedAgreement)
mean(touchAgreement)
mean(tunedPropInTuned)
mean(tunedPropInTouch)
mean(touchPropInTouch)
mean(touchPropInTuned)