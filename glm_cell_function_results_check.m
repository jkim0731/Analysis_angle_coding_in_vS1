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
tic
baseDir = 'C:\JK\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);
%
naive = struct;
for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);
    [naive.glm(ni), naive.ca(ni), naive.spk(ni)] = glm_results_cell_function(mouse, session, baseDir);
end

expert = struct;
for ei = 1 : length(expertInd)    
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);
    [expert.glm(ei), expert.ca(ei), expert.spk(ei)] = glm_results_cell_function(mouse, session, baseDir);
end

L4mice = [70,74,75,76];
L4sessions = [6,4,4,4];
L4 = struct;
for mi = 1 : length(L4mice)
    mouse = L4mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    session = L4sessions(mi);
    [L4.glm(ei), L4.ca(ei), L4.spk(ei)] = glm_results_cell_function(mouse, session, baseDir);
end

save('Y:\Whiskernas\JK\suite2p\cellFunctionDE005.mat', 'naive', 'expert', 'L4')
toc
%% C2 L2/3 touch in 90 degrees. both naive and expert
% props = zeros(length(naive.glm),3); % 1 peronTouch, 2 whisking, 3 mixed
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

%%
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

%%
test = glm_results_cell_function(25, 4, baseDir);

%%

