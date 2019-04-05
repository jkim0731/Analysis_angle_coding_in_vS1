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
baseDir = 'C:\Data\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);
%
naive = struct;
for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    session = sessions{naiveInd(ni)}(1);
    [naive.glm(ni), naive.ca(ni), naive.spk(ni)] = glm_results_cell_function(mouse, session, baseDir);
end

expert = struct;
for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    session = sessions{expertInd(ei)}(2);
    [expert.glm(ei), expert.ca(ei), expert.spk(ei)] = glm_results_cell_function(mouse, session, baseDir);
end

L4mice = [70,74,75,76];
L4sessions = [6,4,4,4];
L4 = struct;
for mi = 1 : length(L4mice)
    mouse = L4mice(mi);
    session = L4sessions(mi);
    [L4.glm(ei), L4.ca(ei), L4.spk(ei)] = glm_results_cell_function(mouse, session, baseDir);
end
