baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% % sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

% mice = [27,36,41,52];
% sessions = {[3,9,10,16,17],[1,17,18],[3],[3,4,21,22,26,27]};  
% Many of these sessions do not have noise-correction.
% sessions = {[3,10],[1,17],[3],[3,21]}; 

for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir, mouse))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        savefn = sprintf('JK%03dS%02dglm_cell_function_lasso_NC',mouse, session);
        glm = glm_results_cell_function_touch(mouse, session, baseDir);
        save(savefn, 'glm')
    end
end

%%
clear
tic
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);


for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);
    load(sprintf('JK%03dS%02dglm_cell_function_lasso_NC',mouse,session))
    naive(ni) = glm;
end

for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
%     expert(ei) = glm_results_cell_function(mouse, session, baseDir);
    load(sprintf('JK%03dS%02dglm_cell_function_lasso_NC',mouse,session))
    expert(ei) = glm;
end
% 
% % L4mice = [70,74,75,76];
% % L4sessions = [6,4,4,4];
% L4mice = [70];
% L4sessions = [6];
% 
% % L4 = struct;
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_cell_function(mouse, session, baseDir);
% end

save('Y:\Whiskernas\JK\suite2p\cellFunctionLasso_NC.mat', 'naive', 'expert')
toc