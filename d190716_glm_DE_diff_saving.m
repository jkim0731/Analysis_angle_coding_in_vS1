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
    savefn = sprintf('glm_DE_cellFunction_JK%03dS%02d',mouse, session);
    glmDEfunction = glm_results_dev_exp_power(mouse, session, baseDir);
    save(savefn, 'glmDEfunction')
    naive(ni) = glmDEfunction;
end

for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
    savefn = sprintf('glm_DE_cellFunction_JK%03dS%02d',mouse, session);
    glmDEfunction = glm_results_dev_exp_power(mouse, session, baseDir);    
    save(savefn, 'glmDEfunction')
    expert(ei) = glmDEfunction;
end

% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
% 
% 
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_dev_exp_power(mouse, session, baseDir);
% end

% save('Y:\Whiskernas\JK\suite2p\glm_results_responseType.mat', 'naive', 'expert', 'L4')
save('Y:\Whiskernas\JK\suite2p\glm_results_responseType.mat', 'naive', 'expert')
toc % takes about 8 min in 36 core 
%%
% clear
% tic
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
% 
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 
% 
% naiveInd = 1:length(mice);
% expertInd = find(cellfun(@length, sessions)==2);
% 
% 
% for ni = 1 : length(naiveInd)
%     mouse = mice(naiveInd(ni));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{naiveInd(ni)}(1);    
%     savefn = sprintf('glm_DE_whisker_touchCell_JK%03dS%02d',mouse, session);
%     glmDEwhiskerTouchCell = glm_results_dev_exp_power(mouse, session, baseDir);    
%     save(savefn,'glmDEwhiskerTouchCell')
%     naive(ni) = glmDEwhiskerTouchCell;
% end
% 
% % expert = struct;
% for ei = 1 : length(expertInd)
%     mouse = mice(expertInd(ei));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{expertInd(ei)}(2);    
%     savefn = sprintf('glm_DE_whisker_touchCell_JK%03dS%02d',mouse, session);
%     glmDEwhiskerTouchCell = glm_results_dev_exp_power(mouse, session, baseDir);
%     save(savefn,'glmDEwhiskerTouchCell')
%     expert(ei) = glmDEwhiskerTouchCell;
% end
% 
% % L4mice = [70,74,75,76];
% % L4sessions = [6,4,4,4];
% % % L4mice = [70];
% % % L4sessions = [6];
% % 
% % % L4 = struct;
% % for mi = 1 : length(L4mice)
% %     mouse = L4mice(mi);
% %     cd(sprintf('%s%03d',baseDir,mouse))
% %     session = L4sessions(mi);    
% %     L4(mi) = glm_results_dev_exp_power(mouse, session, baseDir);
% % end
% 
% % save('Y:\Whiskernas\JK\suite2p\glm_results_withWTV_expert.mat', 'naive', 'expert', 'L4')
% save('Y:\Whiskernas\JK\suite2p\glm_results_WKV.mat', 'naive', 'expert')
% toc % takes about 8 min in 36 core
