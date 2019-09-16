clear
tic
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);

%%
% for ni = 1 : length(naiveInd)
  for ni = 6  
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);    
    naive(ni) = glm_results_dev_exp_wkv_touchCell(mouse, session, baseDir);
end
%%
% expert = struct;
for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
    expert(ei) = glm_results_dev_exp_wkv_touchCell(mouse, session, baseDir);
end

% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
% % L4mice = [70];
% % L4sessions = [6];
% 
% % L4 = struct;
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_dev_exp_power(mouse, session, baseDir);
% end

% save('Y:\Whiskernas\JK\suite2p\glm_results_withWTV_expert.mat', 'naive', 'expert', 'L4')
save('Y:\Whiskernas\JK\suite2p\glmResults_devExp_WKV_touchCell_NC.mat', 'naive', 'expert')
toc