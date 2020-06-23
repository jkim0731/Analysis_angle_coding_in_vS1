% baseDir = 'Y:\Whiskernas\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
% 
% saveFn = 'modelAngleTuning_NC';
% loadFnBase = 'angle_tuning_model_touchCell_NC_preAnswer_perTouch_';
% for mi = 1 : length(mice)
%     mouse = mice(mi);
%     session = sessions{mi}(1);
%     naive(mi) = load(sprintf('%s%03d\\%sJK%03dS%02d', baseDir, mouse, loadFnBase, mouse, session));
% end
% 
% expertInds = find(cellfun(@length, sessions)==2);
% for ei = 1 : length(expertInds)
%     mi = expertInds(ei);
%     mouse = mice(mi);
%     session = sessions{mi}(2);
%     expert(ei) = load(sprintf('%s%03d\\%sJK%03dS%02d', baseDir, mouse, loadFnBase, mouse, session));
% end
% %%
% save(sprintf('%s%s',baseDir, saveFn), 'naive', 'expert')
% 
% 

%% Second round

% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

saveFn = 'modelAngleTuning_NC_combinations';
loadFnBase = 'angle_tuning_model_touchCell_NC_preAnswer_perTouch_';
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    naive(mi) = load(sprintf('%s%03d\\%sJK%03dS%02d_2ndRound', baseDir, mouse, loadFnBase, mouse, session));
end

expertInds = find(cellfun(@length, sessions)==2);
for ei = 1 : length(expertInds)
    mi = expertInds(ei);
    mouse = mice(mi);
    session = sessions{mi}(2);
    expert(ei) = load(sprintf('%s%03d\\%sJK%03dS%02d_2ndRound', baseDir, mouse, loadFnBase, mouse, session));
end

save(sprintf('%s%s',baseDir, saveFn), 'naive', 'expert')



%%

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

saveFn = 'modelAngleTuning_NC_touchOnly';
loadFnBase = 'angle_tuning_model_touchCell_NC_preAnswer_perTouch_';
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    naive(mi) = load(sprintf('%s%03d\\%sJK%03dS%02d_4thRound', baseDir, mouse, loadFnBase, mouse, session));
end

expertInds = find(cellfun(@length, sessions)==2);
for ei = 1 : length(expertInds)
    mi = expertInds(ei);
    mouse = mice(mi);
    session = sessions{mi}(2);
    expert(ei) = load(sprintf('%s%03d\\%sJK%03dS%02d_4thRound', baseDir, mouse, loadFnBase, mouse, session));
end

save(sprintf('%s%s',baseDir, saveFn), 'naive', 'expert')