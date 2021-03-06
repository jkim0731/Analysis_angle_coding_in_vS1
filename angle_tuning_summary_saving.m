%% from predecision touches

clear
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@length, sessions)==2);

for ni = 1 : length(mice)
    mouse = mice(ni);
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{ni}(1);
%     load(sprintf('JK%03dS%02dangle_tuning_lasso_predecision_NC',mouse,session), 'spk')
    load(sprintf('JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC',mouse,session), 'spk')
    load(sprintf('UberJK%03dS%02d_NC',mouse,session), 'u')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
        naive(ni).(fieldnames{fi}) = spk.(fieldnames{fi});
    end    
    naive(ni).depth = u.cellDepths(find(ismember(u.cellNums, spk.touchID))); % spk.,touchID is sorted
    naive(ni).isC2 = u.isC2(find(ismember(u.cellNums, spk.touchID)));
    naive(ni).noise = u.noise(find(ismember(u.cellNums, spk.touchID)));
    
    wrongInds = intersect(find(naive(ni).multimodal), find(naive(ni).unimodalBroad));
    naive(ni).multimodal(wrongInds) = 0;
end

for ei = 1 : length(expertInds)
    mouse = mice(expertInds(ei));
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{expertInds(ei)}(2);
%     load(sprintf('JK%03dS%02dangle_tuning_lasso_predecision_NC',mouse,session), 'spk')
    load(sprintf('JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC',mouse,session), 'spk')
    load(sprintf('UberJK%03dS%02d_NC',mouse,session), 'u')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
        expert(ei).(fieldnames{fi}) = spk.(fieldnames{fi});
    end
    expert(ei).depth = u.cellDepths(find(ismember(u.cellNums, spk.touchID))); % spk.,touchID is sorted
    expert(ei).isC2 = u.isC2(find(ismember(u.cellNums, spk.touchID)));
    expert(ei).noise = u.noise(find(ismember(u.cellNums, spk.touchID)));
    
    wrongInds = intersect(find(expert(ei).multimodal), find(expert(ei).unimodalBroad));
    expert(ei).multimodal(wrongInds) = 0;
end

cd(baseDir)
save('angle_tuning_summary_preAnswer_perTouch_NC.mat','naive','expert')





