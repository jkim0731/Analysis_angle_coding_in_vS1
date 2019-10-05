mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

baseDir =  'Y:\Whiskernas\JK\suite2p\';

for mi =  1 : length(mice)
% for mi =  1 : 8
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        savefn = sprintf('glmResult_WKV_touchCell_NC_JK%03dS%02d',mouse, session);
        glmWKV = glm_results_WKV_exclusion(mouse, session, baseDir);
        save(sprintf('%s%03d\\%s',baseDir,mouse,savefn), 'glmWKV')
    end
end

%%

clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';
savefn = 'glmResults_WKV_touchCell_exclusion_NC';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

for mi = 1 : length(mice)
    load(sprintf('%s%03d\\glmResult_WKV_NC_JK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)))
    naive(mi) = glmWKV;
end

expertInd = find(cellfun(@(x) length(x) == 2, sessions));
for ei = 1 : length(expertInd)
    load(sprintf('%s%03d\\glmResult_WKV_NC_JK%03dS%02d', baseDir, mice(expertInd(ei)), mice(expertInd(ei)), sessions{expertInd(ei)}(2)))
    expert(ei) = glmWKV;
end

save(sprintf('%s%s',baseDir,savefn), 'naive', 'expert')