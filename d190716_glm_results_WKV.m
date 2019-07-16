mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

baseDir =  'Y:\Whiskernas\JK\suite2p\';

for mi =  1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        savefn = sprintf('glmResult_WKV_JK%03dS%02d',mouse, session);
        glmWKV = glm_results_WKV_exclusion(mouse, session, baseDir);
        save([baseDir,savefn], 'glmWKV')
    end
end