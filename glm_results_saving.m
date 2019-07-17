baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir, mouse))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        savefn = sprintf('JK%03dS%02dglm_cell_function_lasso',mouse, session);
        glm = glm_results_cell_function(mouse, session, baseDir);
        save(savefn, 'glm')
    end
end
