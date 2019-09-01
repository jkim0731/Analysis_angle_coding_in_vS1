%% First, change suite2p files (only the ones containing 'final' has C2 info)
clear
% mice = [25,27,30,36,37,39,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  

% mice = [70, 74,75,76];
% sessions = {[4], [4],[4],[4]};  
mice = [27,36,41,52];
sessions = {[3,9,10,16,17],[1,17,18],[3],[3,4,21,22,26,27]};  

% baseD = 'Y:\Whiskernas\JK\suite2p\';
baseD = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseD,mouse));
    load(sprintf('JK%03dC2.mat',mouse))
%     flist = dir('F_0*_final.mat');
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        flist = dir(sprintf('F_%03d_%03d_plane*_proc_final*.mat', mouse, session));    
        for fi = 1 : length(flist)
            fprintf('JK%03d S%02d plane%d processing\n', mouse, session, fi)
            load(flist(fi).name)
            inds = find([dat.stat.iscell]);
            dat.isC2 = zeros(length(inds),1);
            for i = 1 : length(inds)
                if inpolygon(dat.ops.useY(dat.ops.yrange(round(dat.stat(inds(i)).med(1)))), dat.ops.useX(dat.ops.xrange(round(dat.stat(inds(i)).med(2)))), ypoints, xpoints)
                    dat.isC2(i) = 1;
                end
            end
            dat.c2ypoints = ypoints;
            dat.c2xpoints = xpoints;
            save(flist(fi).name,'dat', '-append')
        end
    end
end

%% change all uber files (re-build uber files)
mice = [27,36,41,52];
% sessions = {[3,9,10,16,17],[1,17,18],[3],[3,4,21,22,26,27]};  
% Many of these sessions do not have noise-correction.
sessions = {[3,10],[1,17],[3],[3,21]};  

for mi = 1 : length(mice)
    for si = 1 : length(sessions{mi})
        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));
    end
end

%% touch glm results 

glm_results_saving_touch

%% fix angle_tuning_lasso_predecision 
% only change the c2y and c2x points in 'info'
mice = [27,36,41,52];
sessions = {[3,10],[1,17],[3],[3,21]}; 
baseD = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseD,mouse));
    load(sprintf('JK%03dC2.mat',mouse))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        fn = sprintf('JK%03dS%02dangle_tuning_lasso_predecision_NC.mat', mouse, session);
        load(fn,'info')
        info.c2xpoints = xpoints;
        info.c2ypoints = ypoints;
        save(fn,'info', '-append')        
    end
end

