clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
baseDir = 'C:\Data\suite2p\';
response90 = zeros(length(mice),1);

for mi = 1 : length(mice)
    mouse = mice(mi);
%     for si = 1 : length(sessions{mi})
    for si = 1
        session = sessions{mi}(si);
        
        resultsfn = sprintf('JK%03dS%02dsingleCell_anova_FT_1secB_absthresh_perm_sharp.mat', mouse, session);

        targetDir = sprintf('%s%03d',baseDir,mouse);
%         ufn = sprintf('UberJK%03dS%02d.mat',mouse,session);

        cd(targetDir)
        load(resultsfn)
%         load(ufn)

        clear cell
        
        resp = 0;
        for ci = 1 : length(cellsTotal)
            dF = dFtotal{ci}{4}; % 90 degrees
            temp = nanmean(dF(:,baseFrameNum+1:end),2);
            if mean(temp) > 0.2 && ttest(temp)
                resp = resp + 1;
            end
        end
        response90(mi) = resp / length(cellsTotal);
        
    end
end