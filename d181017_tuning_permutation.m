clear
mice = [25,27,30,36,37,39,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3],[6],[4],[4],[4]};  

numResampling = 10000;

baseD = 'Y:\Whiskernas\JK\suite2p\';
for mi = 1 : length(mice)
    cd(sprintf('%s%03d',baseD,mice(mi)));
    for si = 1 : length(sessions{mi})
        fprintf('Processing JK%03dS%02d\n', mice(mi), sessions{mi}(si));
        fn = sprintf('JK%03dS%02dsingleCell_anova_calcium.mat',mice(mi), sessions{mi}(si));
        savefn = sprintf('JK%03dS%02dpermtest_calcium.mat',mice(mi), sessions{mi}(si));
        dat = load(fn);
        fail = zeros(length(dat.cellsTuned),1);
        anovaP = cell(length(dat.cellsTuned),1);
        permmaxmod = cell(length(dat.cellsTuned),1);
        parfor ci = 1 : length(dat.cellsTuned)
            fprintf('%d/%d\n', ci, length(dat.cellsTuned))
            dF = dat.dFtotal{dat.cellsTotal==dat.cellsTuned(ci)};
            
            groupNums = cellfun(@(x) size(x,1), dF);
            groups = [];
            for i = 1 : length(groupNums)
                groups = [groups; ones(groupNums(i),1)*i];
            end
            
            meanF = cell2mat(cellfun(@(x) nanmean(x(:, dat.baseFrameNum+1 : dat.baseFrameNum+dat.afterFrameNum),2), dF, 'uniformoutput', false));
            [~, ~, anovaStats] = anova1(meanF, groups, 'off');            
            maxmod = max(anovaStats.means) - min(anovaStats.means);
            anovaP{ci} = zeros(numResampling,1);
            permmaxmod{ci} = zeros(numResampling,1);
            for ri = 1 : numResampling
                tempG = groups(randperm(length(groups),length(groups)));
                [anovaP{ci}(ri), ~, anovaStats] = anova1(meanF, tempG, 'off');
                permmaxmod{ci}(ri) = max(anovaStats.means) - min(anovaStats.means);
            end
            
            if length(find(anovaP{ci} < 0.05)) > 0.05 * numResampling
                if length(find(permmaxmod{ci} > maxmod)) > 0.05 * numResampling
                    fail(ci) = 1;
                end
            end            
        end
        save(savefn, 'anovaP', 'permmaxmod', 'fail')
    end
    
    for si = 1 : length(sessions{mi})
        fprintf('Processing JK%03dS%02d\n', mice(mi), sessions{mi}(si));
        fn = sprintf('JK%03dS%02dsingleCell_anova_spk.mat',mice(mi), sessions{mi}(si));
        savefn = sprintf('JK%03dS%02dpermtest_spk.mat',mice(mi), sessions{mi}(si));
        dat = load(fn);
        fail = zeros(length(dat.cellsTuned),1);
        anovaP = cell(length(dat.cellsTuned),1);
        permmaxmod = cell(length(dat.cellsTuned),1);
        parfor ci = 1 : length(dat.cellsTuned)
            fprintf('%d/%d\n', ci, length(dat.cellsTuned))
            infspk = dat.spkTotal{dat.cellsTotal==dat.cellsTuned(ci)};
            
            groupNums = cellfun(@(x) size(x,1), infspk);
            groups = [];
            for i = 1 : length(groupNums)
                groups = [groups; ones(groupNums(i),1)*i];
            end
            
            meanF = cell2mat(cellfun(@(x) nanmean(x(:, dat.baseFrameNum+1 : dat.baseFrameNum+dat.afterFrameNum),2), infspk, 'uniformoutput', false));
            [~, ~, anovaStats] = anova1(meanF, groups, 'off');            
            maxmod = max(anovaStats.means) - min(anovaStats.means);
            anovaP{ci} = zeros(numResampling,1);
            permmaxmod{ci} = zeros(numResampling,1);
            for ri = 1 : numResampling
                tempG = groups(randperm(length(groups),length(groups)));
                [anovaP{ci}(ri), ~, anovaStats] = anova1(meanF, tempG, 'off');
                permmaxmod{ci}(ri) = max(anovaStats.means) - min(anovaStats.means);
            end
            
            if length(find(anovaP{ci} < 0.05)) > 0.05 * numResampling
                if length(find(permmaxmod{ci} > maxmod)) > 0.05 * numResampling
                    fail(ci) = 1;
                end
            end            
        end
        save(savefn, 'anovaP', 'permmaxmod', 'fail')
    end
end

% %% test again with just anovaP
% clear
% mice = [25,27,30,36,37,39,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};
% baseD = 'C:\Data\suite2p\';
% pctAnova = [];
% pctMod = [];
% pctAnovaNMod = [];
% pctfail = [];
% numtuned = [];
% for mi = 1 : length(mice)
%     cd(sprintf('%s%03d',baseD,mice(mi)));
%     for si = 1 : length(sessions{mi})
% %     for si = 1  
%         load(sprintf('JK%03dS%02dpermtest.mat',mice(mi), sessions{mi}(si)));
%         fn = sprintf('JK%03dS%02dsingleCell_anova_FT_1secB.mat',mice(mi), sessions{mi}(si));
%         dat = load(fn);
%         
%         indAnova = find(cellfun(@(x) length(find(x < 0.05))>500, anovaP));
%         pctAnova = [pctAnova, length(indAnova) / length(fail)];
%         pctMod = [pctMod, length(find(cellfun(@(x) prctile(x,95), permmaxmod) > dat.tuneModulationMaxmin)) / length(fail)];
%         pctAnovaNMod = [pctAnovaNMod, length(find(cellfun(@(x) prctile(x,95), permmaxmod(indAnova)) > dat.tuneModulationMaxmin(indAnova))) / length(fail)];
%         numtuned = [numtuned, length(fail)];
%         pctfail = [pctfail, sum(fail)/length(fail)];
%     end
% end
% 
% %%
% pctAnova ./ numtuned