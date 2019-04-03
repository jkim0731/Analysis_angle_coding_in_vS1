mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};

behaviorDir = 'Y:\Whiskernas\JK\SoloData\';
whiskerDir = 'Y:\Whiskernas\JK\whisker\tracked\';

duration = 8; % in sec
inflation = 100; % for graphical clearance

trialNums = cell(length(mice),1);
performance = cell(length(mice),1);
correctLicks = cell(length(mice),1);
wrongLicks = cell(length(mice),1);
theta = cell(length(mice),1);
time = cell(length(mice),1);
%%
for mi = 1 : length(mice)
% for mi = 3
    mouseName = sprintf('JK%03d',mice(mi));
    cd([behaviorDir, mouseName])
    load(['behavior_', mouseName])
    trialNums{mi} = cell(length(sessions{mi}),1);
    performance{mi} = cell(length(sessions{mi}),1);
    correctLicks{mi} = cell(length(sessions{mi}),1);
    wrongLicks{mi} = cell(length(sessions{mi}),1);
    theta{mi} = cell(length(sessions{mi}),1);
    time{mi} = cell(length(sessions{mi}),1);
%     
    for si = 1 : length(sessions{mi})
%     for si = 8
        sessionName = sprintf('S%02d',sessions{mi}(si));
        bind = find(cellfun(@(x) strcmp(x.sessionName, sessionName), b));
        currB = b{bind};
%         lastTn = find(currB.trialCorrects ~= -1, 1, 'last');
        catchInd = find(cellfun(@(x) strcmp(x.trialType,'oo'), currB.trials));
        missInd = find(currB.trialCorrects == -1);
        behaviorTns = currB.trialNums(setdiff(1:length(currB.trials),union(missInd, catchInd)));
        
        cd([whiskerDir,mouseName,sessionName])
        wla = Whisker.WhiskerTrialLite_2padArray([whiskerDir,mouseName,sessionName]);
        whiskerTns = wla.trialNums;
        trialNums{mi}{si} = intersect(behaviorTns, whiskerTns);

        btninds = find(ismember(currB.trialNums, trialNums{mi}{si}));
        wtninds = find(ismember(wla.trialNums, trialNums{mi}{si}));
        performance{mi}{si} = currB.trialCorrects(btninds);
        
        correctLicks{mi}{si} = cell(length(trialNums{mi}{si}),1);
        wrongLicks{mi}{si} = cell(length(trialNums{mi}{si}),1);
        theta{mi}{si} = cell(length(trialNums{mi}{si}),1);
        time{mi}{si} = cell(length(trialNums{mi}{si}),1);
        
        for ti = 1 : length(trialNums{mi}{si})
            if currB.trials{btninds(ti)}.trialType(1) == 'r'
                wrongLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesLeft;
                correctLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesRight;
            elseif currB.trials{btninds(ti)}.trialType(1) == 'l'
                wrongLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesRight;
                correctLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesLeft;
            else
                error('trialType is wrong')
            end
            
            theta{mi}{si}{ti} = wla.trials{wtninds(ti)}.theta{1};
            time{mi}{si}{ti} = wla.trials{wtninds(ti)}.time{1};
        end
    end
% end
% %% Licking plot
% 
%     totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}));
%     bg = zeros(max(totalNumTrialsBySessions),duration * inflation); % 8 sec limit
%     totalNumTrialsBySessions = totalNumTrialsBySessions(1:end-1);
%     totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
%         
%     figure('units', 'normalized', 'position', [mi/20, 0, 0.4, 1]), 
%     imshow(bg), axis off, hold on,
%     for si = 1 : length(totalNumTrialsBySessions)
%         for j = 1 : length(trialNums{mi}{si})
%             click = correctLicks{mi}{si}{j} * inflation;
%             click = click(click < duration * inflation);
%             click = unique(click);
%             wlick = wrongLicks{mi}{si}{j} * inflation;
%             wlick = wlick(wlick < duration * inflation);
%             wlick = unique(wlick);
%             plot(click, ones(length(click),1) * totalNumTrialsBySessions(si)+j, 'g.', 'markersize',0.1)
%             plot(wlick, ones(length(wlick),1) * totalNumTrialsBySessions(si)+j, 'r.', 'markersize',0.1)
%         end
%     end
%     title([mouseName, ' Licking'])
% %% Whisking plot
% %     maxTheta = prctile(cell2mat(cellfun(@(x) x, cellfun(@(x) x,theta{mi}{si}, 'uniformoutput',false), 'uniformoutput', false)),90);
% %     minTheta = prctile(cell2mat(cellfun(@(x) x, cellfun(@(x) x,theta{mi}{si}, 'uniformoutput',false), 'uniformoutput', false)),10);
%     times = 0:1/inflation:duration;
%     
%     w1 = zeros(max(totalNumTrialsBySessions),duration * inflation);    
%     
%     for si = 1 : length(totalNumTrialsBySessions)
%         for j = 1 : length(trialNums{mi}{si})
%             [~, amplitude, filteredSignal] =  jkWhiskerDecomposition(theta{mi}{si}{j});
%             for k = 2 : length(times)
%                 tempInds = find(time{mi}{si}{j} > times(k-1) & time{mi}{si}{j} < times(k));
%                 w1(totalNumTrialsBySessions(si)+j, k-1) = mean(amplitude(tempInds));
%             end
%         end
%     end
% 
% %     figure('units', 'normalized', 'position', [mi/20, 0, 0.6, 1]),
% %     imagesc(w1, [-10 10]), axis off
%     figure('units', 'normalized', 'position', [mi/20, 0, 0.6, 1]),
%     imagesc(w1, [-10 10]), axis off   
%     hold on, plot([99 99], [0, size(bg,1)], 'k--', 'linewidth', 2), plot([199 199], [0, size(bg,1)], 'k--', 'linewidth', 2)
%     
%     title([mouseName, 'Whisking'])

end
%%
cd(whiskerDir)
save('WholeSessionBehaviorNWhisking.mat', 'trialNums', 'performance', 'correctLicks', 'wrongLicks', 'theta', 'time', 'duration', 'inflation', 'mice', 'sessions')

% %%
% duration = 8; % in sec
% inflation = 100; % for graphical clearance
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% %% Licking plot
% 
% mi = 1;
% mouseName = sprintf('JK%03d',mice(mi));
% totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}(1:end-0)));
% 
% bg = zeros(max(totalNumTrialsBySessions),(duration) * inflation); % 8 sec limit
% totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
% 
% c = gray(100); 
% 
% figure('units', 'normalized', 'position', [mi/20, 0, 0.4, 1]), 
% imshow(bg), axis off, axis image, hold on,
% for si = 1 : length(totalNumTrialsBySessions)-1
%     for j = 1 : length(trialNums{mi}{si})
%         click = correctLicks{mi}{si}{j} * inflation;
%         click = click(click < duration * inflation);
%         click = unique(click);
%         wlick = wrongLicks{mi}{si}{j} * inflation;
%         wlick = wlick(wlick < duration * inflation);
%         wlick = unique(wlick);
%         plot(click, ones(length(click),1) * totalNumTrialsBySessions(si)+j, '.', 'markersize',0.1,'color',[86,180,233]/255)
%         plot(wlick, ones(length(wlick),1) * totalNumTrialsBySessions(si)+j, '.', 'markersize',0.1,'color',[230,159,0]/255)
%     end
%     cr = length(find(performance{mi}{si}))/length(find(performance{mi}{si}>-1));
% %     cind = min(max(round((round(cr,2) - 0.4) * 100),1),50) + 20;
%     cind = round(cr*100);
%     plot([-inflation, -inflation], [totalNumTrialsBySessions(si), totalNumTrialsBySessions(si+1)], 'color', c(cind,:), 'linewidth', 10)
%     plot([-inflation *1.5, -inflation *0.5], [totalNumTrialsBySessions(si), totalNumTrialsBySessions(si)], 'k-', 'linewidth', 2)
% end
% plot([-inflation *1.5, -inflation *0.5], [totalNumTrialsBySessions(end), totalNumTrialsBySessions(end)], 'k-', 'linewidth', 2)
% patch([-inflation*1.5, -inflation*0.5, -inflation*0.5,-inflation*1.5], [1, 1, size(bg,1), size(bg,1)], 'k-', 'facealpha',0, 'linewidth', 2)
% xlim([-inflation *1.5, size(bg,2)])
% % title([mouseName, ' Licking'])
% 
% %%
% mi = 10;
% mouseName = sprintf('JK%03d',mice(mi));
% totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}(1:end-0)));
% 
% bg = zeros(max(totalNumTrialsBySessions),(duration) * inflation); % 8 sec limit
% totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
% 
% c = gray(100); 
% 
% figure('units', 'normalized', 'position', [mi/20, 0, 0.4, 1]), 
% imshow(bg), axis off, axis image, hold on,
% for si = 1 : length(totalNumTrialsBySessions)-1
%     for j = 1 : length(trialNums{mi}{si})
%         click = correctLicks{mi}{si}{j} * inflation;
%         click = click(click < duration * inflation);
%         click = unique(click);
%         wlick = wrongLicks{mi}{si}{j} * inflation;
%         wlick = wlick(wlick < duration * inflation);
%         wlick = unique(wlick);
%         plot(click, ones(length(click),1) * totalNumTrialsBySessions(si)+j, '.', 'markersize',0.1,'color',[86,180,233]/255)
%         plot(wlick, ones(length(wlick),1) * totalNumTrialsBySessions(si)+j, '.', 'markersize',0.1,'color',[230,159,0]/255)
%     end
%     cr = length(find(performance{mi}{si}))/length(find(performance{mi}{si}>-1));
% %     cind = min(max(round((round(cr,2) - 0.4) * 100),1),50) + 20;
%     cind = round(cr*100);
%     plot([-inflation, -inflation], [totalNumTrialsBySessions(si), totalNumTrialsBySessions(si+1)], 'color', c(cind,:), 'linewidth', 7)
%     plot([-inflation *1.5, -inflation *0.5], [totalNumTrialsBySessions(si), totalNumTrialsBySessions(si)], 'k-', 'linewidth', 2)
% end
% plot([-inflation *1.5, -inflation *0.5], [totalNumTrialsBySessions(end), totalNumTrialsBySessions(end)], 'k-', 'linewidth', 2)
% patch([-inflation*1.5, -inflation*0.5, -inflation*0.5,-inflation*1.5], [1, 1, size(bg,1), size(bg,1)], 'k-', 'facealpha',0, 'linewidth', 2)
% xlim([-inflation *1.5, size(bg,2)])
% % title([mouseName, ' Licking'])
% 
% %%
% figure, hold on
% c = gray(100); 
% for i = 1 : 10 : size(c,1)
%     plot([i i+1], [i i], 'color', c(i,:))
% end
% colormap(gray), colorbar('westoutside', 'ticks', []), set(gca, 'linewidth', 3, 'fontweight', 'bold', 'fontsize', 12)
% %% Whisking plot
% mi = 8;
% inflation = 100;
% duration = 4;
% times = 0:1/inflation:duration;
% mirrorAngleOffset = 13; % in degrees
% totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}(1:end)));
% totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
% w1 = zeros(max(totalNumTrialsBySessions),duration * inflation);
% 
% for si = 1 : length(totalNumTrialsBySessions)-1
%     for j = 1 : length(trialNums{mi}{si})
%         [~, amplitude, filteredSignal] =  jkWhiskerDecomposition(theta{mi}{si}{j});
%         for k = 2 : length(times)
%             tempInds = find(time{mi}{si}{j} > times(k-1) & time{mi}{si}{j} < times(k));
%             w1(totalNumTrialsBySessions(si)+j, k-1) = mean(amplitude(tempInds));
% %             w1(totalNumTrialsBySessions(si)+j, k-1) = nanmean(theta{mi}{si}{j}(tempInds))+mirrorAngleOffset;
%         end
%     end
% end
% %
% % nanind = zeros(size(w1,1),1);
% % for wi = 1 : size(w1,1)
% %     nanind(wi) = prod(isnan(w1(wi,:)));
% % end
% % w1(nanind,:) = [];
% %
% figure('units', 'normalized', 'position', [mi/20, 0, 0.2, 1]),
% % figure,
% imagesc(w1, [0 10]), axis off
% hold on, plot([99 99], [0, size(w1,1)], 'k--', 'linewidth', 2), plot([199 199], [0, size(w1,1)], 'k--', 'linewidth', 2)
% 
% % title([mouseName, 'Whisking'])
% % title('Non-learner')
% %% Graph for licking and whisking, averaged in each session in each mouse
% % first, look at their time series. Then, compare stereotypy by correlation
% % analysis.
% % for licking, it will be smoothed lick rate, on right/wrong and left/right
% % % (the latter to check their bias)
% close all
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% step = 0.1; % lick histogram step, in s
% maxTime = 8; % in s
% bins = 0:step:maxTime;
% for mi = 1 : length(learner)
%     mouse = mice(learner(mi));
%     naiveWrongLick = zeros(length(bins)-1,1);
%     naiveCorrectLick = zeros(length(bins)-1,1);
%     expertWrongLick = zeros(length(bins)-1,1);
%     expertCorrectLick = zeros(length(bins)-1,1);
%     
%     for si = 1 : 3
%         temp = wrongLicks{learner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveWrongLick = naiveWrongLick + tempwl/3;
%         temp = correctLicks{learner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveCorrectLick = naiveCorrectLick + tempcl/3;
%     end
%     
%     for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
%         temp = wrongLicks{learner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertWrongLick = expertWrongLick + tempwl/3;
%         temp = correctLicks{learner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertCorrectLick = expertCorrectLick + tempcl/3;
%     end
%     
%     figure,
%     subplot(211), plot(bins(1:end-1)+step/2, naiveCorrectLick), hold on, plot(bins(1:end-1)+step/2, naiveWrongLick), title(sprintf('JK%03d', mouse));
%     plot([1 1], [0, max(max(naiveCorrectLick), max(naiveWrongLick))], 'k--')
%     plot([2 2], [0, max(max(naiveCorrectLick), max(naiveWrongLick))], 'k--')
%     subplot(212), plot(bins(1:end-1)+step/2, expertCorrectLick), hold on, plot(bins(1:end-1)+step/2, expertWrongLick), xlabel('Time (s)')
%     plot([1 1], [0, max(max(expertCorrectLick), max(expertWrongLick))], 'k--')
%     plot([2 2], [0, max(max(expertCorrectLick), max(expertWrongLick))], 'k--')
% end
% 
% %%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% step = 0.1; % lick histogram step, in s
% maxTime = 8; % in s
% bins = 0:step:maxTime;
% for mi = 1 : length(nonlearner)
%     mouse = mice(nonlearner(mi));
%     naiveWrongLick = zeros(length(bins)-1,1);
%     naiveCorrectLick = zeros(length(bins)-1,1);
%     expertWrongLick = zeros(length(bins)-1,1);
%     expertCorrectLick = zeros(length(bins)-1,1);
%     
%     for si = 1 : 3
%         temp = wrongLicks{nonlearner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveWrongLick = naiveWrongLick + tempwl/3;
%         temp = correctLicks{nonlearner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveCorrectLick = naiveCorrectLick + tempcl/3;
%     end
%     
%     for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
%         temp = wrongLicks{nonlearner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertWrongLick = expertWrongLick + tempwl/3;
%         temp = correctLicks{nonlearner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertCorrectLick = expertCorrectLick + tempcl/3;
%     end
%     
%     figure,
%     subplot(211), plot(bins(1:end-1)+step/2, naiveCorrectLick), hold on, plot(bins(1:end-1)+step/2, naiveWrongLick), title(sprintf('JK%03d', mouse));
%     plot([1 1], [0, max(max(naiveCorrectLick), max(naiveWrongLick))], 'k--')
%     plot([2 2], [0, max(max(naiveCorrectLick), max(naiveWrongLick))], 'k--')
%     subplot(212), plot(bins(1:end-1)+step/2, expertCorrectLick), hold on, plot(bins(1:end-1)+step/2, expertWrongLick), xlabel('Time (s)')
%     plot([1 1], [0, max(max(expertCorrectLick), max(expertWrongLick))], 'k--')
%     plot([2 2], [0, max(max(expertCorrectLick), max(expertWrongLick))], 'k--')
% end
%     
% %%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% step = 0.1; % lick histogram step, in s
% maxTime = 8; % in s
% bins = 0:step:maxTime;
% naiveWrongLick = zeros(length(bins)-1,length(learner));
% naiveCorrectLick = zeros(length(bins)-1,length(learner));
% expertWrongLick = zeros(length(bins)-1,length(learner));
% expertCorrectLick = zeros(length(bins)-1,length(learner));
% firstWrongLick = zeros(length(bins)-1,length(nonlearner));
% firstCorrectLick = zeros(length(bins)-1,length(nonlearner));
% lastWrongLick = zeros(length(bins)-1,length(nonlearner));
% lastCorrectLick = zeros(length(bins)-1,length(nonlearner));
% for mi = 1 : length(learner)
%     for si = 1 : 3
%         temp = wrongLicks{learner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveWrongLick(:,mi) = naiveWrongLick(:,mi) + tempwl/3;
%         temp = correctLicks{learner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         naiveCorrectLick(:,mi) = naiveCorrectLick(:,mi) + tempcl/3;
%     end
%     
%     for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
%         temp = wrongLicks{learner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertWrongLick(:,mi) = expertWrongLick(:,mi) + tempwl/3;
%         temp = correctLicks{learner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         expertCorrectLick(:,mi) = expertCorrectLick(:,mi) + tempcl/3;
%     end
% end
% 
% for mi = 1 : length(nonlearner)
%     for si = 1 : 3
%         temp = wrongLicks{nonlearner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         firstWrongLick(:,mi) = firstWrongLick(:,mi) + tempwl/3;
%         temp = correctLicks{nonlearner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         firstCorrectLick(:,mi) = firstCorrectLick(:,mi) + tempcl/3;
%     end
%     
%     for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
%         temp = wrongLicks{nonlearner(mi)}{si};
%         tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         lastWrongLick(:,mi) = lastWrongLick(:,mi) + tempwl/3;
%         temp = correctLicks{nonlearner(mi)}{si};
%         tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
%         lastCorrectLick(:,mi) = lastCorrectLick(:,mi) + tempcl/3;
%     end
% end
% % %%
% % figure,
% % subplot(211), errorbar(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), std(naiveCorrectLick,[],2), 'linewidth', 5, 'color', [86,180,223]/255), hold on,
% % errorbar(bins(1:end-1)+step/2, mean(naiveWrongLick,2), std(naiveWrongLick,[],2), 'linewidth', 5, 'color', [230,159,0]/255), 
% % 
% % errorbar(bins(1:end-1)+step/2, mean(expertCorrectLick,2), std(expertCorrectLick,[],2), 'linewidth', 5, 'color', [0,114,178]/255), hold on,
% % errorbar(bins(1:end-1)+step/2, mean(expertWrongLick,2), std(expertWrongLick,[],2), 'linewidth', 5, 'color', [213,94,0]/255), 
% % 
% % plot([1 1], [0 max(max([expertWrongLick;expertCorrectLick;naiveWrongLick;naiveCorrectLick]))], 'k--')
% % plot([2 2], [0 max(max([expertWrongLick;expertCorrectLick;naiveWrongLick;naiveCorrectLick]))], 'k--')
% % xlabel('Time (s)'), title('Learners')
% % 
% % %%
% % figure,
% % subplot(211), errorbar(bins(1:end-1)+step/2, mean(firstCorrectLick,2), std(firstCorrectLick,[],2)), hold on,
% % errorbar(bins(1:end-1)+step/2, mean(firstWrongLick,2), std(firstWrongLick,[],2)), title('Learners')
% % plot([1 1], [0 max(max(max(firstWrongLick)), max(max(firstCorrectLick)))], 'k--')
% % plot([2 2], [0 max(max(max(firstWrongLick)), max(max(firstCorrectLick)))], 'k--')
% % title('Non-learners')
% % 
% % subplot(212), errorbar(bins(1:end-1)+step/2, mean(lastCorrectLick,2), std(lastCorrectLick,[],2)), hold on,
% % errorbar(bins(1:end-1)+step/2, mean(lastWrongLick,2), std(lastWrongLick,[],2)), hold on,
% % plot([1 1], [0 max(max(max(lastWrongLick)), max(max(lastCorrectLick)))], 'k--')
% % plot([2 2], [0 max(max(max(lastWrongLick)), max(max(lastCorrectLick)))], 'k--')
% % xlabel('Time (s)')
% 
% %%
% c = parula(10);
% c1 = c(4,:);
% c2 = c(2,:);
% w1 = c(10,:);
% w2 = c(8,:);
% step = 0.1; % lick histogram step, in s
% maxTime = 8; % in s
% bins = 0:step:maxTime;
% figure, subplot(211), 
% plot(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
% hold on,
% plot(bins(1:end-1)+step/2, mean(expertCorrectLick,2), 'color', c1, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(naiveWrongLick,2), '--', 'color', w2, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(expertWrongLick,2), 'color', w2, 'linewidth', 3);
% 
% boundedline(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), std(naiveCorrectLick, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(bins(1:end-1)+step/2, mean(expertCorrectLick,2), std(expertCorrectLick, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(bins(1:end-1)+step/2, mean(naiveWrongLick,2), std(naiveWrongLick, [], 2), '--', 'cmap', w2, 'alpha', 'transparency', 0.2);
% boundedline(bins(1:end-1)+step/2, mean(expertWrongLick,2), std(expertWrongLick, [], 2), 'cmap', w2, 'alpha', 'transparency', 0.2);
% 
% plot(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(expertCorrectLick,2), 'color', c1, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(naiveWrongLick,2), '--', 'color', w2, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(expertWrongLick,2), 'color', w2, 'linewidth', 3);
% % plot([1 1], [0 max(max([expertWrongLick;expertCorrectLick;naiveWrongLick;naiveCorrectLick]))], 'k--')
% % plot([2 2], [0 max(max([expertWrongLick;expertCorrectLick;naiveWrongLick;naiveCorrectLick]))], 'k--')
% plot([1 1], [0 6], 'k--')
% plot([2 2], [0 6], 'k--')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% title('Learners'), ylabel('Lick rate (Hz)'), ylim([0 6]), xlim([0 6]), xticklabels({})
% legend({'Correct (first 3 sessions)'; 'Correct (last 3 sessions)'; 'Wrong (first 3 sessions)'; 'Wrong (last 3 sessions)'}, 'box', 'off')
% 
% subplot(212), boundedline(bins(1:end-1)+step/2, mean(firstCorrectLick,2), std(firstCorrectLick, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
% hold on,
% % set(h1,'linewidth',5)
% boundedline(bins(1:end-1)+step/2, mean(lastCorrectLick,2), std(lastCorrectLick, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(bins(1:end-1)+step/2, mean(firstWrongLick,2), std(firstWrongLick, [], 2), '--', 'cmap', w2, 'alpha', 'transparency', 0.2);
% boundedline(bins(1:end-1)+step/2, mean(lastWrongLick,2), std(lastWrongLick, [], 2), 'cmap', w2, 'alpha', 'transparency', 0.2);
% 
% plot(bins(1:end-1)+step/2, mean(firstCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(lastCorrectLick,2), 'color', c1, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(firstWrongLick,2), '--','color', w2, 'linewidth', 3);
% plot(bins(1:end-1)+step/2, mean(lastWrongLick,2), 'color', w2, 'linewidth', 3);
% % plot([1 1], [0 max(max([lastWrongLick;lastCorrectLick;firstWrongLick;firstCorrectLick]))], 'k--')
% % plot([2 2], [0 max(max([lastWrongLick;lastCorrectLick;firstWrongLick;firstCorrectLick]))], 'k--')
% plot([1 1], [0 6], 'k--')
% plot([2 2], [0 6], 'k--')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')
% title('Non-learners'), ylabel('Lick rate (Hz)'), ylim([0 6]), xlim([0 6]), xlabel('Time (s)')
% 
% %%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% step = 0.1; % lick histogram step, in s
% maxTime = 8; % in s
% bins = 0:step:maxTime;
% learnerlicktemplate = zeros(length(bins)-1,length(learner));
% nonlearnerlicktemplate = zeros(length(bins)-1,length(nonlearner));
% for mi = 1 : length(learner)
%     for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
%         tempwl = wrongLicks{learner(mi)}{si};
%         tempcl = correctLicks{learner(mi)}{si};
%         temphc = smooth(histcounts(union(cell2mat(tempwl), cell2mat(tempcl)), bins)) / step / length(tempwl);        
%         learnerlicktemplate(:,mi) = learnerlicktemplate(:,mi) + temphc/3;
%     end
% end
% 
% for mi = 1 : length(nonlearner)
%     for si = length(sessions{nonlearner(mi)}) - 4 : length(sessions{nonlearner(mi)}) - 2
%         tempwl = wrongLicks{nonlearner(mi)}{si};
%         tempcl = correctLicks{nonlearner(mi)}{si};
%         temphc = smooth(histcounts(union(cell2mat(tempwl), cell2mat(tempcl)), bins)) / step / length(tempwl);        
%         nonlearnerlicktemplate(:,mi) = nonlearnerlicktemplate(:,mi) + temphc/3;
%     end
% end
% %%
% learnerLickCorr = cell(length(learner),1);
% for mi = 1 : length(learner)
%     learnerLickCorr{mi} = zeros(length(sessions{learner(mi)})-2,1);
%     for si = 1 : length(sessions{learner(mi)})-2
%         tempwl = wrongLicks{learner(mi)}{si};
%         tempcl = correctLicks{learner(mi)}{si};
%         temphc = smooth(histcounts(union(cell2mat(tempwl), cell2mat(tempcl)), bins)) / step / length(tempwl);
%         learnerLickCorr{mi}(si) = corr(temphc, learnerlicktemplate(:,mi));
%     end
% end
% %%
% nonlearnerLickCorr = cell(length(nonlearner),1);
% for mi = 1 : length(nonlearner)
% % for mi = 3
%     nonlearnerLickCorr{mi} = zeros(length(sessions{nonlearner(mi)})-2,1);
%     for si = 1 : length(sessions{nonlearner(mi)})
% %     for si = 1
%         tempwl = wrongLicks{nonlearner(mi)}{si};
%         tempcl = correctLicks{nonlearner(mi)}{si};
%         temphc = smooth(histcounts(union(cell2mat(tempwl), cell2mat(tempcl)), bins)) / step / length(tempwl);
%         nonlearnerLickCorr{mi}(si) = corr(temphc, nonlearnerlicktemplate(:,mi));
%     end
% end
% 
% %%
% mi = 3;
% si = 8;
% tempwl = wrongLicks{learner(mi)}{si};
% tempcl = correctLicks{learner(mi)}{si};
% temphc = smooth(histcounts(union(cell2mat(tempwl), cell2mat(tempcl)), bins)) / step / length(tempwl);
% figure, plot(bins(1:end-1)+step/2, temphc, 'r-', 'linewidth', 3), hold on
% plot(bins(1:end-1)+step/2, learnerlicktemplate(:,mi), 'k-', 'linewidth', 3)
% ylabel('Lick rate (Hz)'), ylim([0 7]), xlim([0 6]), xlabel('Time (s)')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% title(sprintf('JK%03d S%02d VS template',mice(learner(mi)), sessions{learner(mi)}(si)))
% 
% %%
% figure,
% hold on
% % subplot(211), hold on
% for i = 1 : length(learner)
% %     plot(-length(learnerLickCorr{i})+1+3:0, learnerLickCorr{i}(1:end-3), 'k-', 'linewidth', 3)
%     plot(1:length(learnerLickCorr{i}), learnerLickCorr{i}(1:end), 'k-', 'linewidth', 3)
% end
% % plot([-max(cellfun(@(x) length(x), learnerLickCorr))+1+3, 0], [0.8, 0.8], 'k--', 'linewidth', 1)
% % title('Learners'), ylabel('Correlation'), xlabel('Session before learning')
% % set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% % subplot(212), hold on
% for i = 1 : length(nonlearner)
%     plot(1:length(nonlearnerLickCorr{i}), nonlearnerLickCorr{i}(1:end), 'color', [0.7 0.7 0.7], 'linewidth', 3)
% end
% % plot([-max(cellfun(@(x) length(x), nonlearnerLickCorr))+1, 0], [0.8, 0.8], 'k--', 'linewidth', 1)
% 
% title('Lick pattern stereotypy'), ylabel('Correlation'), xlabel('Session')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% 
% %% Comparison between templates
% 
% corrmat = zeros(length(mi), length(mi));
% for i = 1 : length(learner)
%     for j = i : length(learner)
%         corrmat(i,j) = corr(learnerlicktemplate(:,i), learnerlicktemplate(:,j));
%     end
%     for j = 1 : length(nonlearner)
%         corrmat(i, j+length(learner)) = corr(learnerlicktemplate(:,i), nonlearnerlicktemplate(:,j));
%     end
% end
% for i = 1 : length(nonlearner)
%     for j = i : length(nonlearner)
%         corrmat(i+length(learner), j + length(learner)) = corr(nonlearnerlicktemplate(:,i), nonlearnerlicktemplate(:,j));
%     end
% end
% figure, imagesc(corrmat),
% title('Between templates')
% %%
% for i = 1 : length(mi)
%     corrmat(i,i) = 0;
% end
% co = corrmat(:);
% mean(co(find(co>0.1)))
% 
% %%
% %%
% %%  Whisking
% %%
% %%
% %%
% 
% mi = 10;
% times = 0:1/inflation:duration;
% totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}(1:end-0)));
% 
% totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
% 
% w1 = nan(max(totalNumTrialsBySessions),duration * inflation);    
% 
% for si = 1 : length(totalNumTrialsBySessions)-1
%     if ~isempty(theta{mi}{si})
%         for j = 1 : length(trialNums{mi}{si})
%             [~, amplitude, filteredSignal] =  jkWhiskerDecomposition(theta{mi}{si}{j});
%             for k = 2 : length(times)
%                 tempInds = find(time{mi}{si}{j} > times(k-1) & time{mi}{si}{j} < times(k));
%                 w1(totalNumTrialsBySessions(si)+j, k-1) = mean(amplitude(tempInds));
%             end
%         end
%     end
% end
% %
% figure('units', 'normalized', 'position', [mi/20, 0, 0.6, 1]),
% imagesc(w1, [-10 10]), axis off   
% hold on, plot([99 99], [0, size(w1,1)], 'k--', 'linewidth', 2), plot([199 199], [0, size(w1,1)], 'k--', 'linewidth', 2)
% 
% 
% %% before and after training (first 3 and last 3)
% behaviorDir = 'Y:\Whiskernas\JK\SoloData\';
% whiskerDir = 'Y:\Whiskernas\JK\whisker\tracked\';
% videoFreq = 311.24;
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% learner = [1,2,3,4,7,9];
% nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
% maxLengthTime = 4; % in sec
% maxLength = ceil(maxLengthTime * videoFreq);
% naiveTheta = zeros(maxLength,length(learner));
% expertTheta = zeros(maxLength,length(learner));
% firstTheta = zeros(maxLength,length(nonlearner));
% lastTheta = zeros(maxLength,length(nonlearner));
% naiveMidpoint = zeros(maxLength,length(nonlearner));
% expertMidpoint = zeros(maxLength,length(nonlearner));
% firstMidpoint = zeros(maxLength,length(nonlearner));
% lastMidpoint = zeros(maxLength,length(nonlearner));
% naiveAmplitude = zeros(maxLength,length(nonlearner));
% expertAmplitude = zeros(maxLength,length(nonlearner));
% firstAmplitude = zeros(maxLength,length(nonlearner));
% lastAmplitude = zeros(maxLength,length(nonlearner));
% 
% %%
% for mi = 1 : length(learner)
%     for si = 1 : 3
%         temp = theta{learner(mi)}{si};
%         tempTheta = nan(maxLength, length(temp));
%         tempMidpoint = nan(maxLength, length(temp));
%         tempAmplitude = nan(maxLength, length(temp));
%         for i = 1 : length(temp)
%             indLength = min(maxLength, length(temp{i}));
%             tempTheta(1:indLength,i) = temp{i}(1:indLength);
%             [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
%         end
%         naiveTheta(:,mi) = naiveTheta(:,mi) + nanmean(tempTheta,2)/3;
%         naiveMidpoint(:,mi) = naiveMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
%         naiveAmplitude(:,mi) = naiveAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
%     end
%     
%     for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
%         temp = theta{learner(mi)}{si};
%         tempTheta = nan(maxLength, length(temp));
%         for i = 1 : length(temp)
%             indLength = min(maxLength, length(temp{i}));
%             tempTheta(1:indLength, i) = temp{i}(1:indLength);
%             [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
%         end
%         expertTheta(:,mi) = expertTheta(:,mi) + nanmean(tempTheta,2)/3;
%         expertMidpoint(:,mi) = expertMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
%         expertAmplitude(:,mi) = expertAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
%     end
% end
% 
% for mi = 1 : length(nonlearner)
%     for si = 1 : 3
%         temp = theta{nonlearner(mi)}{si};
%         tempTheta = nan(maxLength, length(temp));
%         for i = 1 : length(temp)
%             indLength = min(maxLength, length(temp{i}));
%             tempTheta(1:indLength, i) = temp{i}(1:indLength);
%             [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
%         end
%         firstTheta(:,mi) = firstTheta(:,mi) + nanmean(tempTheta,2)/3;
%         firstMidpoint(:,mi) = firstMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
%         firstAmplitude(:,mi) = firstAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
%     end
%     
%     for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
%         temp = theta{nonlearner(mi)}{si};
%         tempTheta = nan(maxLength, length(temp));
%         for i = 1 : length(temp)
%             indLength = min(maxLength, length(temp{i}));
%             tempTheta(1:indLength, i) = temp{i}(1:indLength);
%             [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
%         end
%         lastTheta(:,mi) = lastTheta(:,mi) + nanmean(tempTheta,2)/3;
%         lastMidpoint(:,mi) = lastMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
%         lastAmplitude(:,mi) = lastAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
%     end
% end
% %%
% c = parula(10);
% c1 = c(4,:);
% c2 = c(8,:);
% times = 0:1/videoFreq:1/videoFreq*(maxLength-1);
% mirrorAngleOffset = 13; % in degrees
% figure, 
% plot(times, mean(naiveTheta,2)+mirrorAngleOffset, '--', 'color', c1, 'linewidth', 3);
% hold on,
% plot(times, mean(expertTheta,2)+mirrorAngleOffset, 'color', c1, 'linewidth', 3);
% plot(times, mean(firstTheta,2)+mirrorAngleOffset, '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastTheta,2)+mirrorAngleOffset, 'color', c2, 'linewidth', 3);
% 
% boundedline(times, mean(naiveTheta,2)+mirrorAngleOffset, std(naiveTheta, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(expertTheta,2)+mirrorAngleOffset, std(expertTheta, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(firstTheta,2)+mirrorAngleOffset, std(firstTheta, [], 2), '--', 'cmap', c2, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(lastTheta,2)+mirrorAngleOffset, std(lastTheta, [], 2), 'cmap', c2, 'alpha', 'transparency', 0.2);
% 
% plot(times, mean(naiveTheta+mirrorAngleOffset,2), '--', 'color', c1, 'linewidth', 3);
% plot(times, mean(expertTheta+mirrorAngleOffset,2), 'color', c1, 'linewidth', 3);
% plot(times, mean(firstTheta+mirrorAngleOffset,2), '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastTheta+mirrorAngleOffset,2), 'color', c2, 'linewidth', 3);
% plot([1 1], [-25 5], 'k--')
% plot([2 2], [-25 5], 'k--')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% title('Average whisker angle'), ylabel('Theta (\circ)'), ylim([-25 5]), xlim([0 4]), xlabel('Time (s)')
% legend({'Learned (first 3 sessions)'; 'Learned (last 3 sessions)'; 'Not-learned (first 3 sessions)'; 'Not-learned (last 3 sessions)'}, 'box', 'off')
% 
% 
% %%
% 
% mirrorAngleOffset = 13; % in degrees
% figure, 
% plot(times, mean(naiveMidpoint,2)+mirrorAngleOffset, '--', 'color', c1, 'linewidth', 3);
% hold on,
% plot(times, mean(expertMidpoint,2)+mirrorAngleOffset, 'color', c1, 'linewidth', 3);
% plot(times, mean(firstMidpoint,2)+mirrorAngleOffset, '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastMidpoint,2)+mirrorAngleOffset, 'color', c2, 'linewidth', 3);
% 
% boundedline(times, mean(naiveMidpoint,2)+mirrorAngleOffset, std(naiveMidpoint, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(expertMidpoint,2)+mirrorAngleOffset, std(expertMidpoint, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(firstMidpoint,2)+mirrorAngleOffset, std(firstMidpoint, [], 2), '--', 'cmap', c2, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(lastMidpoint,2)+mirrorAngleOffset, std(lastMidpoint, [], 2), 'cmap', c2, 'alpha', 'transparency', 0.2);
% 
% plot(times, mean(naiveMidpoint+mirrorAngleOffset,2), '--', 'color', c1, 'linewidth', 3);
% plot(times, mean(expertMidpoint+mirrorAngleOffset,2), 'color', c1, 'linewidth', 3);
% plot(times, mean(firstMidpoint+mirrorAngleOffset,2), '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastMidpoint+mirrorAngleOffset,2), 'color', c2, 'linewidth', 3);
% plot([1 1], [-25 5], 'k--')
% plot([2 2], [-25 5], 'k--')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% title('Average whisker midpoint'), ylabel('Midpoint (\circ)'), ylim([-25 5]), xlim([0 4]), xlabel('Time (s)')
% legend({'Learned (first 3 sessions)'; 'Learned (last 3 sessions)'; 'Not-learned (first 3 sessions)'; 'Not-learned (last 3 sessions)'}, 'box', 'off')
% 
% 
% %%
% figure, 
% plot(times, mean(naiveAmplitude,2), '--', 'color', c1, 'linewidth', 3);
% hold on,
% plot(times, mean(expertAmplitude,2), 'color', c1, 'linewidth', 3);
% plot(times, mean(firstAmplitude,2), '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastAmplitude,2), 'color', c2, 'linewidth', 3);
% 
% boundedline(times, mean(naiveAmplitude,2), std(naiveAmplitude, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(expertAmplitude,2), std(expertAmplitude, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(firstAmplitude,2), std(firstAmplitude, [], 2), '--', 'cmap', c2, 'alpha', 'transparency', 0.2);
% boundedline(times, mean(lastAmplitude,2), std(lastAmplitude, [], 2), 'cmap', c2, 'alpha', 'transparency', 0.2);
% 
% plot(times, mean(naiveAmplitude,2), '--', 'color', c1, 'linewidth', 3);
% plot(times, mean(expertAmplitude,2), 'color', c1, 'linewidth', 3);
% plot(times, mean(firstAmplitude,2), '--', 'color', c2, 'linewidth', 3);
% plot(times, mean(lastAmplitude,2), 'color', c2, 'linewidth', 3);
% plot([1 1], [0 10], 'k--')
% plot([2 2], [0 10], 'k--')
% set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
% title('Average whisking amplitude'), ylabel('Amplitude (\circ)'), ylim([0 10]), xlim([0 4]), xlabel('Time (s)')
% legend({'Learned (first 3 sessions)'; 'Learned (last 3 sessions)'; 'Not-learned (first 3 sessions)'; 'Not-learned (last 3 sessions)'}, 'box', 'off')
% 


%%
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};

behaviorDir = 'Y:\Whiskernas\JK\SoloData\';
whiskerDir = 'Y:\Whiskernas\JK\whisker\tracked\';

touch = cell(length(mice),1);

for mi = 1 : length(mice)
    mouseName = sprintf('JK%03d',mice(mi));
    cd([behaviorDir, mouseName])
    load(['behavior_', mouseName])

    touch{mi} = cell(length(sessions{mi}),1);

    for si = 1 : length(sessions{mi})
%     for si = 8
        sessionName = sprintf('S%02d',sessions{mi}(si));
        bind = find(cellfun(@(x) strcmp(x.sessionName, sessionName), b));
        currB = b{bind};
%         lastTn = find(currB.trialCorrects ~= -1, 1, 'last');
        catchInd = find(cellfun(@(x) strcmp(x.trialType,'oo'), currB.trials));
        missInd = find(currB.trialCorrects == -1);
        behaviorTns = currB.trialNums(setdiff(1:length(currB.trials),union(missInd, catchInd)));
        
        cd([whiskerDir,mouseName,sessionName])
        wla = Whisker.WhiskerTrialLite_2padArray([whiskerDir,mouseName,sessionName]);
        wtninds = find(ismember(wla.trialNums, trialNums{mi}{si}));
        touch{mi}{si} = cell(length(trialNums{mi}{si}),1);
        for ti = 1 : length(trialNums{mi}{si})
            touch{mi}{si}{ti} = union(wla.trials{wtninds(ti)}.protractionTouchFrames, wla.trials{wtninds(ti)}.retractionTouchFrames);
        end
    end
end

%%
cd(whiskerDir)
save('allTouches.mat', 'touch', 'mice', 'sessions')


%% Drawing all touches

behaviorDir = 'Y:\Whiskernas\JK\SoloData\';
whiskerDir = 'Y:\Whiskernas\JK\whisker\tracked\';
videoFreq = 311.24;
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
learner = [1,2,3,4,7,9];
nonlearner = [5,6,8,10,11,12];
% sessions = {[1:19,22],[1:14],[1:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
sessions = {[1:18],[1:12],[1:7,9:20],[1:16],[1:8,10:23],[1:30],[1:21],[1:29],[1:20],[1:20],[1:24],[1:13]};

maxLengthTime = 4; % in sec
maxLength = ceil(videoFreq*maxLengthTime);

for mi = 1 : length(mice)
% for mi = 5 
    totalNumTrialsBySessions = cumsum(cellfun(@(x) length(x),trialNums{mi}(sessions{mi})));
    totalNumTrialsBySessions = [0; totalNumTrialsBySessions];
    touchMat = zeros(totalNumTrialsBySessions(end),maxLength)+0.7;
    for si = 1 : length(sessions{mi})
        for ti = 1 : length(touch{mi}{sessions{mi}(si)})
            currTouch = touch{mi}{sessions{mi}(si)}{ti};
            row = totalNumTrialsBySessions(si)+ti;
            touchMat(row,currTouch(currTouch<=maxLength)) = deal(0);
        end
    end
    figure, imshow(touchMat),title(sprintf('JK%03d',mice(mi)))
end


%% before and after training (first 3 and last 3)
behaviorDir = 'Y:\Whiskernas\JK\SoloData\';
whiskerDir = 'Y:\Whiskernas\JK\whisker\tracked\';
videoFreq = 311.24;
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
learner = [1,2,3,4,7,9];
nonlearner = [5,6,8,10,11,12];
sessions = {[1:19,22],[1:14],[1:7,9:22],[1:18],[1:10,12:24],[1:22,24:31],[1:25],[1:19,21:30],[1:26],[1:3,5:21],[1:24],[1:13]};
maxLengthTime = 4; % in sec
maxLength = ceil(maxLengthTime * videoFreq);
% first touches timepoints
naiveFT = zeros(maxLength,length(learner));
expertFT = zeros(maxLength,length(learner));
firstFT = zeros(maxLength,length(nonlearner));
lastFT = zeros(maxLength,length(nonlearner));
% # of touches
naiveNumT = zeros(maxLength,length(nonlearner));
expertNumT = zeros(maxLength,length(nonlearner));
firstNumT = zeros(maxLength,length(nonlearner));
lastNumT = zeros(maxLength,length(nonlearner));
naiveAmplitude = zeros(maxLength,length(nonlearner));
expertAmplitude = zeros(maxLength,length(nonlearner));
firstAmplitude = zeros(maxLength,length(nonlearner));
lastAmplitude = zeros(maxLength,length(nonlearner));

%%
for mi = 1 : length(learner)
    for si = 1 : 3
        temp = theta{learner(mi)}{si};
        tempTheta = nan(maxLength, length(temp));
        tempMidpoint = nan(maxLength, length(temp));
        tempAmplitude = nan(maxLength, length(temp));
        for i = 1 : length(temp)
            indLength = min(maxLength, length(temp{i}));
            tempTheta(1:indLength,i) = temp{i}(1:indLength);
            [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
        end
        naiveTheta(:,mi) = naiveTheta(:,mi) + nanmean(tempTheta,2)/3;
        naiveMidpoint(:,mi) = naiveMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
        naiveAmplitude(:,mi) = naiveAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
    end
    
    for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
        temp = theta{learner(mi)}{si};
        tempTheta = nan(maxLength, length(temp));
        for i = 1 : length(temp)
            indLength = min(maxLength, length(temp{i}));
            tempTheta(1:indLength, i) = temp{i}(1:indLength);
            [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
        end
        expertTheta(:,mi) = expertTheta(:,mi) + nanmean(tempTheta,2)/3;
        expertMidpoint(:,mi) = expertMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
        expertAmplitude(:,mi) = expertAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
    end
end

for mi = 1 : length(nonlearner)
    for si = 1 : 3
        temp = theta{nonlearner(mi)}{si};
        tempTheta = nan(maxLength, length(temp));
        for i = 1 : length(temp)
            indLength = min(maxLength, length(temp{i}));
            tempTheta(1:indLength, i) = temp{i}(1:indLength);
            [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
        end
        firstTheta(:,mi) = firstTheta(:,mi) + nanmean(tempTheta,2)/3;
        firstMidpoint(:,mi) = firstMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
        firstAmplitude(:,mi) = firstAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
    end
    
    for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
        temp = theta{nonlearner(mi)}{si};
        tempTheta = nan(maxLength, length(temp));
        for i = 1 : length(temp)
            indLength = min(maxLength, length(temp{i}));
            tempTheta(1:indLength, i) = temp{i}(1:indLength);
            [~, tempAmplitude(1:indLength,i), ~, tempMidpoint(1:indLength,i)] = jkWhiskerDecomposition(tempTheta(1:indLength,i));
        end
        lastTheta(:,mi) = lastTheta(:,mi) + nanmean(tempTheta,2)/3;
        lastMidpoint(:,mi) = lastMidpoint(:,mi) + nanmean(tempMidpoint,2)/3;
        lastAmplitude(:,mi) = lastAmplitude(:,mi) + nanmean(tempAmplitude,2)/3;
    end
end
%%
c = parula(10);
c1 = c(4,:);
c2 = c(8,:);
times = 0:1/videoFreq:1/videoFreq*(maxLength-1);
mirrorAngleOffset = 13; % in degrees
figure, 
plot(times, mean(naiveTheta,2)+mirrorAngleOffset, '--', 'color', c1, 'linewidth', 3);
hold on,
plot(times, mean(expertTheta,2)+mirrorAngleOffset, 'color', c1, 'linewidth', 3);
plot(times, mean(firstTheta,2)+mirrorAngleOffset, '--', 'color', c2, 'linewidth', 3);
plot(times, mean(lastTheta,2)+mirrorAngleOffset, 'color', c2, 'linewidth', 3);

boundedline(times, mean(naiveTheta,2)+mirrorAngleOffset, std(naiveTheta, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
boundedline(times, mean(expertTheta,2)+mirrorAngleOffset, std(expertTheta, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
boundedline(times, mean(firstTheta,2)+mirrorAngleOffset, std(firstTheta, [], 2), '--', 'cmap', c2, 'alpha', 'transparency', 0.2);
boundedline(times, mean(lastTheta,2)+mirrorAngleOffset, std(lastTheta, [], 2), 'cmap', c2, 'alpha', 'transparency', 0.2);

plot(times, mean(naiveTheta+mirrorAngleOffset,2), '--', 'color', c1, 'linewidth', 3);
plot(times, mean(expertTheta+mirrorAngleOffset,2), 'color', c1, 'linewidth', 3);
plot(times, mean(firstTheta+mirrorAngleOffset,2), '--', 'color', c2, 'linewidth', 3);
plot(times, mean(lastTheta+mirrorAngleOffset,2), 'color', c2, 'linewidth', 3);
plot([1 1], [-25 5], 'k--')
plot([2 2], [-25 5], 'k--')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold', 'box', 'off')
title('Average whisker angle'), ylabel('Theta (\circ)'), ylim([-25 5]), xlim([0 4]), xlabel('Time (s)')
legend({'Learned (first 3 sessions)'; 'Learned (last 3 sessions)'; 'Not-learned (first 3 sessions)'; 'Not-learned (last 3 sessions)'}, 'box', 'off')


