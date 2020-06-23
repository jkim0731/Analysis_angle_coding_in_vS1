% About information accumulation.
% To figure out if mice might have discriminated angles with multiple touches.
% Are they accumulating information from multiple touches to discriminate angles?
% For example, accumulating touch points along the pole.

%% (1) Number of touch before first lick - line 17
%% Results: Mean # of touch are 4-5 before first lick, in both naive and expert mice - line 102
%% What about if I do not divide by whisking? - line 104
%% Results: similar to "byWhisking", just a little bit less. - line 187

%% (2) Divide into different # of touches before first lick, and look at the correct rate. - line 190
%% Since # of touch were not different between angles, just use the all angles, except for 90 degrees.
%% Also test in 2 angles expert and radial distance test.
%% Also, compare with non-touch trials (did mice guess correctly without touch? Possible with whisking only)

%% First, by whisking - line 195
%% Then, by touch chunks - line 255

%% Results: need at least 2 touches to answer correctly in 7 angles, but in 2 angles and radial distance, they could guess correctly even with no touch before the first lick.
%% This is be possible by 1) whisking downward, where touch happens only for 135 degrees but not for 45 degrees pole most of the time (touch trial proporion is lower in 45 degrees)
%% and 2) using lick-induced whisking to touch the pole. There are behavioral variations between mice. - line 324

%% What about before answer lick? - line 328
%% Use touch chunks from now on, since there were not so much different between touch chunks and bywhisking, and instantaneous touch feature makes more sense with touch chunks.

%% Results: Mice need at least 3 touches in general, but their number of trials with <= 3 touches are low. 
%% Some mice simply do not have enough number of trials in this bin (at least 2 - JK036 & Jk039) - line 418

% -> Moved on to d190805_touch_evidence_accumulation. 




%% (1) Number of touch before first lick (among the answered trials)
% baseDir = 'C:\Data\suite2p\';
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
% expertInds = [1,2,3,9]; % excluding 36 and 39
angles = 45:15:135;
numTouch = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Only when there's touch
meanTouchNum = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Including non-touch trials
propTouchTrial = zeros(length(expertInds),length(angles),2); % Pre-decision. Proportion of trials with touch, when the animal responded.
for ei = 1 : length(expertInds)
% for ei = 1
    mouse = mice(expertInds(ei));
    for si = 1 : 2 % 1 for naive, 2 for expert
%     for si = 1
        session = sessions{expertInds(ei)}(si);
        ufn = sprintf('UberJK%03dS%02d', mouse, session);
        dn = sprintf('%s%03d\\',baseDir,mouse);
        load([dn,ufn],'u')
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), u.trials));
        protractionTouchTrialInds = find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
        touchTrialInds = intersect(answerTrialInds, protractionTouchTrialInds);        
        
        for ai = 1 : length(angles)
            angleInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            angleTouchTrialInds = intersect(angleInds, touchTrialInds);
            allLickTime = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(angleTouchTrialInds), 'uniformoutput', false);
            firstLickTime = cellfun(@(x,y) x(find(x>y.poleUpOnsetTime,1)), allLickTime, u.trials(angleTouchTrialInds), 'uniformoutput', false);
            numTouch(ei,ai,si) = mean(cellfun(@(x,z) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= z)), x.protractionTouchChunksByWhisking)), u.trials(angleTouchTrialInds), firstLickTime));
            
            angleAnswerTrialInds = intersect(angleInds, answerTrialInds);
            
            allLickTime = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(angleAnswerTrialInds), 'uniformoutput', false);
            firstLickTime = cellfun(@(x,y) x(find(x>y.poleUpOnsetTime,1)), allLickTime, u.trials(angleAnswerTrialInds), 'uniformoutput', false);
            
            meanTouchNum(ei,ai,si) = mean(cellfun(@(x,z) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= z)), x.protractionTouchChunksByWhisking)), u.trials(angleAnswerTrialInds), firstLickTime));
            propTouchTrial(ei,ai,si) = length(find(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunksByWhisking)), u.trials(angleAnswerTrialInds)))) / length(angleAnswerTrialInds);
        end        
    end
end

%%

figure, 
subplot(231), hold on
shadedErrorBar(angles, mean(numTouch(:,:,1)), std(numTouch(:,:,1))/sqrt(size(numTouch,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(numTouch(:,:,2)), std(numTouch(:,:,2))/sqrt(size(numTouch,1)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(from touch trials only)'})

subplot(232), hold on
shadedErrorBar(angles, mean(meanTouchNum(:,:,1)), std(meanTouchNum(:,:,1))/sqrt(size(meanTouchNum,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(meanTouchNum(:,:,2)), std(meanTouchNum(:,:,2))/sqrt(size(meanTouchNum,1)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(including non-touch trials)'})

subplot(233), hold on
shadedErrorBar(angles, mean(propTouchTrial(:,:,1)), std(propTouchTrial(:,:,1))/sqrt(size(propTouchTrial,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(propTouchTrial(:,:,2)), std(propTouchTrial(:,:,2))/sqrt(size(propTouchTrial,1)), 'lineprop', 'r-')
xticks(angles), title({'Proportion of touch trials'; '(from all answered trials)'})
legend({'Naive', 'Expert'}, 'location', 'southeast', 'box', 'off')
% Too much variation between mice?
% look at the difference between expert and naive

subplot(234), hold on
temp = numTouch;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')

subplot(235), hold on
temp = meanTouchNum;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')

subplot(236), hold on
temp = propTouchTrial;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')
legend('Expert - Naive', 'location', 'northeast', 'box', 'off')


%% Results: Mean # of touch are 4-5 before first lick, in both naive and expert mice

%% What about if I do not divide by whisking?
% baseDir = 'C:\Data\suite2p\';
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
% expertInds = [1,2,3,9]; % excluding 36 and 39
angles = 45:15:135;
numTouch = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Only when there's touch
meanTouchNum = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Including non-touch trials
propTouchTrial = zeros(length(expertInds),length(angles),2); % Pre-decision. Proportion of trials with touch, when the animal responded.
for ei = 1 : length(expertInds)
% for ei = 1
    mouse = mice(expertInds(ei));
    for si = 1 : 2 % 1 for naive, 2 for expert
%     for si = 1
        session = sessions{expertInds(ei)}(si);
        ufn = sprintf('UberJK%03dS%02d', mouse, session);
        dn = sprintf('%s%03d\\',baseDir,mouse);
        load([dn,ufn],'u')
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), u.trials));
        protractionTouchTrialInds = find(cellfun(@(x) length(x.protractionTouchChunks), u.trials));
        touchTrialInds = intersect(answerTrialInds, protractionTouchTrialInds);        
        
        for ai = 1 : length(angles)
            angleInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            angleTouchTrialInds = intersect(angleInds, touchTrialInds);
            allLickTime = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(angleTouchTrialInds), 'uniformoutput', false);
            firstLickTime = cellfun(@(x,y) x(find(x>y.poleUpOnsetTime,1)), allLickTime, u.trials(angleTouchTrialInds), 'uniformoutput', false);
            numTouch(ei,ai,si) = mean(cellfun(@(x,z) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= z)), x.protractionTouchChunks)), u.trials(angleTouchTrialInds), firstLickTime));
            
            angleAnswerTrialInds = intersect(angleInds, answerTrialInds);
            
            allLickTime = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(angleAnswerTrialInds), 'uniformoutput', false);
            firstLickTime = cellfun(@(x,y) x(find(x>y.poleUpOnsetTime,1)), allLickTime, u.trials(angleAnswerTrialInds), 'uniformoutput', false);
            
            meanTouchNum(ei,ai,si) = mean(cellfun(@(x,z) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= z)), x.protractionTouchChunks)), u.trials(angleAnswerTrialInds), firstLickTime));
            propTouchTrial(ei,ai,si) = length(find(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunks)), u.trials(angleAnswerTrialInds)))) / length(angleAnswerTrialInds);
        end        
    end
end

%
figure, 
subplot(231), hold on
shadedErrorBar(angles, mean(numTouch(:,:,1)), std(numTouch(:,:,1))/sqrt(size(numTouch,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(numTouch(:,:,2)), std(numTouch(:,:,2))/sqrt(size(numTouch,1)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(from touch trials only)'})

subplot(232), hold on
shadedErrorBar(angles, mean(meanTouchNum(:,:,1)), std(meanTouchNum(:,:,1))/sqrt(size(meanTouchNum,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(meanTouchNum(:,:,2)), std(meanTouchNum(:,:,2))/sqrt(size(meanTouchNum,1)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(including non-touch trials)'})

subplot(233), hold on
shadedErrorBar(angles, mean(propTouchTrial(:,:,1)), std(propTouchTrial(:,:,1))/sqrt(size(propTouchTrial,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(propTouchTrial(:,:,2)), std(propTouchTrial(:,:,2))/sqrt(size(propTouchTrial,1)), 'lineprop', 'r-')
xticks(angles), title({'Proportion of touch trials'; '(from all answered trials)'})
legend({'Naive', 'Expert'}, 'location', 'southeast', 'box', 'off')
% Too much variation between mice?
% look at the difference between expert and naive

subplot(234), hold on
temp = numTouch;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')

subplot(235), hold on
temp = meanTouchNum;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')

subplot(236), hold on
temp = propTouchTrial;
tempDiff = temp(:,:,2) - temp(:,:,1);
shadedErrorBar(angles, mean(tempDiff), std(tempDiff), 'lineprop', 'k-')
xticks(angles), xlabel('Angle (\circ)')
legend('Expert - Naive', 'location', 'northeast', 'box', 'off')

%% Results: similar to "byWhisking", just a little bit less. 


%% (2) Divide into different # of touches before first lick, and look at the correct rate.
%% Since # of touch were not different between angles, just use the all angles, except for 90 degrees.
%% Also test in 2 angles expert and radial distance test.
%% Also, compare with non-touch trials (did mice guess correctly without touch? Possible with whisking only)

%% First, by whisking
wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
bDir = 'Y:\Whiskernas\JK\SoloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
numTouchBin = [0:7,100]; % 7th bin will have # of touches >= 7 
% correctRate = zeros(length(mice),length(numTouchBin)-1,3);
numTrial = zeros(length(mice),length(numTouchBin)-1,3);
numCorrect = zeros(length(mice),length(numTouchBin)-1,3);
for mi = 1 : length(mice)
% for ei = 1
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%     for si = 1
        session = sessions{mi}(si);
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));
        angleTrialInds = find(cellfun(@(x) x.poleAngle ~= 90, wfa.trials));
        testInds = intersect(answerTrialInds, angleTrialInds);
        
        testB = wB(testInds);
        testW = wfa.trials(testInds);
        
        firstLickTime = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), testB, 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)) < y)), testW, firstLickTime);
        correctInd = find(cellfun(@(x) x.trialCorrect, testB));

        numTrial(mi,:,si) = histcounts(numTouch,numTouchBin);
        numCorrect(mi,:,si) = histcounts(numTouch(correctInd),numTouchBin);        

    end
end
%
correctRate = numCorrect./numTrial;
propTrial = numTrial./sum(numTrial,2);
titles = {'2 angles', '7 angles', 'Radial distance'};
figure,
for i = 1 : 3    
    subplot(2,3,i), errorbar(0:7, mean(correctRate(:,:,i)), std(correctRate(:,:,i))/sqrt(length(mice)))
    title(titles(i))
    xlim([-1 8])
    xticks([0:7])
    xticklabels({'0','1','2','3','4','5','6','>7'})
    ylabel('Correct rate')
    ylim([0.4 1])
    subplot(2,3,i+3), errorbar(0:7, mean(propTrial(:,:,i)), std(propTrial(:,:,i))/sqrt(length(mice)))
    xlim([-1 8])
    xticks([0:7])
    xticklabels({'0','1','2','3','4','5','6','>7'})
    ylabel('Proportion of trials')
    ylim([0 0.4])
    if i == 2
        xlabel('# of touch before first lick')
    end
end

%% Then, by touch chunks
wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
bDir = 'Y:\Whiskernas\JK\SoloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
numTouchBin = [0:7,100]; % 7th bin will have # of touches >= 7 
% correctRate = zeros(length(mice),length(numTouchBin)-1,3);
numTrial = zeros(length(mice),length(numTouchBin)-1,3);
numCorrect = zeros(length(mice),length(numTouchBin)-1,3);
for mi = 1 : length(mice)
% for ei = 1
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%     for si = 1
        session = sessions{mi}(si);
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));
        angleTrialInds = find(cellfun(@(x) x.poleAngle ~= 90, wfa.trials));
        testInds = intersect(answerTrialInds, angleTrialInds);
        
        testB = wB(testInds);
        testW = wfa.trials(testInds);
        
        firstLickTime = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), testB, 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), testW, firstLickTime);
        correctInd = find(cellfun(@(x) x.trialCorrect, testB));

        numTrial(mi,:,si) = histcounts(numTouch,numTouchBin);
        numCorrect(mi,:,si) = histcounts(numTouch(correctInd),numTouchBin);        

    end
end
%
correctRate = numCorrect./numTrial;
propTrial = numTrial./sum(numTrial,2);
titles = {'2 angles', '7 angles', 'Radial distance'};
figure,
for i = 1 : 3    
    subplot(2,3,i), errorbar(0:7, mean(correctRate(:,:,i)), std(correctRate(:,:,i))/sqrt(length(mice)))
    title(titles(i))
    xlim([-1 8])
    xticks([0:7])
    xticklabels({'0','1','2','3','4','5','6','>7'})
    ylabel('Correct rate')
    ylim([0.4 1])
    subplot(2,3,i+3), errorbar(0:7, mean(propTrial(:,:,i)), std(propTrial(:,:,i))/sqrt(length(mice)))
    xlim([-1 8])
    xticks([0:7])
    xticklabels({'0','1','2','3','4','5','6','>7'})
    ylabel('Proportion of trials')
    ylim([0 0.4])
    if i == 2
        xlabel('# of touch before first lick')
    end
end

%% Results: need at least 2 touches to answer correctly in 7 angles, but in 2 angles and radial distance, they could guess correctly even with no touch before the first lick.
%% This is be possible by 1) whisking downward, where touch happens only for 135 degrees but not for 45 degrees pole most of the time (touch trial proporion is lower in 45 degrees)
%% and 2) using lick-induced whisking to touch the pole. There are behavioral variations between mice. 

%% What about before answer lick?
%% Use touch chunks from now on, since there were not so much different between touch chunks and bywhisking, and instantaneous touch feature makes more sense with touch chunks.

wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
bDir = 'Y:\Whiskernas\JK\SoloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
numTouchBin = [0:10,100]; % 10th bin will have # of touches >= 10 
% correctRate = zeros(length(mice),length(numTouchBin)-1,3);
numTrial = zeros(length(mice),length(numTouchBin)-1,3);
numCorrect = zeros(length(mice),length(numTouchBin)-1,3);
for mi = 1 : length(mice)
% for ei = 1
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%     for si = 1
        session = sessions{mi}(si);
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));
        angleTrialInds = find(cellfun(@(x) x.poleAngle ~= 90, wfa.trials));
        testInds = intersect(answerTrialInds, angleTrialInds);
        
        testB = wB(testInds);
        testW = wfa.trials(testInds);
        
        answerLickTime = cellfun(@(x) x.answerLickTime, testB, 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), testW, answerLickTime);
        correctInd = find(cellfun(@(x) x.trialCorrect, testB));

        numTrial(mi,:,si) = histcounts(numTouch,numTouchBin);
        numCorrect(mi,:,si) = histcounts(numTouch(correctInd),numTouchBin);        

    end
end
%%
totalCR = squeeze(sum(numCorrect,2) ./ sum(numTrial,2));
correctRate = numCorrect./numTrial;
pvalChance = zeros(size(correctRate,3),size(correctRate,2));
pvalBest = zeros(size(correctRate,3),size(correctRate,2));
for i = 1 : size(correctRate,3)
    [~,pvalChance(i,:)] = ttest(correctRate(:,:,i)-0.5);
    [~,pvalBest(i,:)] = ttest(correctRate(:,:,i) - totalCR(:,i));
end
%%
propTrial = numTrial./sum(numTrial,2);
titles = {'2 angles', '7 angles', 'Radial distance'};
figure,
for i = 1 : 3    
    subplot(2,3,i), errorbar(0:10, nanmean(correctRate(:,:,i)), nanstd(correctRate(:,:,i))/sqrt(length(mice)), '-k'), hold on
    tempInd = find(pvalChance(i,:)<0.01);
    xval = 0:10;
    errorbar(xval(tempInd), nanmean(correctRate(:,tempInd,i)), nanstd(correctRate(:,tempInd,i))/sqrt(length(mice)), '.r')
    title(titles(i))
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Correct rate')
    ylim([0.4 1])
    subplot(2,3,i+3), errorbar(0:10, mean(propTrial(:,:,i)), std(propTrial(:,:,i))/sqrt(length(mice)), '-k')
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Proportion of trials')
    ylim([0 0.4])
    if i == 2
        xlabel('# of touch before answer lick')
    end
end

%% spread of num trials in each category
titles = {'2 angles', '7 angles', 'Radial distance'};
figure, 
for i = 1 : 3
    subplot(1,3,i), errorbar(0:10, mean(numTrial(:,:,i)), std(numTrial(:,:,i)))
    title(titles(i))
    xlim([-1 11]), xticks([0:2:10]), xticklabels({'0','2','4','6','8','>10'})
    xlabel('# of touch before answer lick')
    ylabel('# of trials')
end


%% Results: Mice need at least 3 touches in general, but their number of trials with <= 3 touches are low. 
%% Some mice simply do not have enough number of trials in this bin (at least 2 - JK036 & Jk039)






%% Re-visited 2019/12/09
% look at having jsut 1 touch, lasting at least 10 ms, before the answer lick

wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
bDir = 'Y:\Whiskernas\JK\SoloData\';
mice = [25,27,30,36,39,52];
% sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
numTouchBin = [0:10,100]; % 10th bin will have # of touches >= 10 
% correctRate = zeros(length(mice),length(numTouchBin)-1,3);
numTrial = zeros(length(mice),length(numTouchBin)-1,3);
numCorrect = zeros(length(mice),length(numTouchBin)-1,3);
for mi = 1 : length(mice)
% for ei = 1
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%     for si = 1
        session = sessions{mi}(si);
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));
        angleTrialInds = find(cellfun(@(x) x.poleAngle ~= 90, wfa.trials));
        touchTrialInds = find(cellfun(@(x) length(x.protractionTouchDuration), wfa.trials));
        testInds = intersect(intersect(answerTrialInds, angleTrialInds), touchTrialInds);
        
        testBtemp = wB(testInds);
        testWtemp = wfa.trials(testInds);
        
        firstTouchLong = find(cellfun(@(x) x.protractionTouchDuration(1) > 0.015, testWtemp));
        
        testB = testBtemp(firstTouchLong);
        testW = testWtemp(firstTouchLong);
        
        answerLickTime = cellfun(@(x) x.answerLickTime, testB, 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), testW, answerLickTime);
        correctInd = find(cellfun(@(x) x.trialCorrect, testB));

        numTrial(mi,:,si) = histcounts(numTouch,numTouchBin);
        numCorrect(mi,:,si) = histcounts(numTouch(correctInd),numTouchBin);        

    end
end
%%
totalCR = squeeze(sum(numCorrect,2) ./ sum(numTrial,2));
correctRate = numCorrect./numTrial;
pvalChance = zeros(size(correctRate,3),size(correctRate,2));
pvalBest = zeros(size(correctRate,3),size(correctRate,2));
for i = 1 : size(correctRate,3)
    [~,pvalChance(i,:)] = ttest(correctRate(:,:,i)-0.5);
    [~,pvalBest(i,:)] = ttest(correctRate(:,:,i) - totalCR(:,i));
end
%%
propTrial = numTrial./sum(numTrial,2);
titles = {'2 angles', '7 angles', 'Radial distance'};
figure,
for i = 1 : 3    
    subplot(2,3,i), errorbar(0:10, nanmean(correctRate(:,:,i)), nanstd(correctRate(:,:,i))/sqrt(length(mice)), '-k'), hold on
    tempInd = find(pvalChance(i,:)<0.01);
    xval = 0:10;
    errorbar(xval(tempInd), nanmean(correctRate(:,tempInd,i)), nanstd(correctRate(:,tempInd,i))/sqrt(length(mice)), '.r')
    title(titles(i))
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Correct rate')
    ylim([0.4 1])
    subplot(2,3,i+3), errorbar(0:10, mean(propTrial(:,:,i)), std(propTrial(:,:,i))/sqrt(length(mice)), '-k')
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Proportion of trials')
    ylim([0 0.4])
    if i == 2
        xlabel('# of touch before answer lick')
    end
end


%% Results: it's still the same. They need at least 2 touches for better than chance performance, and more than 3 touches perform similarly.
% what about before the first lick? maybe many of 1 touch before answer
% lick were impulsive lickings, start licking even before touches


%%
wDir = 'Y:\Whiskernas\JK\whisker\tracked\';
bDir = 'Y:\Whiskernas\JK\SoloData\';
mice = [25,27,30,36,39,52];
% sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 
numTouchBin = [0:10,100]; % 10th bin will have # of touches >= 10 
% correctRate = zeros(length(mice),length(numTouchBin)-1,3);
numTrial = zeros(length(mice),length(numTouchBin)-1,3);
numCorrect = zeros(length(mice),length(numTouchBin)-1,3);
for mi = 1 : length(mice)
% for ei = 1
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%     for si = 1
        session = sessions{mi}(si);
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));
        angleTrialInds = find(cellfun(@(x) x.poleAngle ~= 90, wfa.trials));
        touchTrialInds = find(cellfun(@(x) length(x.protractionTouchDurationByWhisking), wfa.trials));
        testInds = intersect(intersect(answerTrialInds, angleTrialInds), touchTrialInds);
        
        testBtemp = wB(testInds);
        testWtemp = wfa.trials(testInds);
        
        firstLickTime = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), testBtemp, 'uniformoutput', false); % no empty cell because these are from answer trials
        
        indFirstTouchBeforeLick = find(cellfun(@(x,y) x.time(x.protractionTFchunksByWhisking{1}(end)) < y, testWtemp, firstLickTime));
        indFirstTouchLong = find(cellfun(@(x) x.protractionTouchDuration(1) > 0.015, testWtemp));
        
        testB = testBtemp(intersect(indFirstTouchBeforeLick, indFirstTouchLong));
        testW = testWtemp(intersect(indFirstTouchBeforeLick, indFirstTouchLong));
        
        answerLickTime = cellfun(@(x) x.answerLickTime, testB, 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)) < y)), testW, firstLickTime(intersect(indFirstTouchBeforeLick, indFirstTouchLong)));
        correctInd = find(cellfun(@(x) x.trialCorrect, testB));

        numTrial(mi,:,si) = histcounts(numTouch,numTouchBin);
        numCorrect(mi,:,si) = histcounts(numTouch(correctInd),numTouchBin);        

    end
end
%
totalCR = squeeze(sum(numCorrect,2) ./ sum(numTrial,2));
correctRate = numCorrect./numTrial;
pvalChance = zeros(size(correctRate,3),size(correctRate,2));
pvalBest = zeros(size(correctRate,3),size(correctRate,2));
for i = 1 : size(correctRate,3)
    [~,pvalChance(i,:)] = ttest(correctRate(:,:,i)-0.5);
    [~,pvalBest(i,:)] = ttest(correctRate(:,:,i) - totalCR(:,i));
end
%
propTrial = numTrial./sum(numTrial,2);
titles = {'2 angles', '7 angles', 'Radial distance'};
figure,
for i = 1 : 3    
    subplot(2,3,i), errorbar(0:10, nanmean(correctRate(:,:,i)), nanstd(correctRate(:,:,i))/sqrt(length(mice)), '-k'), hold on
    tempInd = find(pvalChance(i,:)<0.01);
    xval = 0:10;
    errorbar(xval(tempInd), nanmean(correctRate(:,tempInd,i)), nanstd(correctRate(:,tempInd,i))/sqrt(length(mice)), '.r')
    title(titles(i))
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Correct rate')
    ylim([0.4 1])
    subplot(2,3,i+3), errorbar(0:10, mean(propTrial(:,:,i)), std(propTrial(:,:,i))/sqrt(length(mice)), '-k')
    xlim([-1 11])
    xticks([0:2:10])
    xticklabels({'0','2','4','6','8','>10'})
    ylabel('Proportion of trials')
    ylim([0 0.4])
    if i == 2
        xlabel('# of touch before answer lick')
    end
end

%% Results: well now, in 2 angles and radial distance session, just 1 touch makes it perform better than change (> 70 %)

