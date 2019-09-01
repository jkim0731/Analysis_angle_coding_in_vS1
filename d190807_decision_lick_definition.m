% Purpose: To see how many touches are there (either by chunks and by whisking)
% Distribution of # of touch before (1) first lick, (2) decision lick, (3) answer lick
% Before that, first need to define decision lick. 
% Even before that, need to see how to define decision lick.

%% (1) How much do mice alternate between first lick and answer lick?
%% First show first lick time and answer lick time.
%% Then look for time difference, # of licks inbetween, # of alternation, 
%% Compared to correct rate, angles (especially in 7 angles trials)
%% See distribution, not the mean values. In each mouse.

bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

answerLickTime = cell(length(mice),3);
firstLickTime = cell(length(mice),3);
lickOverlapTrials = cell(length(mice),3);

lickPatternGroup = cell(length(mice),3);
lickTimeGroup = cell(length(mice),3);

numLick = cell(length(mice),3);
numAlt = cell(length(mice),3);
correct = cell(length(mice),3);
objectAngle = cell(length(mice),3);
for mi = 1 : length(mice)
%%
% for mi = 4:6
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%%
%     for si = 2:3
        session = sessions{mi}(si);
        fprintf('Processing JK%03d S%02d\n', mouse, session)
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        
        answerTrialInds = setdiff(find(cellfun(@(x) length(x.answerLickTime), currB.trials)), find(cellfun(@(x) strcmp(x.trialType, 'oo'), currB.trials)));
        answerLickTime{mi,si} = cellfun(@(x) x.answerLickTime, currB.trials(answerTrialInds)); % no empty cell because these are from answer trials
        firstLickTime{mi,si} = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), currB.trials(answerTrialInds)); % no empty cell because these are from answer trials
        
        lickPattern = cell(1,length(answerTrialInds));
        lickTime = cell(1,length(answerTrialInds));
        lickOverlap = zeros(1,length(answerTrialInds));
        parfor ti = 1 : length(answerTrialInds)
            tempB = currB.trials{answerTrialInds(ti)};
            lickTimes = union(tempB.beamBreakTimes, tempB.answerLickTime);
            firstInd = find(lickTimes > tempB.poleUpOnsetTime,1);
            answerInd = find(lickTimes == tempB.answerLickTime);
            lickTime{ti} = lickTimes(firstInd:answerInd);
            tempLicks = zeros(1,length(lickTimes(firstInd:answerInd)));
            leftLicks = ismember(lickTimes(firstInd:answerInd), tempB.beamBreakTimesLeft);
            rightLicks = ismember(lickTimes(firstInd:answerInd), tempB.beamBreakTimesRight); 
            tempLicks(leftLicks) = 1;
            tempLicks(rightLicks) = 2; % there can be overlaps, but for now, just treat them to R lick. Count the trials with overlapping liks though.
            if tempB.trialType(1) == 'l'
                if tempB.trialCorrect == 1
                    tempLicks(end) = 1;
                elseif tempB.trialCorrect == 0
                    tempLicks(end) = 2;
                else
                    error('no answer trial, #%d',tempB.trialNum)
                end
            elseif tempB.trialType(1) == 'r'
                if tempB.trialCorrect == 1
                    tempLicks(end) = 2;
                elseif tempB.trialCorrect == 0
                    tempLicks(end) = 1;
                else
                    error('no answer trial, #%d',tempB.trialNum)
                end
            end
            if find(tempLicks == 0)
                error('There is no lick side determined in trial #%d', tempB.trialNum)
            end
            lickPattern{ti} = tempLicks;
            lickOverlap(ti) = length(intersect(leftLicks, rightLicks));            
        end
        lickOverlapTrials{mi,si} = lickOverlap;
        
        lickPatternGroup{mi, si} = lickPattern;
        lickTimeGroup{mi, si} = lickTime;
        
        numLick{mi,si} = cellfun(@(x) length(x)-1, lickPattern);
        numAlt{mi,si} = cellfun(@(x) sum(abs(diff(x))), lickPattern);
        
        correct{mi,si} = cellfun(@(x) x.trialCorrect, currB.trials(answerTrialInds));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, currB.trials(answerTrialInds));
    end
end

%% distribution of lick times
for si = 1 : 3
    figure,
    lickRange = [1:0.1:5.1];
    for mi = 1 : length(mice)
        subplot(2,3,mi), hold on
        plot(lickRange(1:end-1), histcounts(answerLickTime{mi,si}, lickRange))
        plot(lickRange(1:end-1), histcounts(firstLickTime{mi,si}, lickRange))
        xlim([1 5])
        if mi == 5
            xlabel('Time (s)')
        end
        ylabel('Count')
        if mi == 3
            legend({'Answer lick', 'First lick'})
        end
        title(sprintf('JK%03d', mice(mi)))
    end
end

%% distribution of answer-first lick times
figure, 
titleList = {'2 angles', '7 angles', 'Radial distance'};
lickDiffRange = [0:0.1:1.6];
for si = 1 : 3
    subplot(2,3,si), hold on
    for mi = 1 : length(mice)
        plot(lickDiffRange(1:end-1), histcounts(answerLickTime{mi,si} - firstLickTime{mi,si}, lickDiffRange, 'normalization', 'probability'))
    end    
    if si == 3
        legendList = cell(1,length(mice));
        for li = 1 : length(mice)
            legendList{li} = sprintf('JK%03d',mice(li));
        end
        legend(legendList)
    end
    if si == 1
        ylabel('Proportion')
    end
    title(titleList{si})
    
    subplot(2,3,si+3), hold on
    for mi = 1 : length(mice)
        plot(lickDiffRange(1:end-1), histcounts(answerLickTime{mi,si} - firstLickTime{mi,si}, lickDiffRange, 'normalization', 'cdf'))
    end
    plot(lickDiffRange(1:end-1), ones(1,length(lickDiffRange)-1) * 0.5, 'k--')
    xlabel('Time (s)')
    if si == 1
        ylabel('Cummulative proportion')
    end
end

%% distribution of # of licks before the answer lick
figure,
titleList = {'2 angles', '7 angles', 'Radial distance'};
lickRange = [0:10,20];
for si = 1 : 3
    subplot(2,3,si), hold on
    for mi = 1 : length(mice)
        plot(lickRange(1:end-1), histcounts(numLick{mi,si}, lickRange, 'normalization', 'probability'))
    end
    if si == 3
        legendList = cell(1,length(mice));
        for li = 1 : length(mice)
            legendList{li} = sprintf('JK%03d',mice(li));
        end
        legend(legendList)
    end
    if si == 1
        ylabel('Proportion')
    end
    title(titleList{si})
    
    subplot(2,3,si+3), hold on
    for mi = 1 : length(mice)
        plot(lickRange(1:end-1), histcounts(numLick{mi,si}, lickRange, 'normalization', 'cdf'))
    end
    plot(lickRange(1:end-1), ones(1,length(lickRange)-1) * 0.5, 'k--')
    xlabel('# of licks before the answer lick')
    if si == 1
        ylabel('Cummulative proportion')
    end
end

%% distribution of # of alternating licks till the answer lick
figure,
titleList = {'2 angles', '7 angles', 'Radial distance'};
altRange = [0:4,10];
for si = 1 : 3
    subplot(2,3,si), hold on
    for mi = 1 : length(mice)
        plot(altRange(1:end-1), histcounts(numAlt{mi,si}, altRange, 'normalization', 'probability'))
    end
    if si == 3
        legendList = cell(1,length(mice));
        for li = 1 : length(mice)
            legendList{li} = sprintf('JK%03d',mice(li));
        end
        legend(legendList)
    end
    if si == 1
        ylabel('Proportion')
    end
    title(titleList{si})
    
    subplot(2,3,si+3), hold on
    for mi = 1 : length(mice)
        plot(altRange(1:end-1), histcounts(numAlt{mi,si}, altRange, 'normalization', 'cdf'))
    end
    if si == 2
        xlabel('# of alternation until the answer lick')
    end
    if si == 1
        ylabel('Cummulative proportion')
    end
end

%% Correct rate vs alternating licks (0, 1, and 2)
figure,
titleList = {'2 angles', '7 angles', 'Radial distance'};
for si = 1 : 3
    subplot(1,3,si), hold on
    for mi = 1 : length(mice)
        cr = zeros(1,3); % correct rate
        for i = 1 : 3
            altind = find(numAlt{mi,si} == i-1);
            cr(i) = sum(correct{mi,si}(altind)) / length(altind);
        end        
        plot(cr)
    end
    xlim([0 4])
    xticks([1:3])
    xticklabels([0:2])
    xlabel('# of alternation')
    ylim([0.5 1])
    if si == 1
        ylabel('Correct rate')
    end
    title(titleList{si})
    if si == 3
        legendList = cell(1,length(mice));
        for li = 1 : length(mice)
            legendList{li} = sprintf('JK%03d',mice(li));
        end
        legend(legendList)
    end
end
        
    
%% Distribution of first lick vs # of alternating licks (0, 1, and 2)
for si = 1 : 3
    figure,
    lickRange = [1:0.1:5.1];
    for mi = 1 : length(mice)    
        subplot(2,3,mi), hold on    
        for i = 1 : 3
            plot(lickRange(1:end-1), histcounts(firstLickTime{mi,si}(find(numAlt{mi,si}==i-1)), lickRange))
        end
        xlim([1 3]) 
        if mi == 5
            xlabel('Time (s)')
        end
        ylabel('Count')
        if mi == 3
            legend({'0', '1', '2'})
        end
        title(sprintf('JK%03d', mice(mi)))
    end
end

%% Max lick interval VS first lick time
for si = 1 : 3
    figure
    for mi = 1 : length(mice)
        subplot(2,3,mi), hold on
        intervalInd = find(cellfun(@(x) length(x)>1, lickTimeGroup{mi,si}));
        altInd = find(numAlt{mi,si}(intervalInd) == 0);
        tempFirstLick = firstLickTime{mi,si}(intervalInd(altInd));
        tempMaxInterval = cellfun(@(x) max(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
        scatter(tempFirstLick, tempMaxInterval, 10, 'k', 'filled')
%         tempMeanInterval = cellfun(@(x) mean(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
%         scatter(tempFirstLick, tempMeanInterval, 10, 'k', 'filled')
        
        altInd = find(numAlt{mi,si}(intervalInd) > 0);
        tempFirstLick = firstLickTime{mi,si}(intervalInd(altInd));
        tempMaxInterval = cellfun(@(x) max(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
        scatter(tempFirstLick, tempMaxInterval, 10, 'r', 'filled')
%         tempMeanInterval = cellfun(@(x) mean(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
%         scatter(tempFirstLick, tempMeanInterval, 10, 'k', 'filled')        
    end
end
        
%% Distribution of max lick interval in different alternation (can I define a threshold for lick bout interval?)
intervalRange = [0:0.1:2,5];
for si = 1 : 3
    figure
    for mi = 1 : length(mice)
        subplot(2,3,mi), hold on
        intervalInd = find(cellfun(@(x) length(x)>1, lickTimeGroup{mi,si}));
        altInd = find(numAlt{mi,si}(intervalInd) == 0);        
        tempMaxInterval = cellfun(@(x) max(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
        plot(intervalRange(1:end-1), histcounts(tempMaxInterval,intervalRange, 'normalization', 'probability'))
        
        altInd = find(numAlt{mi,si}(intervalInd) > 0);
        tempFirstLick = firstLickTime{mi,si}(intervalInd(altInd));
        tempMaxInterval = cellfun(@(x) mean(diff(x)), lickTimeGroup{mi,si}(intervalInd(altInd)));
        plot(intervalRange(1:end-1), histcounts(tempMaxInterval,intervalRange, 'normalization', 'probability'))
        if mi == 3
            legend({'No alt', '>0 alt'})
        end
        if mi == 5
            xlabel('Time (s)')
        end
        if mi == 1 || mi == 4
            ylabel('Proportion')
        end
        title(sprintf('JK%03d', mice(mi)))
    end
end

%% Distribution of lick interval (can I define a threshold for lick bout interval?)
titleList = {'2 angles', '7 angles', 'Radial distance'};
intervalRange = [0:0.05:2,5];
for si = 1 : 3
    subplot(1,3,si), hold on
    for mi = 1 : length(mice)        
        intervalInd = find(cellfun(@(x) length(x)>1, lickTimeGroup{mi,si}));        
        intervalVal = cell2mat(cellfun(@(x) diff(x)', lickTimeGroup{mi,si}(intervalInd), 'uniformoutput', false));
        plot(intervalRange(1:end-1), histcounts(intervalVal, intervalRange, 'normalization', 'probability'))
    end
    xlabel('Time (s)')
    xlim([0 0.6])
    if si == 1
        ylabel('Proportion')
    end
    if si == 3
        legendList = cell(1,length(mice));
        for li = 1 : length(mice)
            legendList{li} = sprintf('JK%03d',mice(li));
        end
        legend(legendList)
    end
    title(titleList{si})
end

%% Distribution of time difference between alternation
%% Compared to consecutive same-side licking

intervalRange = [0:0.05:2,5];
for si = 1 : 3
    figure,
    for mi = 1 : length(mice)
        subplot(2,3,mi), hold on
        tempInd = find(cellfun(@(x) length(x) > 1, lickTimeGroup{mi,si}));
        altInd = cell2mat(cellfun(@(x) abs(diff(x)) == 1, lickPatternGroup{mi,si}(tempInd), 'uniformoutput', false));
        noAltInd = cell2mat(cellfun(@(x) abs(diff(x)) == 0, lickPatternGroup{mi,si}(tempInd), 'uniformoutput', false));
        altTimeDiff = cell2mat(cellfun(@(x) diff(x)', lickTimeGroup{mi,si}(tempInd), 'uniformoutput', false));
        noAltDist = histcounts(altTimeDiff(noAltInd), intervalRange);
        altDist = histcounts(altTimeDiff(altInd), intervalRange);
        plot(intervalRange(1:end-1), noAltDist, 'k-')
        plot(intervalRange(1:end-1), altDist, 'r-')
        if mi == 3
            legend({'No alt', 'Alt'})
        end
        if mi == 5
            xlabel('Time (s)')
        end
        if mi == 1 || mi == 4
            ylabel('Proportion')
        end
        title(sprintf('JK%03d', mice(mi)))
        xlim([0 0.6])
    end
end











%% Proportion or number of trials with different # of touches before the first lick

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
numTouchBeforeFirstLick = cell(length(mice),3); 
objectAngle = cell(length(mice),3);
for mi = 1 : length(mice)
%%
% for mi = 4:6
    mouse = mice(mi);
    load(sprintf('%sJK%03d\\behavior_JK%03d',bDir,mouse,mouse));
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
%%
%     for si = 2:3
        session = sessions{mi}(si);
        fprintf('Processing JK%03d S%02d\n', mouse, session)
        currB = b{find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b))};
        wfa = Whisker.WhiskerFinal_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));

        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
%         answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        firstLickTime = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouchBeforeFirstLick{mi,si} = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)) < y)), wfa.trials(answerTrialInds), firstLickTime);
        
        
        correct{mi,si} = cellfun(@(x) x.trialCorrect, wB(answerTrialInds));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, wB(answerTrialInds));
    end
end

%%

cellfun(@(x) length(find(x>2)) / length(x), numTouchBeforeFirstLick)
