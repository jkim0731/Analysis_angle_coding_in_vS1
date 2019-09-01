%% To see if mice used touch evidence accumulation to solve the task.
%% Are differences in multiple touches correlated with correct rate?

%% In 3 different expert sessions - 2 angles, 7 angles, radial distance.
%% In each angle. Either touch onset or middle point (accumulation of accumulation).
%% (1) Using all touches, from trials >=5 touches before the answer lick
%% (2) First 3 touches from trials >= 3 touches before the answer lick
%% (3) First 5 touches from trials >= 5 touches before the answer lick
%% (4) First 7 touches from trials >= 7 touches before the answer lick
%% All from touch chunks
%% 

%% (1) using all touches, from trials >=5 touches before the answer lick
%% Onset contact point

numTouchThreshold = 5;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
%         answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        answerLickTime = cellfun(@(x) x.beamBreakTimes(find(x.beamBreakTimes > x.poleUpOnsetTime,1)), wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        firstTouchIndsW3 = cellfun(@(x,y) cellfun(@(z) z(1), x(y)), touchIndsW3, touchIndsW3Inds, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, firstTouchIndsW3) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) x.intersectPoint(y,[1,3]), testW3(finalInd), firstTouchIndsW3(finalInd), 'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactDists{mi,si} = cellfun(@(x) sqrt(sum((x - x(1,:)).^2,2)) .* sign(sum(x(:,2)-x(1,2),2)) / pixPerMm, poleContactPoints{mi,si}, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end
%% Correlation between # of touches and spread of contact points
testAngles = [45,135];
groupWiseCorr = zeros([size(objectAngle),length(testAngles)]);
groupWiseP = zeros([size(objectAngle),length(testAngles)]);
totalCorr = zeros(1,size(objectAngle,2),length(testAngles));
totalP = zeros(1,size(objectAngle,2),length(testAngles));
for ai = 1 : 2
    for j = 1 : size(objectAngle,2)
        tempNum = [];
        tempStd = [];
        for i = 1 : size(objectAngle,1)
            angleInd = find(objectAngle{i,j} == testAngles(ai));
            [groupWiseCorr(i,j,ai), groupWiseP(i,j,ai)] = corr(cellfun(@length, poleContactDists{i,j}(angleInd))', contactSpread{i,j}(angleInd)');
            tempNum = [tempNum; cellfun(@length, poleContactDists{i,j}(angleInd))'];
            tempStd = [tempStd; contactSpread{i,j}(angleInd)'];
        end
        [totalCorr(1,j,ai), totalP(1,j,ai)] = corr(tempNum, tempStd);
    end
end
%%
offset = 0.2;
figure, hold on
for i = 1 : 3
    for j = 1 : 6
        if groupWiseP(j,i,1) < 0.05
            plot(i-offset, groupWiseCorr(j,i,1), 'r.', 'markersize', 10)
        else
            plot(i-offset, groupWiseCorr(j,i,1), 'k.', 'markersize', 10)
        end
        if groupWiseP(j,i,2) < 0.05
            plot(i+offset, groupWiseCorr(j,i,2), 'r.', 'markersize', 10)
        else
            plot(i+offset, groupWiseCorr(j,i,2), 'k.', 'markersize', 10)
        end
    end
end
%% Distribution of spread of contact points 

tempDist = cell2mat(cellfun(@(x) x', contactSpread(:,1), 'uniformoutput', false));
max(tempDist)
tempDist = cell2mat(cellfun(@(x) x', contactSpread(:,2), 'uniformoutput', false));
max(tempDist)
tempDist = cell2mat(cellfun(@(x) x', contactSpread(:,3), 'uniformoutput', false));
max(tempDist)

%% max dist = 2.8

%% 2 angles

histRange = 0:0.1:2.8;
% histRange = [0:0.05:0.8,3];

testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,1} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,1}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,1}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.3 0.7]),
subplot(211), hold on
for mi = 1 : 6
    plot(histRange(1:end-1), eachDist(mi,:,1), 'b-');
    plot(histRange(1:end-1), eachDist(mi,:,2), 'r-');
end
xlabel('Std touch points (mm)')
ylabel('Proportion')
legend({'45\circ', '135\circ'})
title('2 angles')
subplot(212), hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'b')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'r')
% plot(histRange(1:end-1), totalDist(1,:,1), 'b-', 'linewidth', 3)
% plot(histRange(1:end-1), totalDist(1,:,2), 'r-', 'linewidth', 3)
xlabel('Std touch points (mm)')
ylabel('Proportion')
% legend({'45\circ group', '135\circ group','45\circ total', '135\circ total'})
legend({'45\circ', '135\circ'})
%% 7 angles

% histRange = 0:0.1:2.8;
histRange = [0:0.05:0.8,3];

testAngles = 45:15:135;
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 7
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,2} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,2}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,2}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.3 0.7]), 
subplot(211), hold on
colorList = jet(7);
for mi = 1 : 6
    for ai = 1 : 7
        plot(histRange(1:end-1), eachDist(mi,:,ai), '-', 'color', colorList(ai,:));
    end
end
xlabel('Std touch points (mm)')
ylabel('Proportion')
legend({'45\circ', '60\circ', '75\circ', '90\circ', '105\circ', '120\circ', '135\circ'})
title('7 angles')
subplot(212), hold on
for ai = 1 : 7
    shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,ai)), std(eachDist(:,:,ai))/sqrt(6), 'lineprops', {'color',colorList(ai,:)})
%     plot(histRange(1:end-1), totalDist(1,:,ai), '--', 'linewidth', 3, 'color', colorList(ai,:))
end
xlabel('Std touch points (mm)')
ylabel('Proportion')


%% Radial distance

% histRange = 0:0.1:2.8;
histRange = [0:0.05:0.8,3];

testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,3} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,3}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,3}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.3 0.7]),
subplot(211), hold on
for mi = 1 : 6
    plot(histRange(1:end-1), eachDist(mi,:,1), 'c-');
    plot(histRange(1:end-1), eachDist(mi,:,2), 'm-');
end
xlabel('Std touch points (mm)')
ylabel('Proportion')
legend({'45\circ', '135\circ'})
title('Radial distance')
subplot(212), hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'c')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'm')
% plot(histRange(1:end-1), totalDist(1,:,1), 'b-', 'linewidth', 3)
% plot(histRange(1:end-1), totalDist(1,:,2), 'r-', 'linewidth', 3)
xlabel('Std touch points (mm)')
ylabel('Proportion')
legend({'45\circ', '135\circ'})

%% correlation between touch spread and correct rate

%% 2 angles
titleName = '2 angles';
groupInd = 1;

spreadRange = [0:0.3:1.5,3];

tempAnswer = cell2mat(cellfun(@(x) x', correct(:,groupInd), 'uniformoutput', false));
tempAngle = cell2mat(cellfun(@(x) x', objectAngle(:,groupInd), 'uniformoutput', false));
tempSpread = cell2mat(cellfun(@(x) x', contactSpread(:,groupInd), 'uniformoutput', false));

testAngles = unique(tempAngle);
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(testAngles), length(spreadRange)-1);

for ai = 1 : length(testAngles)
    angleInd = find(tempAngle == testAngles(ai));
    for ri = 1 : length(spreadRange)-1
        rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
        correctRate(ai,ri) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), correctRate(ai,:), 'color', colorList(ai,:));
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles
titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.3:1.5,3];

tempAnswer = cell2mat(cellfun(@(x) x', correct(:,groupInd), 'uniformoutput', false));
tempAngle = cell2mat(cellfun(@(x) x', objectAngle(:,groupInd), 'uniformoutput', false));
tempSpread = cell2mat(cellfun(@(x) x', contactSpread(:,groupInd), 'uniformoutput', false));

testAngles = unique(tempAngle);
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(testAngles), length(spreadRange)-1);

for ai = 1 : length(testAngles)
    angleInd = find(tempAngle == testAngles(ai));
    for ri = 1 : length(spreadRange)-1
        rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
        correctRate(ai,ri) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), correctRate(ai,:), 'color', colorList(ai,:));
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% radial distance
titleName = 'Radial distance';
groupInd = 3;

spreadRange = [0:0.3:1.5,3];

tempAnswer = cell2mat(cellfun(@(x) x', correct(:,groupInd), 'uniformoutput', false));
tempAngle = cell2mat(cellfun(@(x) x', objectAngle(:,groupInd), 'uniformoutput', false));
tempSpread = cell2mat(cellfun(@(x) x', contactSpread(:,groupInd), 'uniformoutput', false));

testAngles = unique(tempAngle);
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(testAngles), length(spreadRange)-1);

for ai = 1 : length(testAngles)
    angleInd = find(tempAngle == testAngles(ai));
    for ri = 1 : length(spreadRange)-1
        rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
        correctRate(ai,ri) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), correctRate(ai,:), 'color', colorList(ai,:));
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% Group-wise analysis

%% 2 angles
titleName = '2 angles';
groupInd = 1;

% spreadRange = [0:0.3:1.5,3];
spreadRange = [0:0.1:0.6,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles

titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.3:1.5,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), 'color',colorList(ai,:));
end

title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles - 45 and 135 only

titleName = '7 angles';
groupInd = 2;

% spreadRange = [0:0.3:1.5,3];
spreadRange = [0:0.05:0.4,3];

testAngles = [60,120];
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([2,6],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% Radial distance
titleName = 'Radial distance';
groupInd = 3;

% spreadRange = [0:0.3:1.5,3];
spreadRange = [0:0.1:0.6,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)


%% first 3 touches from trials >= 3 touches

numTouchThreshold = 3;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        firstTouchIndsW3 = cellfun(@(x,y) cellfun(@(z) z(1), x(y)), touchIndsW3, touchIndsW3Inds, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, firstTouchIndsW3) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) x.intersectPoint(y(1:numTouchThreshold),[1,3]), testW3(finalInd), firstTouchIndsW3(finalInd), 'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactDists{mi,si} = cellfun(@(x) sqrt(sum((x - x(1,:)).^2,2)) .* sign(sum(x(:,2)-x(1,2),2)) / pixPerMm, poleContactPoints{mi,si}, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end

%% first 5 touches from trials >= 5 touches

numTouchThreshold = 5;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        firstTouchIndsW3 = cellfun(@(x,y) cellfun(@(z) z(1), x(y)), touchIndsW3, touchIndsW3Inds, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, firstTouchIndsW3) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) x.intersectPoint(y(1:numTouchThreshold),[1,3]), testW3(finalInd), firstTouchIndsW3(finalInd), 'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactDists{mi,si} = cellfun(@(x) sqrt(sum((x - x(1,:)).^2,2)) .* sign(sum(x(:,2)-x(1,2),2)) / pixPerMm, poleContactPoints{mi,si}, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end


%% first 7 touches from trials >= 7 touches

numTouchThreshold = 7;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        firstTouchIndsW3 = cellfun(@(x,y) cellfun(@(z) z(1), x(y)), touchIndsW3, touchIndsW3Inds, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, firstTouchIndsW3) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) x.intersectPoint(y(1:numTouchThreshold),[1,3]), testW3(finalInd), firstTouchIndsW3(finalInd), 'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactDists{mi,si} = cellfun(@(x) sqrt(sum((x - x(1,:)).^2,2)) .* sign(sum(x(:,2)-x(1,2),2)) / pixPerMm, poleContactPoints{mi,si}, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end













%% touch mid point, not the onset point







%% from trials >= 3 touches

numTouchThreshold = 3;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, touchIndsW3Inds) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,[1,3]), y, 'uniformoutput', false), testW3(finalInd), touchIndsW3(finalInd), 'uniformoutput', false);
        poleContactInds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), poleContactPoints{mi,si}, 'uniformoutput', false);
        poleContactOrigin = cellfun(@(x,y) x{y(1)}(1,:), poleContactPoints{mi,si}, poleContactInds,'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactPointsLinear = cellfun(@(x,y,w) cellfun(@(z) sqrt(sum((z-y).^2,2)) .* sign(sum(z(:,2)-y(2),2)) / pixPerMm, x(w),'uniformoutput', false), poleContactPoints{mi,si}, poleContactOrigin, poleContactInds, 'uniformoutput', false);
        %%
        poleContactDists{mi,si} = cellfun(@(x) cellfun(@(y) (max(y)-min(y))/2, x), poleContactPointsLinear, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end



%% from trials >= 7 touches

numTouchThreshold = 7;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, touchIndsW3Inds) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,[1,3]), y, 'uniformoutput', false), testW3(finalInd), touchIndsW3(finalInd), 'uniformoutput', false);
        poleContactInds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), poleContactPoints{mi,si}, 'uniformoutput', false);
        poleContactOrigin = cellfun(@(x,y) x{y(1)}(1,:), poleContactPoints{mi,si}, poleContactInds,'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactPointsLinear = cellfun(@(x,y,w) cellfun(@(z) sqrt(sum((z-y).^2,2)) .* sign(sum(z(:,2)-y(2),2)) / pixPerMm, x(w),'uniformoutput', false), poleContactPoints{mi,si}, poleContactOrigin, poleContactInds, 'uniformoutput', false);
        %%
        poleContactDists{mi,si} = cellfun(@(x) cellfun(@(y) (max(y)-min(y))/2, x), poleContactPointsLinear, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end

































%% Summary figure
%% Onset contact point

numTouchThreshold = 5;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        firstTouchIndsW3 = cellfun(@(x,y) cellfun(@(z) z(1), x(y)), touchIndsW3, touchIndsW3Inds, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, firstTouchIndsW3) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) x.intersectPoint(y,[1,3]), testW3(finalInd), firstTouchIndsW3(finalInd), 'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactDists{mi,si} = cellfun(@(x) sqrt(sum((x - x(1,:)).^2,2)) .* sign(sum(x(:,2)-x(1,2),2)) / pixPerMm, poleContactPoints{mi,si}, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end

%%
histRange = 0:0.1:2.8;

testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,1} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,1}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,1}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure, hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'b')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'r')
xlabel('Std touch points (mm)'), ylabel('Proportion'), legend({'45\circ', '135\circ'}), title('2 angles')

%%
testAngles = 45:15:135;
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 7
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,2} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,2}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,2}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end


figure, hold on
colorList = jet(7);
for ai = 1 : 7
    shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,ai)), std(eachDist(:,:,ai))/sqrt(6), 'lineprops', {'color',colorList(ai,:)})
end
for ai = 1 : 7
    plot(histRange(1:end-1), mean(eachDist(:,:,ai)), 'color',colorList(ai,:))
end
xlabel('Std touch points (mm)'), ylabel('Proportion'), title('7 angles')
legend({'45\circ', '60\circ', '75\circ', '90\circ', '105\circ', '120\circ', '135\circ'})
%%
histRange = 0:0.1:2.8;

testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,3} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,3}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,3}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure, hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'b')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'r')
plot(histRange(1:end-1), mean(eachDist(:,:,1)), 'b-')
plot(histRange(1:end-1), mean(eachDist(:,:,2)), 'r-')
xlabel('Std touch points (mm)'), ylabel('Proportion'), legend({'45\circ', '135\circ'}), title('Radial distance')

%% 2 angles
titleName = '2 angles';
groupInd = 1;

spreadRange = [0:0.2:1.4,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), 'color',colorList(ai,:));
end

title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles

titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.2:1.4,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), 'color',colorList(ai,:));
end

title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles - 45 and 135 only

titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.2:1.4,3];
% spreadRange = [0:0.05:0.4,3];

testAngles = [45,135];
tempColors = jet(7);
if length(testAngles) == 2
    clrinds = find(ismember(45:15:135, testAngles));
    colorList = tempColors(clrinds,:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% Radial distance
titleName = 'Radial distance';
groupInd = 3;

spreadRange = [0:0.2:1.4,3];
% spreadRange = [0:0.1:0.6,3];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)










%% touch mid point
numTouchThreshold = 5;

wDir = 'D:\TPM\JK\tracked\';
bDir = 'D:\TPM\JK\soloData\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
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
        w3a = Whisker.Whisker3D_2padArray(sprintf('%sJK%03dS%02d',wDir,mouse,session));
        wB = currB.trials(find(ismember(currB.trialNums, wfa.trialNums)));
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), wB));        
        answerLickTime = cellfun(@(x) x.answerLickTime, wB(answerTrialInds), 'uniformoutput', false); % no empty cell because these are from answer trials
        numTouch = cellfun(@(x,y) length(find(x.time(cellfun(@(z) z(1), x.protractionTFchunks)) < y)), wfa.trials(answerTrialInds), answerLickTime, 'uniformoutput', false);
        
        testInds = find(cell2mat(numTouch) >= numTouchThreshold);
        testB = wB(answerTrialInds(testInds));
        testWF = wfa.trials(answerTrialInds(testInds));
        testW3 = w3a.trials(answerTrialInds(testInds));
        
        touchTimes = cellfun(@(x,z) cellfun(@(y) x.time(y), x.protractionTFchunks(1:z), 'uniformoutput', false), testWF, numTouch(testInds), 'uniformoutput', false);
        touchIndsW3 = cellfun(@(x,y) cellfun(@(z) find(ismember(x.time,z)), y, 'uniformoutput', false), testW3, touchTimes, 'uniformoutput', false);
        touchIndsW3Inds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), touchIndsW3, 'uniformoutput', false);
        
        finalInd = find(cellfun(@length, touchIndsW3Inds) >= numTouchThreshold);
        %%
        poleContactPoints{mi,si} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,[1,3]), y, 'uniformoutput', false), testW3(finalInd), touchIndsW3(finalInd), 'uniformoutput', false);
        poleContactInds = cellfun(@(x) find(cellfun(@(y) ~isempty(y), x)), poleContactPoints{mi,si}, 'uniformoutput', false);
        poleContactOrigin = cellfun(@(x,y) x{y(1)}(1,:), poleContactPoints{mi,si}, poleContactInds,'uniformoutput', false);
        pixPerMm = testW3{1}.pxPerMm;
        poleContactPointsLinear = cellfun(@(x,y,w) cellfun(@(z) sqrt(sum((z-y).^2,2)) .* sign(sum(z(:,2)-y(2),2)) / pixPerMm, x(w),'uniformoutput', false), poleContactPoints{mi,si}, poleContactOrigin, poleContactInds, 'uniformoutput', false);
        %%
        poleContactDists{mi,si} = cellfun(@(x) cellfun(@(y) (max(y)-min(y))/2, x), poleContactPointsLinear, 'uniformoutput', false);
        contactSpread{mi,si} = cellfun(@std, poleContactDists{mi,si});
        correct{mi,si} = cellfun(@(x) x.trialCorrect, testB(finalInd));
        objectAngle{mi,si} = cellfun(@(x) x.servoAngle, testB(finalInd));
    end
end


%%
histRange = [0:0.05:1,2];

testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,1} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,1}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,1}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure, hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'b')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'r')
xlabel('Std touch points (mm)'), ylabel('Proportion'), legend({'45\circ', '135\circ'}), title('2 angles')
yy = ylim;
ylim([0, yy(2)])
%%
testAngles = 45:15:135;
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 7
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,2} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,2}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,2}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end


figure, hold on
colorList = jet(7);
for ai = 1 : 7
    shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,ai)), std(eachDist(:,:,ai))/sqrt(6), 'lineprops', {'color',colorList(ai,:)})
end
for ai = 1 : 7
    plot(histRange(1:end-1), mean(eachDist(:,:,ai)), 'color',colorList(ai,:))
end
xlabel('Std touch points (mm)'), ylabel('Proportion'), title('7 angles')
legend({'45\circ', '60\circ', '75\circ', '90\circ', '105\circ', '120\circ', '135\circ'})
yy = ylim;
ylim([0, yy(2)])
%%
testAngles = [45,135];
eachDist = zeros(6,length(histRange)-1, length(testAngles));
totalDist = zeros(1,length(histRange)-1, length(testAngles));
for ai = 1 : 2
    tempDist = [];
    for mi = 1 : 6
        angleInd = find(objectAngle{mi,3} == testAngles(ai));
        eachDist(mi,:,ai) = histcounts(contactSpread{mi,3}(angleInd),histRange,'normalization','probability');
        tempDist = [tempDist, contactSpread{mi,3}(angleInd)];
    end
    totalDist(1,:,ai) = histcounts(tempDist,histRange,'normalization','probability');
end

figure, hold on
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,1)), std(eachDist(:,:,1))/sqrt(6), 'lineprops', 'b')
shadedErrorBar(histRange(1:end-1), mean(eachDist(:,:,2)), std(eachDist(:,:,2))/sqrt(6), 'lineprops', 'r')
plot(histRange(1:end-1), mean(eachDist(:,:,1)), 'b-')
plot(histRange(1:end-1), mean(eachDist(:,:,2)), 'r-')
xlabel('Std touch points (mm)'), ylabel('Proportion'), legend({'45\circ', '135\circ'}), title('Radial distance')
yy = ylim;
ylim([0, yy(2)])


%% 2 angles
titleName = '2 angles';
groupInd = 1;

spreadRange = [0:0.1:0.6,2];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), 'color',colorList(ai,:));
end

title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles

titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.1:0.6,2];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
for ai = 1 : length(testAngles)
    plot(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), 'color',colorList(ai,:));
end

title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% 7 angles - 45 and 135 only

titleName = '7 angles';
groupInd = 2;

spreadRange = [0:0.05:0.5,3];

testAngles = [75,105];
tempColors = jet(7);
if length(testAngles) == 2
    clrinds = find(ismember(45:15:135, testAngles));
    colorList = tempColors(clrinds,:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)

%% Radial distance
titleName = 'Radial distance';
groupInd = 3;

spreadRange = [0:0.1:0.6,2];

testAngles = unique(objectAngle{1,groupInd});
tempColors = jet(7);
if length(testAngles) == 2
    colorList = tempColors([1,7],:);
else
    colorList = tempColors;
end

legendList = cell(1,length(testAngles));
for ai = 1 : length(testAngles)
    legendList{ai} = [num2str(testAngles(ai)), '\circ'];
end

correctRate = zeros(length(mice), length(spreadRange)-1, length(testAngles));

for mi = 1 : length(mice)
    tempAngle = objectAngle{mi,groupInd};
    tempSpread = contactSpread{mi,groupInd};
    tempAnswer = correct{mi,groupInd};
    for ai = 1 : length(testAngles)
        angleInd = find(tempAngle == testAngles(ai));
        for ri = 1 : length(spreadRange)-1
            rangeInd = find(tempSpread(angleInd) >= spreadRange(ri) & tempSpread(angleInd) < spreadRange(ri+1));
            correctRate(mi,ri,ai) = sum(tempAnswer(angleInd(rangeInd))) / length(rangeInd);
        end
    end
end

figure, hold on
for ai = 1 : length(testAngles)
    shadedErrorBar(spreadRange(1:end-1), nanmean(correctRate(:,:,ai)), nanstd(correctRate(:,:,ai)), 'lineprops', {'color',colorList(ai,:)});
end
title(titleName)
xlabel('Std of touch points (mm)')
ylabel('Correct rate')
legend(legendList)






%% Next question: what does correlate with correct rate? 
%% (1) Mean touch duration, (2) mean slide distance, (3) mean dKv, (4) mean dPhi
%% from touch by whisking
uDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,39,52];
sessions = {[18,19,22],[7,10,14],[20,21,22],[16,17,18],[21,23,24],[20,21,26]}; 

correct = cell(length(mice),3);
poleContactPoints = cell(length(mice),3); % 2d value, considering only x and z dimension (horizontal and vertical)
poleContactDists = cell(length(mice),3); % 1d value, distance from the first contact point. in mm.
contactSpread = cell(length(mice),3); % std of pole contact distances
objectAngle = cell(length(mice),3);
for mi = 1 : length(mice)
    mouse = mice(mi);    
    for si = 1 : 3 % 1 for 2 angles expert, 2 for 7 angles expert, 3 for radial distance
        session = sessions{mi}(si);
        fprintf('Processing JK%03d S%02d\n', mouse, session)
        load(sprintf('%s%03d\\UberJK%03dS%02d_NC',uDir,mouse,mouse,session))
        %%
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), u.trials));        
        touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));        
        testInds = intersect(answerTrialInds, touchTrialInds);
        
        preTouchInds = cellfun(@(x) find(cellfun(@(y) (x.whiskerTime(y(1)) < x.answerLickTime), x.protractionTouchChunksByWhisking)), u.trials(testInds), 'uniformoutput', false);
        
%         duration = 
        
        correct{mi,si} = cellfun(@(x) x.response, u.trials(testInds));
        objectAngle{mi,si} = cellfun(@(x) x.angle, u.trials(testInds));
    end
end
