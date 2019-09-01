%% See if there is motor strategy change before and after learning
% (1) touch
% # of touches
% change in whisker variables during touch
% (2) whisking
% # of whisking
% amplitude, midpoint, theta
% (3) response
% response time


%% (1) touch - 1) # of touches
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
            numTouch(ei,ai,si) = mean(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunksByWhisking)), u.trials(angleTouchTrialInds)));
            angleAnswerTrialInds = intersect(angleInds, answerTrialInds);
            meanTouchNum(ei,ai,si) = mean(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunksByWhisking)), u.trials(angleAnswerTrialInds)));
            propTouchTrial(ei,ai,si) = length(find(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunksByWhisking)), u.trials(angleAnswerTrialInds)))) / length(angleAnswerTrialInds);
        end        
    end
end

%

figure, 
subplot(231), hold on
shadedErrorBar(angles, mean(numTouch(:,:,1)), std(numTouch(:,:,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(numTouch(:,:,2)), std(numTouch(:,:,2)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(from touch trials only)'})

subplot(232), hold on
shadedErrorBar(angles, mean(meanTouchNum(:,:,1)), std(meanTouchNum(:,:,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(meanTouchNum(:,:,2)), std(meanTouchNum(:,:,2)), 'lineprop', 'r-')
xticks(angles), title({'Mean number of touches'; '(including non-touch trials)'})

subplot(233), hold on
shadedErrorBar(angles, mean(propTouchTrial(:,:,1)), std(propTouchTrial(:,:,1)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(propTouchTrial(:,:,2)), std(propTouchTrial(:,:,2)), 'lineprop', 'r-')
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


%% (1) touch - 2) change in whisker variables during touch
%% touch chunks by whisking

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
dTheta = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Only in answer and touch trials. Mean of max of each trial.
dKappaH = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. 
dPhi = zeros(length(expertInds),length(angles),2); % Pre-decision. 
dKappaV = zeros(length(expertInds),length(angles),2); % Pre-decision. 
touchDuration = zeros(length(expertInds),length(angles),2); % Pre-decision. 
slideDistance = zeros(length(expertInds),length(angles),2); % Pre-decision. 
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
            numPreDecisionTouch = cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunksByWhisking)), u.trials(angleTouchTrialInds), 'uniformoutput', false);
            preDecisionTouchTrialInds = angleTouchTrialInds(find(cell2mat(numPreDecisionTouch)));
            dTheta(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchDThetaByWhisking(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dKappaH(ei,ai,si) = mean(cellfun(@(x,y) min(x.protractionTouchDKappaHByWhisking(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dPhi(ei,ai,si) = mean(cellfun(@(x,y) x.protractionTouchDPhiByWhisking(find(abs(x.protractionTouchDPhiByWhisking(1:y)) == max(abs(x.protractionTouchDPhiByWhisking(1:y))),1)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dKappaV(ei,ai,si) = mean(cellfun(@(x,y) x.protractionTouchDKappaVByWhisking(find(abs(x.protractionTouchDKappaVByWhisking(1:y)) == max(abs(x.protractionTouchDKappaVByWhisking(1:y))),1)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            touchDuration(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchDurationByWhisking(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            slideDistance(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchSlideDistanceByWhisking(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
        end        
    end
end

figure, 
for i = 1 : 6
    subplot(2,3,i), hold on
    switch i
        case 1
            temp = dTheta;
            titleName = 'Max \Delta\theta';
        case 2
            temp = dKappaH;
            titleName = 'Max \Delta\kappa_H';
        case 3
            temp = dPhi;
            titleName = 'Max \Delta\phi';
        case 4
            temp = dKappaV;
            titleName = 'Max \Delta\kappa_V';
        case 5
            temp = touchDuration;
            titleName = 'Max touch duration';
        case 6
            temp = slideDistance;
            titleName = 'Max slide distance';
    end
    
    shadedErrorBar(angles, mean(temp(:,:,1)), std(temp(:,:,1)), 'lineprop', 'b-')
    shadedErrorBar(angles, mean(temp(:,:,2)), std(temp(:,:,2)), 'lineprop', 'r-')
    xticks(angles), title(titleName)
    if i > 3
        xlabel('Angle (\circ)')
    end
end

%% touch chunks itself

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
dTheta = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. Only in answer and touch trials. Mean of max of each trial.
dKappaH = zeros(length(expertInds),length(angles),2); % for both naive and expert. Predecision. 
dPhi = zeros(length(expertInds),length(angles),2); % Pre-decision. 
dKappaV = zeros(length(expertInds),length(angles),2); % Pre-decision. 
touchDuration = zeros(length(expertInds),length(angles),2); % Pre-decision. 
slideDistance = zeros(length(expertInds),length(angles),2); % Pre-decision. 
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
            numPreDecisionTouch = cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunks)), u.trials(angleTouchTrialInds), 'uniformoutput', false);
            preDecisionTouchTrialInds = angleTouchTrialInds(find(cell2mat(numPreDecisionTouch)));
            dTheta(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchDTheta(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dKappaH(ei,ai,si) = mean(cellfun(@(x,y) min(x.protractionTouchDKappaH(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dPhi(ei,ai,si) = mean(cellfun(@(x,y) x.protractionTouchDPhi(find(abs(x.protractionTouchDPhi(1:y)) == max(abs(x.protractionTouchDPhi(1:y))),1)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            dKappaV(ei,ai,si) = mean(cellfun(@(x,y) x.protractionTouchDKappaV(find(abs(x.protractionTouchDKappaV(1:y)) == max(abs(x.protractionTouchDKappaV(1:y))),1)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            touchDuration(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchDuration(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
            slideDistance(ei,ai,si) = mean(cellfun(@(x,y) max(x.protractionTouchSlideDistance(1:y)), u.trials(preDecisionTouchTrialInds), numPreDecisionTouch(find(cell2mat(numPreDecisionTouch)))));
        end        
    end
end
%%
figure, 
for i = 1 : 6
    subplot(2,3,i), hold on
    switch i
        case 1
            temp = dTheta;
            titleName = 'Max \Delta\theta';
        case 2
            temp = dKappaH;
            titleName = 'Max \Delta\kappa_H';
        case 3
            temp = dPhi;
            titleName = 'Max \Delta\phi';
        case 4
            temp = dKappaV;
            titleName = 'Max \Delta\kappa_V';
        case 5
            temp = touchDuration;
            titleName = 'Max touch duration';
        case 6
            temp = slideDistance;
            titleName = 'Max slide distance';
    end
    
    shadedErrorBar(angles, mean(temp(:,:,1)), std(temp(:,:,1)), 'lineprop', 'b-')
    shadedErrorBar(angles, mean(temp(:,:,2)), std(temp(:,:,2)), 'lineprop', 'r-')
    xticks(angles), title(titleName)
    if i > 3
        xlabel('Angle (\circ)')
    end
end

%% (2) whisking - 1) # of whisking

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
numWhisking = zeros(length(expertInds),length(angles),2); % for both naive and expert. From pole onset to pre-decision. 
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
        catchTrialInds = find(cellfun(@(x) isempty(x.poleUpTime), u.trials));
        for ai = 1:length(angles)            
            angleTrialInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            trialInds = setdiff(intersect(angleTrialInds, answerTrialInds), catchTrialInds);
            testFrames = cellfun(@(x) find(x.whiskerTime > x.poleUpTime(1),1,'first'):find(x.whiskerTime < x.answerLickTime, 1, 'last'), u.trials(trialInds), 'uniformoutput', false);
            numWhisking(ei,ai,si) = mean(cellfun(@(x,y) length(jkWhiskerOnsetNAmplitude(x.theta(y),2.5)), u.trials(trialInds), testFrames));
        end
    end
end

figure,
subplot(121), hold on
shadedErrorBar(angles, mean(squeeze(numWhisking(:,:,1))), std(squeeze(numWhisking(:,:,1))), 'lineprop', 'b-')
shadedErrorBar(angles, mean(squeeze(numWhisking(:,:,2))), std(squeeze(numWhisking(:,:,2))), 'lineprop', 'r-')
xticks(angles), title('Number of whiskings'), xlim([angles(1) angles(end)]), xlabel('Angle (\circ)')
legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')
subplot(122),
shadedErrorBar(angles, mean(squeeze(numWhisking(:,:,2) - numWhisking(:,:,1))), std(squeeze(numWhisking(:,:,2) - numWhisking(:,:,1))), 'lineprop', 'k-')
xticks(angles), title('Number of whiskings'), xlim([angles(1) angles(end)]), xlabel('Angle (\circ)')
legend('Expert - Naive', 'location', 'northeast', 'box', 'off')


%% (2) whisking - 2) amplitude, midpoint, and theta
% angle-degrees image
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
ampDegrees = 0:1:10;
midDegrees = -30:1:20;
thetaDegrees = -30:1:20;

amplitude = zeros(length(expertInds),length(angles),length(ampDegrees)-1,2); % for both naive and expert. From pole onset to pre-decision. pdf 
midpoint = zeros(length(expertInds),length(angles),length(midDegrees)-1,2); 
theta = zeros(length(expertInds),length(angles),length(thetaDegrees)-1,2);


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
        catchTrialInds = find(cellfun(@(x) isempty(x.poleUpTime), u.trials));        
        for ai = 1 : length(angles)
            angleTrialInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            trialInds = setdiff(intersect(angleTrialInds, answerTrialInds), catchTrialInds);
            testFrames = cellfun(@(x) find(x.whiskerTime > x.poleUpTime(1),1,'first'):find(x.whiskerTime < x.answerLickTime, 1, 'last'), u.trials(trialInds), 'uniformoutput', false);
            
            tempTheta = zeros(length(testFrames),length(thetaDegrees)-1);
            tempAmp = zeros(length(testFrames),length(ampDegrees)-1);
            tempMid = zeros(length(testFrames),length(midDegrees)-1);
            

            for ti = 1 : length(trialInds)
                th = u.trials{trialInds(ti)}.theta(testFrames{ti});
                [~,amp,mid] = jkWhiskerOnsetNAmplitude(th);
                
                tempTheta(ti,:) = histcounts(th, thetaDegrees, 'normalization','pdf');
                tempAmp(ti,:) = histcounts(amp, ampDegrees, 'normalization','pdf');
                tempMid(ti,:) = histcounts(mid, midDegrees, 'normalization','pdf');
%                 tempTheta(ti,:) = histcounts(th, thetaDegrees, 'normalization','cdf');
%                 tempAmp(ti,:) = histcounts(amp, ampDegrees, 'normalization','cdf');
%                 tempMid(ti,:) = histcounts(mid, midDegrees, 'normalization','cdf');
            end            
            
            theta(ei,ai,:,si) = mean(tempTheta);
            amplitude(ei,ai,:,si) = mean(tempAmp);
            midpoint(ei,ai,:,si) = mean(tempMid);
        end
    end
end

%%
figure,
for i = 1 : 3
    switch i
        case 1
            temp = amplitude;
            titleName = 'Amplitude';
        case 2
            temp = midpoint;
            titleName = 'Midpoint';
        case 3
            temp = theta;
            titleName = 'Theta';
    end
    temp1 = squeeze(mean(squeeze(temp(:,:,:,1)),1));
    temp2 = squeeze(mean(squeeze(temp(:,:,:,2)),1));
    clim = [min([min(min(temp1)), min(min(temp2))]), max([max(max(temp1)), max(max(temp2))])];
    subplot(3,3,i), imagesc(temp1, clim), colorbar
    title([titleName, ' (Naive)']), ylabel('Object angles (\circ)'), yticks(1:length(angles)), yticklabels(angles)
    
    subplot(3,3,i+3), imagesc(temp2, clim), colorbar
    title([titleName, ' (Expert)']), ylabel('Object angles (\circ)'), yticks(1:length(angles)), yticklabels(angles)
    
    subplot(3,3,i+6), imagesc(temp2-temp1), colorbar    
    title([titleName, ' (Expert - Naive)']), ylabel('Object angles (\circ)'), yticks(1:length(angles)), yticklabels(angles)
    xlabel('Degrees (\circ)')
end

%%
figure,
for i = 1 : 3
    switch i
        case 1
            temp = amplitude;
            titleName = 'Amplitude';
            x = ampDegrees;
        case 2
            temp = midpoint;
            titleName = 'Midpoint';
            x = midDegrees;
        case 3
            temp = theta;
            titleName = 'Theta';
            x = thetaDegrees;
    end
    temp1 = squeeze(mean(squeeze(temp(:,:,:,1)),2));
    temp2 = squeeze(mean(squeeze(temp(:,:,:,2)),2));
    
    subplot(2,3,i), 
    shadedErrorBar(x(2:end), mean(temp1), std(temp1), 'lineprop', 'b-')
    shadedErrorBar(x(2:end), mean(temp2), std(temp2), 'lineprop', 'r-')
    title(titleName), xlabel('Degrees (\circ)')
    if i == 1
         ylabel('Proportion')
    end
    if i == 3
        legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')
    end
    subplot(2,3,i+3)
    shadedErrorBar(x(2:end), mean(temp1-temp2), std(temp1-temp2), 'lineprop', 'k-')
    xlabel('Degrees (\circ)')
    if i == 1
         ylabel('Proportion')
    end
    if i == 3
        legend({'Expert - Naive'}, 'location', 'northeast', 'box', 'off')
    end    
end



%% (2) whisking - 2) amplitude, midpoint, and theta
% from all angles
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));

ampDegrees = 0:0.1:10;
midDegrees = -40:1:20;
thetaDegrees = -40:1:20;

meanAmp = zeros(length(expertInds),2);
meanMid = zeros(length(expertInds),2);
meanTheta = zeros(length(expertInds),2);

amplitude = zeros(length(expertInds),length(ampDegrees)-1,2); % for both naive and expert. From pole onset to pre-decision. pdf 
midpoint = zeros(length(expertInds),length(midDegrees)-1,2);
theta = zeros(length(expertInds),length(thetaDegrees)-1,2);

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
        catchTrialInds = find(cellfun(@(x) isempty(x.poleUpTime), u.trials));
        
        trialInds = setdiff(answerTrialInds, catchTrialInds);
        testFrames = cellfun(@(x) find(x.whiskerTime > x.poleUpTime(1),1,'first'):find(x.whiskerTime < x.answerLickTime, 1, 'last'), u.trials(trialInds), 'uniformoutput', false);

        tempMeanTheta = zeros(length(testFrames),1);
        tempMeanAmp = zeros(length(testFrames),1);
        tempMeanMid = zeros(length(testFrames),1);

        tempTheta = zeros(length(testFrames),length(thetaDegrees)-1);
        tempAmp = zeros(length(testFrames),length(ampDegrees)-1);
        tempMid = zeros(length(testFrames),length(midDegrees)-1);

        for ti = 1 : length(trialInds)
            th = u.trials{trialInds(ti)}.theta(testFrames{ti});
            [~,amp,mid] = jkWhiskerOnsetNAmplitude(th);

            tempMeanTheta(ti) = nanmean(th);
            tempMeanAmp(ti) = mean(amp);
            tempMeanMid(ti) = mean(mid);

            tempTheta(ti,:) = histcounts(th, thetaDegrees, 'normalization','pdf');
            tempAmp(ti,:) = histcounts(amp, ampDegrees, 'normalization','pdf');
            tempMid(ti,:) = histcounts(mid, midDegrees, 'normalization','pdf');
%                 tempTheta(ti,:) = histcounts(th, thetaDegrees, 'normalization','cdf');
%                 tempAmp(ti,:) = histcounts(amp, ampDegrees, 'normalization','cdf');
%                 tempMid(ti,:) = histcounts(mid, midDegrees, 'normalization','cdf');
        end
        
        meanAmp(ei,si) = mean(tempMeanAmp);
        meanMid(ei,si) = mean(tempMeanMid);
        meanTheta(ei,si) = mean(tempMeanTheta);

        theta(ei,:,si) = mean(tempTheta);
        amplitude(ei,:,si) = mean(tempAmp);
        midpoint(ei,:,si) = mean(tempMid);
    end
end


%%
figure,
for i = 1 : 3
    switch i
        case 1
            temp = amplitude;
            titleName = 'Amplitude';
            x = ampDegrees;
            tempMean = meanAmp;
            y = 0.3;
        case 2
            temp = midpoint;
            titleName = 'Midpoint';
            x = midDegrees;
            tempMean = meanMid;
            y = 0.09;
        case 3
            temp = theta;
            titleName = 'Theta';
            x = thetaDegrees;
            tempMean = meanTheta;
            y = 0.09;
    end
    temp1 = squeeze(temp(:,:,1));
    temp2 = squeeze(temp(:,:,2));
    
    subplot(2,3,i), hold on
    shadedErrorBar(x(2:end), mean(temp1), std(temp1), 'lineprop', 'b-')    
    shadedErrorBar(x(2:end), mean(temp2), std(temp2), 'lineprop', 'r-')
    mean1 = mean(tempMean(:,1));
    mean2 = mean(tempMean(:,2));
%     [~, x1ind] = min(abs(x(2:end)-mean1));
%     [~, x2ind] = min(abs(x(2:end)-mean2));
%     x1 = x(x1ind+1); x2 = x(x2ind+1);
%     y = max([mean(temp1(:,x1ind))+std(temp1(:,x1ind))*2, mean(temp2(:,x2ind))+std(temp2(:,x2ind))*2]);
    scatter(mean1, y, 'bv', 'filled')
    scatter(mean2, y, 'rv', 'filled')
    title(titleName), xlabel('Degrees (\circ)')
    xlim([x(1) x(end)])
    if i == 1
        ylabel('Proportion')
        legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')
    end
    subplot(2,3,i+3)
    shadedErrorBar(x(2:end), mean(temp2-temp1), std(temp2-temp1), 'lineprop', 'k-')
    xlim([x(1) x(end)])
    xlabel('Degrees (\circ)')
    if i == 1
        ylabel('Proportion')
        legend({'Expert - Naive'}, 'location', 'northeast', 'box', 'off')
    end
end

%%
mean(meanAmp)
[~, h] = ttest(meanAmp(:,1),meanAmp(:,2))

mean(meanMid)
[~, h] = ttest(meanMid(:,1),meanMid(:,2))

mean(meanTheta)
[~, h] = ttest(meanTheta(:,1),meanTheta(:,2))


%% (3) response
% response time - by answer time

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
answerLickTime = zeros(length(expertInds), length(angles), 2);

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
        for ai = 1 : length(angles)
            angleTrialInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            trialInds = intersect(angleTrialInds, answerTrialInds);
            answerLickTime(ei,ai,si) = mean(cellfun(@(x) x.answerLickTime - x.poleUpOnsetTime, u.trials(trialInds)));
        end
    end
end

figure, hold on
shadedErrorBar(angles, mean(answerLickTime(:,:,1)), std(answerLickTime(:,:,1))/sqrt(length(expertInds)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(answerLickTime(:,:,2)), std(answerLickTime(:,:,2))/sqrt(length(expertInds)), 'lineprop', 'r-')

xlabel('Angle (\circ)'), ylabel('Time after pole onset (s)')
xticks(angles)
title('Answer lick time')
legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')

%% (3) response
% response time - by the first lick time after pole up, only when it was
% correct (because of alternation before the answer)

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
firstLickTime = zeros(length(expertInds), length(angles), 2);

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
        for ai = 1 : length(angles)
            angleTrialInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            correctTrialInds = find(cellfun(@(x) x.response == 1, u.trials));
            trialInds = intersect(intersect(angleTrialInds, answerTrialInds), correctTrialInds);
            allLickTimes = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(trialInds), 'uniformoutput', false);
            firstLickInds = cellfun(@(x,y) find(y > x.poleUpOnsetTime, 1), u.trials(trialInds), allLickTimes, 'uniformoutput', false);
            firstLickTime(ei,ai,si) = mean(cellfun(@(x,y,z) y(z) - x.poleUpOnsetTime, u.trials(trialInds), allLickTimes, firstLickInds));
        end
    end
end
%
figure, hold on
shadedErrorBar(angles, mean(firstLickTime(:,:,1)), std(firstLickTime(:,:,1))/sqrt(length(expertInds)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(firstLickTime(:,:,2)), std(firstLickTime(:,:,2))/sqrt(length(expertInds)), 'lineprop', 'r-')

xlabel('Angle (\circ)'), ylabel('Time after pole onset (s)')
xticks(angles)
title('First lick time')
% legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')


%% (3) Response
% what is the average lick rate between first lick and answer lick?

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
angles = 45:15:135;
firstLickTime = zeros(length(expertInds), length(angles), 2);
correctAnswerLickTime = zeros(length(expertInds), length(angles), 2);
lickRate = zeros(length(expertInds), length(angles), 2);
wastedLicks = zeros(length(expertInds), length(angles), 2);
for ei = 1 : length(expertInds)
% for ei = 5
    mouse = mice(expertInds(ei));
    for si = 1 : 2 % 1 for naive, 2 for expert
%     for si = 1
        session = sessions{expertInds(ei)}(si);
        ufn = sprintf('UberJK%03dS%02d', mouse, session);
        dn = sprintf('%s%03d\\',baseDir,mouse);
        load([dn,ufn],'u')
        answerTrialInds = find(cellfun(@(x) length(x.answerLickTime), u.trials));
        for ai = 1 : length(angles)
%         for ai = 6
            angleTrialInds = find(cellfun(@(x) x.angle == angles(ai), u.trials));
            correctTrialInds = find(cellfun(@(x) x.response == 1, u.trials));
            trialInds = intersect(intersect(angleTrialInds, answerTrialInds), correctTrialInds);
            allLickTimes = cellfun(@(x) union(union(x.leftLickTime, x.rightLickTime), x.answerLickTime), u.trials(trialInds), 'uniformoutput', false);
            
            for alti = 1 : length(allLickTimes)
                temp = allLickTimes{alti};
                while min(diff(temp)) < 1/10
                    errorInd = find(diff(temp) < 1/10, 1);
                    temp(errorInd+1) = [];
                end
                allLickTimes{alti} = temp;
            end
            firstLickInds = cellfun(@(x,y) find(y > x.poleUpOnsetTime, 1), u.trials(trialInds), allLickTimes, 'uniformoutput', false);
            answerLickInds = cellfun(@(x,y) find(y > x.poleUpOnsetTime+1,1), u.trials(trialInds), allLickTimes, 'uniformoutput', false);
            emptyInd = find(cellfun(@isempty, answerLickInds));
            if ~isempty(emptyInd)
                answerLickInds(emptyInd) = firstLickInds(emptyInd);
            end
            firstLickTime(ei,ai,si) = mean(cellfun(@(x,y,z) y(z) - x.poleUpOnsetTime, u.trials(trialInds), allLickTimes, firstLickInds));
            correctAnswerLickTime(ei,ai,si) = mean(cellfun(@(x) x.answerLickTime - x.poleUpOnsetTime, u.trials(trialInds)));
            lickRate(ei,ai,si) = nanmean(cellfun(@(x,y,z) (y-z)/(x(y) - x(z)), allLickTimes, answerLickInds, firstLickInds));
            wastedLicks(ei,ai,si) = mean(cell2mat(answerLickInds) - cell2mat(firstLickInds));    
        end
    end
end
%
figure, hold on
shadedErrorBar(angles, mean(lickRate(:,:,1)), std(lickRate(:,:,1))/sqrt(length(expertInds)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(lickRate(:,:,2)), std(lickRate(:,:,2))/sqrt(length(expertInds)), 'lineprop', 'r-')

xlabel('Angle (\circ)'), ylabel('(Hz)')
xticks(angles)
title('Lick rate')
%%
figure, hold on
shadedErrorBar(angles, mean(correctAnswerLickTime(:,:,1)), std(correctAnswerLickTime(:,:,1))/sqrt(length(expertInds)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(correctAnswerLickTime(:,:,2)), std(correctAnswerLickTime(:,:,2))/sqrt(length(expertInds)), 'lineprop', 'r-')

xlabel('Angle (\circ)'), ylabel('Time after pole onset (s)')
xticks(angles)
title('Correct answer lick time')

%% (3) Response
% How many licks are wasted?

figure, hold on
shadedErrorBar(angles, mean(wastedLicks(:,:,1)), std(wastedLicks(:,:,1))/sqrt(length(expertInds)), 'lineprop', 'b-')
shadedErrorBar(angles, mean(wastedLicks(:,:,2)), std(wastedLicks(:,:,2))/sqrt(length(expertInds)), 'lineprop', 'r-')

xlabel('Angle (\circ)'), ylabel('# of licks')
xticks(angles)
title('# of licks between first and answer lick')