% baseDir = 'C:\Data\suite2p\';
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

expertInds = find(cellfun(@(x) length(x) == 2, sessions));
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
            numTouch(ei,ai,si) = mean(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunks)), u.trials(angleTouchTrialInds)));
            angleAnswerTrialInds = intersect(angleInds, answerTrialInds);
            meanTouchNum(ei,ai,si) = mean(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunks)), u.trials(angleAnswerTrialInds)));
            propTouchTrial(ei,ai,si) = length(find(cellfun(@(x) sum(cellfun(@(y) length(find(x.whiskerTime(y(1)) <= x.answerLickTime)), x.protractionTouchChunks)), u.trials(angleAnswerTrialInds)))) / length(angleAnswerTrialInds);
        end        
    end
end

%%

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