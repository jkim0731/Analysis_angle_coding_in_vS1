% Looking for decision neurons 
% inspired by Buetfering et al (Hausser lab)
% First, find out what is the distribution of time difference between first touch and first lick (before answer lick)

%% finding out the distribution of time difference between first touch and first lick
%% single sessions test

baseDir = 'Y:\Whiskernas\JK\suite2p\';
mouse = 52;
session = 25;

cd(sprintf('%s%03d',baseDir,mouse))
load(sprintf('UberJK%03dS%02d',mouse, session))

%% look at time difference between first touch and first lick

% about how much of the trials in a session have touch?
indTouchTrials = find(cellfun(@(x) length(x.protractionTouchChunks) + length(x.retractionTouchChunks), u.trials));
length(indTouchTrials)/length(u.trials)

%%
% about how much of the trials with lick in a session have touch?
indLick = intersect(find(cellfun(@(x) length(x.leftLickTime) + length(x.rightLickTime), u.trials)), find(cellfun(@(x) 1-strcmp(x.trialType, 'oo'), u.trials)));
indTouchTrialsLick = find(cellfun(@(x) length(x.protractionTouchChunks) + length(x.retractionTouchChunks), u.trials(indLick)));
length(indTouchTrialsLick) / length(indLick)


%% What is the correct rate of lick trials without touch?
% kind of confirmation about the touch detection method
indNTTLick = setdiff(indLick,indTouchTrialsLick);
responsesNTT = cellfun(@(x) x.response, u.trials(indNTTLick));
length(find(responsesNTT == 1)) / length(responsesNTT>=0)
% variable results.... high correct rate can come from the fact that
% absence of touch with a specific whisking pattern can discriminate the
% correct lick port quite well.

%% what is the distribution of time difference between first touch and first lick?
% only from lick trials
% the time difference is between pole in position to the first lick
timeFirstTouch = zeros(length(indLick),1);
for i = 1 : length(indLick)
    temptrial = u.trials{indLick(i)};
    tempTouchMat = [cell2mat(temptrial.protractionTouchChunks'); cell2mat(temptrial.retractionTouchChunks')];
    if isempty(tempTouchMat)
        timeFirstTouch(i) = temptrial.poleUpTime(1);
    else
        timeFirstTouch(i) = temptrial.whiskerTime(tempTouchMat(1));
    end
end
timeDistr = cellfun(@(x) min([x.leftLickTime;x.rightLickTime]), u.trials(indLick)) - timeFirstTouch;

%%
figure, histogram(timeDistr)

%% Look at trials of negative time differences
indExample = indLick(find(timeDistr <= -1));
tnExample = cellfun(@(x) x.trialNum, u.trials(indExample));
firstTouchFrames = cellfun(@(x) find(x.whiskerTime == x.poleUpTime(1)), u.trials(indExample));