baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
% sessions = {[19],[3,10],[3,21],[1,17],[7],[2],[22],[3],[3,21],[3],[3],[3]};
for mi = 1
    for si = 1
        mouse = mice(mi);
        session = sessions{mi}(si);
        
        posShiftTouch = 2;
        posShiftSound = 3;
        posShiftReward = 3;
        posShiftWhisking = 4;
        posShiftLicking = 1;
        posShift = 4; % maximum posShift
        negShift = 2;
        testPortion = 0.3; % 30 % test set

        dn = sprintf('%s%03d',baseDir,mouse);
        ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
        cd(dn)
%         if exist(ufn,'file')
            load(ufn)
%         else
%             u = Uber.buildUberArray(mouse, session);
%         end
        frameRate = u.frameRate;
        
        angles = unique(cellfun(@(x) x.angle, u.trials));
        distances = unique(cellfun(@(x) x.distance, u.trials));

        touchGroup = cell(2,1);
        choiceGroup = cell(3,1);
        angleGroup = cell(length(angles),1);
        distanceGroup = cell(length(distances),1);
        timeGroup = cell(3,1); % dividing whole session into 5 different time points

        ptouchGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) length(x.protractionTouchChunks), u.trials))));
        ptouchGroup{2} = setdiff(u.trialNums, ptouchGroup{1});

        rtouchGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) length(x.retractionTouchChunks), u.trials))));
        rtouchGroup{2} = setdiff(u.trialNums, rtouchGroup{1});

        choiceGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 1, u.trials))));
        choiceGroup{2} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 0, u.trials))));
        choiceGroup{3} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == -1, u.trials))));

        for i = 1 : length(angles)
            angleGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.angle == angles(i), u.trials))));
        end

        for i = 1 : length(distances)
            distanceGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.distance == distances(i), u.trials))));
        end
        %% Design matrices
        % standardized using all the trials
        allPredictors = cell(8,1);
        allPredictorsMean = cell(8,1);
        allPredictorsStd = cell(8,1);
        nani = cell(8,1);
%                 trainingPredictorInd = cell(8,1);
%                 testPredictorInd = cell(8,1);
%                 trainingInputMat = cell(8,1);
%                 testInputMat = cell(8,1);

        for cgi = 1:2 % cell group index
            tindcell = find(cellfun(@(x) ismember(1001+(cgi-1)*4000, x.neuindSession), u.trials));

            tind = tindcell;
            for plane = 1 : 4    

%                         trainingPredictorInd{(cgi-1)*4 + plane} = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, trainingTn), u.trials(tind)','uniformoutput',false));
%                         testPredictorInd{(cgi-1)*4 + plane} = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, testTn), u.trials(tind)','uniformoutput',false));

                pTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) x.whiskerTime(y(1)), x.protractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                pTouchFrame = pTouchCount;
                pTouchFrame(pTouchFrame > 0) = 1;
                pTouchFrameAngles = cell(length(angles)+1,1);
                for ai = 1 : length(angles)
                    tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                    pTouchFrameAngles{ai} = pTouchFrame .* tempAngleBinary';
                end
                pTouchFrameAngles{end} = pTouchFrame;
                scPoleup = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                drink = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime:0.03:x.drinkingOnsetTime+1, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                drink(drink>0) = 1;
                drinkAngles = cell(length(angles)+1,1);
                for ai = 1 : length(angles)
                    tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                    drinkAngles{ai} = drink .* tempAngleBinary';
                end
                drinkAngles{end} = drink;
                lLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                rLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                %%
                whiskingOnsetCell = cell(1,length(tind));
                whiskingAmplitudeCell = cell(1,length(tind));
                whiskingMidpointCell = cell(1,length(tind));

                for ti = 1 : length(tind)
                    currTrial = u.trials{tind(ti)};
                    time = [0, currTrial.tpmTime{plane}];
                    wtimes = currTrial.whiskerTime;
                    if iscell(wtimes)
                        wtimes = wtimes{1};
                    end
                    [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
                    whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration; % in s
                    onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
                    tempOnset = histcounts(onsetTimes, time);
                    whiskingOnsetCell{ti} = [nan(1,posShift), tempOnset, nan(1,posShift)];

                    tempMid = zeros(1,length(time)-1);
                    tempAmp = zeros(1,length(time)-1);
                    for i = 1 : length(tempMid)
                        startInd = find(wtimes >= time(i), 1, 'first');
                        endInd = find(wtimes < time(i+1), 1, 'last');
                        tempMid(i) = mean(midpoint(startInd:endInd));
                        tempAmp(i) = max(amplitude(startInd:endInd));
                    end
                    tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
                    tempAmp(isnan(tempAmp)) = deal(mode(tempAmp(isfinite(tempAmp))));
                    whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];
                    whiskingAmplitudeCell{ti} = [nan(1,posShift), tempAmp, nan(1,posShift)];
                end
                whiskingOnset = cell2mat(whiskingOnsetCell);
                whiskingMidpoint = cell2mat(whiskingMidpointCell);
                whiskingAmplitude = cell2mat(whiskingAmplitudeCell);

                %%
                pTouchFrameMat = zeros(length(pTouchFrame), (posShiftTouch + 1) * (length(angles)+1));

                scPoleUpMat = zeros(length(scPoleup), posShiftSound + 1);
                drinkAnglesMat = zeros(length(drink), (posShiftReward + 1) * (length(angles)+1));
                for i = 1 : posShiftTouch + 1
                    for ai = 1 : length(angles) + 1
                        pTouchFrameMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchFrameAngles{ai}, [0 i-1])';
                    end
                end
                for i = 1 : posShiftSound + 1
                    scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                end
                for i = 1 : posShiftReward + 1
                    for ai = 1 : length(angles) + 1
                        drinkAnglesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(drinkAngles{ai}, [0 i-1])';
                    end                                
                end

                whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShiftWhisking + 1);
                whiskingAmplitudeMat = zeros(length(whiskingAmplitude), negShift + posShiftWhisking + 1);
                whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShiftWhisking + 1);

                lLickMat = zeros(length(lLick), negShift + posShiftLicking + 1);
                rLickMat = zeros(length(rLick), negShift + posShiftLicking + 1);

                for i = 1 : negShift + posShiftWhisking + 1
                    whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                    whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
                    whiskingAmplitudeMat(:,i) = circshift(whiskingAmplitude, [0 -negShift + i - 1])';                            
                end

                for i = 1 : negShift + posShiftLicking + 1
                    lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                    rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                end

                touchMat = [pTouchFrameMat];
                soundMat = [scPoleUpMat];
                drinkMat = drinkAnglesMat;
                whiskingMat = [whiskingOnsetMat, whiskingAmplitudeMat, whiskingMidpointMat];
                lickingMat = [lLickMat, rLickMat];
                allPredictors{(cgi-1)*4 + plane} = [touchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                nani{(cgi-1)*4 + plane} = find(nanstd(allPredictors{(cgi-1)*4 + plane})==0);
                allPredictorsMean{(cgi-1)*4 + plane} = nanmean(allPredictors{(cgi-1)*4 + plane});
                allPredictorsStd{(cgi-1)*4 + plane} = nanstd(allPredictors{(cgi-1)*4 + plane});
                % normalization of all predictors
                allPredictors{(cgi-1)*4 + plane} = (allPredictors{(cgi-1)*4 + plane} - nanmean(allPredictors{(cgi-1)*4 + plane})) ./ nanstd(allPredictors{(cgi-1)*4 + plane});
                allPredictors{(cgi-1)*4 + plane}(:,nani{(cgi-1)*4 + plane}) = deal(0);
%                         trainingInputMat{(cgi-1)*4 + plane} = allPredictors{(cgi-1)*4 + plane}(find(trainingPredictorInd{(cgi-1)*4 + plane}),:);
%                         testInputMat{(cgi-1)*4 + plane} = allPredictors{(cgi-1)*4 + plane}(find(testPredictorInd{(cgi-1)*4 + plane}),:);
            end
        end
    end
end

%%
apTest = allPredictors;
load('glmResponseType_JK025S04_m44_R10','allPredictors')

%%
compare(apTest{1}(find(~isnan(apTest{1}))), allPredictors{1}(find(~isnan(allPredictors{1}))))
%%
col = 82;
compare(apTest{1}(:,col), allPredictors{1}(:,col))
%%
load('glmResponseType_JK025S04_m44_R04','allPredictors')
ap1 = allPredictors;
load('glmResponseType_JK025S04_m44_R03','allPredictors')
ap2 = allPredictors;
compare(ap1{1}, ap2{1})

%% The difference between allPredictors came from whisking parameters.
% BUT, if the spikes are the same, then the result should be the same.
    

%% Test GLM results with new uber arrays.
baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
mi = 2;
si = 1;
mouse = mice(mi);
session = sessions{mi}(si);
ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
load(ufn)
% cIDAll = u.cellNums;
%%
numCell = length(cIDAll);
            
tindcellAll = cell(numCell,1);
cindAll = zeros(numCell,1);
planeIndAll = zeros(numCell,1);
for i = 1 : numCell
    tindcellAll{i} = find(cellfun(@(x) ismember(cIDAll(i), x.neuindSession), u.trials));
    cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == cIDAll(i));
    planeIndAll(i) = floor(cIDAll(i)/1000);
end
spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);
%%
cellnum = 1;
cind = cindAll(cellnum);
tindCell = tindcellAll{cellnum};
testInd = find(ismember(u.trialNums,testTn{cellnum}));
iTest = intersect(tindCell, testInd);
planeInd = planeIndAll(cellnum);
plane = mod(planeInd,4);
if plane==0
    plane = 4;
end
                    
testPredictorInd = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, testTn{cellnum}), u.trials(tindCell)','uniformoutput',false));
testInput = allPredictors{planeInd}(find(testPredictorInd),:);

spkTemp = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
spkTemp = spkTemp';
finiteIndTest = intersect(find(isfinite(spkTemp)), find(isfinite(sum(testInput,2))));
spkTest = spkTemp(finiteIndTest)';

% spkTemp2 = 


model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*fitCoeffs{cellnum});
mu = mean(spkTest); % null poisson parameter
nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);

% fitDeviance(cellnum) - devianceFullNull
fitDevExplained(cellnum) - 1 + (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood)
