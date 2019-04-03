% # of trials, # of touches in each angle in 2 different volumes
% Considering just unimodal (or monotonic) tuning for now 2018/06/27 JK

% load u and, if it exists, ANOVA tuning file
clear
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};  
% mice = [39,52,53,54,56];
% sessions = {[22],[3,21],[3],[3],[3]};  

        % settings
        angles = 45:15:135;
%         baseFrameNum = 3;
        baseFrameDuration = 0.5; % in sec
        afterFrameDuration = 0.5; % in sec
        
        thresholdAnovaP = 0.05; 
        thresholdTtestSharpness = 0.05; 
        thresholdTtestResponse = 0.05;
        thresholdCategory = 0.05;
        thresholdResponseTotal = 0.05;
%         thresholdResponsePositive = 1 * afterFrameDuration/frameRate; % spk - base spk
%         thresholdResponseNegative = -1 * afterFrameDuration/frameRate; % spk - base spk
        
        anovactype = 'hsd';
        excludeDrinkingTime = 0;
        onlyBeforeDecision = 0;
%         if onlyBeforeDecision 
            onlyAfterDecision = 0;
%         else
%             onlyAfterDecision = 1;
%         end
        
        onlyFirstTouch = 1;
        allowOverlap = 0;
        
for mi = 1 : length(mice)
% for mi = 3
    for si = 1 : length(sessions{mi})
%     for si = 1

        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));
        
        % still some settings
        frameRate = u.frameRate;
        
        
        thresholdResponsePositive = 1; % spk - base spk
        thresholdResponseNegative = -1; % spk - base spk

        
        baseFrameNum = round(frameRate*baseFrameDuration);
        afterFrameNum = round(frameRate*afterFrameDuration);
        savefn = [u.mouseName,u.sessionName,'singleCell_anova_spk.mat']; %
%         if exist(savefn, 'file')
%             fprintf('JK%03dS%02d file already exists. Abort.\n', mice(mi), sessions{mi}(si))
%             continue
%         end
        % initialization
        cellsTuned = [];
        tuneAngle = []; % Max response angle, in case of multimodal response
        tuneDirection = []; % 1 up, 2 down, 3 bipolar
        tuneAmplitude = []; % response to the tuned angle.
        tuneModulationMaxmin = []; % -2~2. abs > 2 means bipolar. sign shows which direction is stronger, and value shows how strong it's bipolarity (both large increase and large reduction).
        tuneSharpness = []; % based on posthoc analysis. How many neighboring bins are significantly indifferent from the max tuned bin.
        tuneSharpnessFWHM = []; % # of neighbors with response larger than half of the modulation. For now it will be useful for only within-angle comparison
        tuneSharpnessTtest = []; % # of neighbors with same tune direction. For now it will be useful for only within-angle comparison
        tuneResponseProb = []; % proportion of response to the tuned angle
        tuneReliability = []; % To the max response angle. Average of Pearson's correlation to the average time series
        
        tuneSingle = [];
        tuneBroad = []; 
        tuneLOO = []; % leave-one-out
        tuneMM = []; % multimodal. Including bipolar.
        tuneCateg = [] ; % categorical (>= 90 or <= 90)
        tuneRamp = []; % ramping up or down. One may be off.
        tuneRampStrict = []; % from all means

        cellsNTResponse = [];
        cellsNTNR = [];
        NTRAmplitude = [];
        NTRdirection = []; % 1 up, 2 down
        NTRresponseProb = []; % proportion of response to the touch
        NTRreliability = []; % average of Pearson's correlation to the average time series (from all the angles)
        
        tuneMaxResponseTimepoint = [];
        NTRMaxResponseTimepoint = [];
        
        % for confirmation
        spkTotal = cell(length(u.cellNums),1); % for all cells. ordered by cell number        
        cellsTotal = u.cellNums;
        cellDepths = u.cellDepths;
        cellMaps = u.cellmap;
        cellisC2 = u.isC2;
        cellThresholds = zeros(length(u.cellNums),2); % 1st for positive, 2nd for negative
        oneSampleH = zeros(length(u.cellNums), length(angles));
        anovaP = zeros(length(u.cellNums),1);
        responseP = zeros(length(u.cellNums),1);
        meanTotal = zeros(length(u.cellNums),1);
        
                
        for cellid = 1:length(u.cellNums)
            disp(['Processing cell id ', num2str(cellid)])
            cellNum = u.cellNums(cellid);
            touchNumTrial = cell(length(angles),1);
            touchNumChunk = cell(length(angles),1);
            countedTouchs = cell(length(angles),1);
            spk = cell(length(angles),1);
            for ai = 1 : length(angles)
                angle = angles(ai);
                plane = floor(cellNum/1000);
                trialPlaneInd = find(cellfun(@(x) ~isempty(find(x.planes == plane, 1)), u.trials));
                trialAngleInd = find(cellfun(@(x) x.angle == angle, u.trials));
                trialInd = intersect(trialPlaneInd, trialAngleInd);

                touchNum = 0;
                touchNumTrial{ai} = [];
                touchNumChunk{ai} = [];
                for i = 1 : length(trialInd)
                    if ~isempty(u.trials{trialInd(i)}.protractionTouchChunks)
                        if onlyFirstTouch
                            if excludeDrinkingTime && ~isempty(u.trials{trialInds(i)}.drinkingTime)
                                if u.trials{trialInds(i)}.protractionTouchChunks{1}(1) + afterFrameNum / frameRate > u.trials{trialInds(i)}.drinkingTime
                                    continue
                                else
                                    touchChunks = {u.trials{trialInd(i)}.protractionTouchChunks{1}};
                                    chunkIndSave = 1;
                                end
                            else
                                touchChunks = {u.trials{trialInd(i)}.protractionTouchChunks{1}};
                                chunkIndSave = 1;
                            end
                        else                            
                            if excludeDrinkingTime && ~isempty(u.trials{trialInd(i)}.drinkingTime)
                                tempChunksDrinking = u.trials{trialInd(i)}.protractionTouchChunks;
                                chunkIndDrinking = 1:length(tempChunksDrinking);
                                drinkInd = [];
                                for j = 1 : length(tempChunksDrinking)
                                    if tempChunksDrinking{j}(1) > u.trials{trialInd(i)}.drinkingTime(1) - afterFrameDuration || tempChunksDrinking{j}(1) < u.trials{trialInd(i)}.drinkingTime(2) + afterFrameDuration
                                        drinkInd = [drinkInd, j];
                                    end
                                end
                                chunkIndDrinking = setdiff(chunkIndDrinking, drinkInd);
                                touchChunks = cell(length(chunkIndDrinking),1);
                                for j = 1 : length(touchChunks)
                                    touchChunks{j} = tempChunksDrinking{chunkIndDrinking(j)};                        
                                end
                                chunkIndSave = chunkIndDrinking;
                            else
                                touchChunks = u.trials{trialInd(i)}.protractionTouchChunks;
                                chunkIndSave = 1 : length(u.trials{trialInd(i)}.protractionTouchChunks);
                            end

                            if onlyBeforeDecision && ~isempty(u.trials{trialInd(i)}.answerLickTime)
                                decisionChunkInd = length(touchChunks)+1;
                                for j = 1 : length(touchChunks)
                                    if touchChunks{j}(1) > u.trials{trialInd(i)}.answerLickTime
                                        decisionChunkInd = j;
                                        break
                                    end
                                end
                                if decisionChunkInd == 1
                                    touchChunks = {};
                                    chunkIndSave = [];
                                else
                                    decisionInd = 1 : decisionChunkInd - 1;
                                    tempTouchChunks = touchChunks;
                                    touchChunks = cell(length(decisionInd),1);
                                    for j = 1 : length(decisionInd)
                                        touchChunks{j} = tempTouchChunks{decisionInd(j)};
                                    end
                                    chunkIndSave = chunkIndSave(decisionInd);
                                end
                            end

                            if onlyAfterDecision && ~isempty(u.trials{trialInd(i)}.answerLickTime)
                                decisionChunkInd = length(touchChunks)+1;
                                for j = 1 : length(touchChunks)
                                    if touchChunks{j}(1) > u.trials{trialInd(i)}.answerLickTime
                                        decisionChunkInd = j;
                                        break
                                    end
                                end
                                if decisionChunkInd == length(touchChunks)
                                    touchChunks = {};
                                    chunkIndSave = [];
                                else
                                    decisionInd = decisionChunkInd : length(touchChunks);
                                    tempTouchChunks = touchChunks;
                                    touchChunks = cell(length(decisionInd),1);
                                    for j = 1 : length(decisionInd)
                                        touchChunks{j} = tempTouchChunks{decisionInd(j)};
                                    end
                                    chunkIndSave = chunkIndSave(decisionInd);
                                end
                            end
                        end

                        if ~isempty(touchChunks)
                            touchNum = touchNum + length(touchChunks);
                            touchNumTrial{ai} = [touchNumTrial{ai}; ones(length(touchChunks),1) * trialInd(i)];
                            touchNumChunk{ai} = [touchNumChunk{ai}; chunkIndSave'];
                        end
                    end
                end
                %%
                countedTouchs{ai} = 0;
                spk{ai} = nan(touchNum, baseFrameNum + afterFrameNum);                
                for i = 1 : touchNum
                    trial = touchNumTrial{ai}(i);
                    chunk = touchNumChunk{ai}(i);
                    trialSpk = u.trials{trial}.spk;
                    time = u.trials{trial}.protractionTouchChunks{chunk}(1); % for now, just consider the first touch point
                    frameInd = find(u.trials{trial}.tpmTime{mod(plane-1,4)+1} >= time, 1, 'first') - 1; % -1 to adjust 0 s point 2018/12/11 JK
                    if i > 1 && ~allowOverlap
                        if trial == trialOld
                            if frameInd - frameOld < afterFrameNum
                                continue
                            end
                        end
                    end
                    countedTouchs{ai} = countedTouchs{ai} + 1;
                    spkInd = find(u.trials{trial}.neuindSession == cellNum);
                    if frameInd < baseFrameNum
                        tempBaseNum = frameInd;
                        spk{ai}(i,baseFrameNum-frameInd+1:baseFrameNum) = trialSpk(spkInd, 1:frameInd);
                    else
                        spk{ai}(i,1:baseFrameNum) = trialSpk(spkInd, frameInd-baseFrameNum + 1 : frameInd);
                    end

                    maxafternum = max(size(trialSpk,2)-frameInd, afterFrameNum);

                    if size(trialSpk,2)-frameInd < afterFrameNum
                        tempAfter = size(trialSpk,2)-frameInd;
                        spk{ai}( i, baseFrameNum + 1 : baseFrameNum + tempAfter ) = trialSpk(spkInd, frameInd + 1 : frameInd + tempAfter);
                    else
                        spk{ai}( i, baseFrameNum + 1 : end) = trialSpk(spkInd, frameInd + 1 : frameInd + afterFrameNum);
                    end

                    tempBase = nanmean(spk{ai}(i, 1 : baseFrameNum));
                    spk{ai}(i,:) = (spk{ai}(i,:) - tempBase);
                    frameOld = frameInd;
                    trialOld = trial;
                end
                spk{ai}(isnan(nanmean(spk{ai},2)),:) = []; % removing all the NaN rows
            end
            spkTotal{cellid} = spk;
            
            % ANOVA
            touchNumGroups = cellfun(@(x) size(x,1), spk);
            timeAverageSpk = zeros(sum(touchNumGroups),1);
            anovaGroups = zeros(sum(cellfun(@(x) size(x,1), spk)),1);
            frames = 1 : afterFrameNum;
            groupAverageSpk = zeros(length(angles),length(frames));
            for ai = 1 : length(angles)
                timeAverageSpk( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = nanmean(spk{ai}(:,baseFrameNum+frames),2);
                anovaGroups( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = deal(ai);

                groupAverageSpk(ai,:) = nanmean(spk{ai}(:,baseFrameNum+frames),1);
            end
            
            [anovaP(cellid), ~, anovaStat] = anova1(timeAverageSpk, anovaGroups, 'off');
            pairComp = multcompare(anovaStat, 'Ctype', anovactype, 'Display', 'off');                        
            statMeans = cellfun(@(x) mean(nanmean(x(:,baseFrameNum+frames),2)), spk);
            
            oneSampleH(cellid,:) = cellfun(@(x) ttest(nanmean(x(:,baseFrameNum+frames),2)), spk);
            if ~isempty(find(isnan(oneSampleH(cellid,:))))
                oneSampleH(cellid,find(isnan(oneSampleH(cellid,:)))) = deal(0);
            end
            oneSampleInd = find(oneSampleH(cellid,:));
            
            [~, responseP(cellid)] = ttest(timeAverageSpk); % regardless of the angles.
            meanTotal(cellid) = mean(timeAverageSpk);
            tempThreshold = [thresholdResponsePositive, thresholdResponseNegative];
            cellThresholds(cellid,:) = tempThreshold;
            
            if anovaP(cellid) <= thresholdAnovaP && ~isempty(oneSampleInd) 
                if max(statMeans(oneSampleInd)) > tempThreshold(1) || min(statMeans(oneSampleInd)) < tempThreshold(2)                    
                    cellsTuned = [cellsTuned; cellNum];
                    tuneModulationMaxmin = [tuneModulationMaxmin; (max(statMeans) - min(statMeans))];
                    
                    if max(statMeans(oneSampleInd)) > tempThreshold(1) % have to breakdown because the absolute threshold is different 2018/10/19 JK
                        if min(statMeans(oneSampleInd)) < tempThreshold(2) % bipolar % Response to major tuned angle is positive
                            tuneDirection = [tuneDirection; 3]; % bipolar

                            threshinds = find(statMeans(oneSampleInd) > tempThreshold(1) | statMeans(oneSampleInd) < tempThreshold(2));                            
                            [maxAmp, maxInd] = max(abs(statMeans(oneSampleInd(threshinds))));
                            tuneAmplitude = [tuneAmplitude; maxAmp];
                            tunedAngleInd = oneSampleInd(threshinds(maxInd)); % for readability
                            tuneAngle = [tuneAngle; angles(tunedAngleInd)];                            
                            tempTimeseries = groupAverageSpk(tunedAngleInd,:);
                            if statMeans(tunedAngleInd) > tempThreshold(1) % major response positive
                                mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            elseif statMeans(tunedAngleInd) < tempThreshold(2) % majore response negative
                                mt = find(tempTimeseries < min(tempTimeseries)*0.9,1);
                            else % error
                                mt = -1;
                            end
                        else % unidirectional
                            tuneDirection = [tuneDirection; 1]; % positive
                            [maxAmp, maxInd] = max(statMeans(oneSampleInd));
                            tuneAmplitude = [tuneAmplitude; maxAmp];
                            tunedAngleInd = oneSampleInd(maxInd); % for readability
                            tuneAngle = [tuneAngle; angles(tunedAngleInd)];
                            tempTimeseries = groupAverageSpk(tunedAngleInd,:);
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                        end
                    elseif min(statMeans(oneSampleInd)) < tempThreshold(2) % No positive response present. Purely negative response                        
                        tuneDirection = [tuneDirection; 2]; % decrease
                        [maxAmp, maxInd] = min(statMeans(oneSampleInd));
                        tuneAmplitude = [tuneAmplitude; maxAmp];
                        tunedAngleInd = oneSampleInd(maxInd); % for readability
                        tuneAngle = [tuneAngle; angles(tunedAngleInd)];
                        tempTimeseries = groupAverageSpk(tunedAngleInd,:);
                        mt = find(tempTimeseries < min(tempTimeseries)*0.9,1);
                    else % error
                        tuneDirection = [tuneDirection; 0]; % error
                        tuneAmplitude = [tuneAmplitude; 0];
                        tuneAngle = [tuneAngle; 0];
                        mt = -1;
                    end
                    tuneMaxResponseTimepoint = [tuneMaxResponseTimepoint; mt / frameRate];

                    % Response probability
                    spktuned = spk{tunedAngleInd}(:,baseFrameNum+1 : end); % for readability
                    responseNum = 0;
                    for ri = 1:size(spktuned,1)
                        tempTS = spktuned(ri,:);
                        tempTS = tempTS(~isnan(tempTS));
                        [~, p] = ttest(tempTS');
                        if p < thresholdTtestResponse
                            if (statMeans(tunedAngleInd) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                    (statMeans(tunedAngleInd) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                responseNum = responseNum + 1;
                            end
                        end
                    end
                    tuneResponseProb = [tuneResponseProb; responseNum/ri*100];

                    % Response reliability
                    template = nanmean(spk{tunedAngleInd}(:,baseFrameNum:end)); % include the last baseline frames
                    rho = zeros(size(spktuned,1),1);
                    for ri = 1 : length(rho)
                        tempTS = spk{tunedAngleInd}(ri,baseFrameNum:end);
                        inds = find(~isnan(tempTS));
                        rho(ri) = corr(template(inds)',tempTS(inds)');
                    end
                    tuneReliability = [tuneReliability; mean(rho)];

                    % Categorization
                    ind__1 = find(pairComp(:,1) == tunedAngleInd);
                    ind__2 = find(pairComp(:,2) == tunedAngleInd);
                    testInd = union(ind__1, ind__2);
                    insigInd = find(pairComp(testInd,6) >= thresholdCategory);
                    sigInd = find(pairComp(testInd,6) < thresholdCategory);
                    temp = pairComp(testInd(insigInd),1:2);
                    insigIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd, except when there's nothing
                    
                    if isempty(insigIndGroup)
                        tuneSharpness = [tuneSharpness; 1];
                        tuneSingle = [tuneSingle; cellNum];
                    else
                        broadInd = intersect(oneSampleInd,insigIndGroup);
                        if length(broadInd) < 2
                            tuneSharpness = [tuneSharpness; 1];
                            tuneSingle = [tuneSingle; cellNum];
                        else
                            tuneSharpness = [tuneSharpness; length(broadInd)];
                            broadNum = 1;
                            for tunei = tunedAngleInd-1:-1:1
                                if ismember(tunei, broadInd)
                                    broadNum = broadNum + 1;
                                else 
                                    break
                                end
                            end
                            for tunei = tunedAngleInd+1:length(angles)
                                if ismember(tunei, broadInd)
                                    broadNum = broadNum + 1;
                                end
                            end
                            if broadNum > 2
                                tuneBroad = [tuneBroad; cellNum];
                            end
                        end
                    end

                    temp = pairComp(testInd(sigInd),1:2);
                    sigIndGroup = setdiff(temp(:), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.
                    if ~isempty(find(diff(insigIndGroup)>1,1))
                        if sum(oneSampleH(cellid,sigIndGroup))
                            tuneMM = [tuneMM; cellNum]; % multimodal. Including bipolar.
                        end
                        if length(sigIndGroup) == 1 && ... % only one bin is significantly different from the tuned bin. (can't be larger in response because of the way tuned bin is defined)
                                all(oneSampleH(cellid,insigIndGroup)) % and all insignicant indices are different from 0
                            tuneLOO = [tuneLOO; cellNum]; % leave-one-out. Part of multimodal in definition.
                        end                            
                    end
                    
                    center = (length(angles)+1) / 2;
                    compInd = union(find(pairComp(:,1) == tunedAngleInd), find(pairComp(:,2) == tunedAngleInd));
                    indMat = pairComp(compInd,1:2);
                    if tunedAngleInd < center
                        withinInd = unique(mod( setdiff( find(indMat < center), find(indMat == tunedAngleInd) ) , size(indMat,1)));
                        withinInd(withinInd==0) = size(indMat,1);
                        betweenInd = unique(mod( find(indMat > center) , size(indMat,1) ));
                        betweenInd(betweenInd==0) = size(indMat,1);
                    else
                        withinInd = unique(mod( setdiff( find(indMat > center), find(indMat == tunedAngleInd) ) , size(indMat,1)));
                        withinInd(withinInd==0) = size(indMat,1);
                        betweenInd = unique(mod( find(indMat < center) , size(indMat,1) ));
                        betweenInd(betweenInd==0) = size(indMat,1);
                    end
                    if isempty(find(pairComp(compInd(withinInd),6) < thresholdCategory, 1)) && ... % nothing within the same half is different from the max ind
                            isempty(find(pairComp(compInd(betweenInd),6) >= thresholdCategory, 1)) % nothing between different half is same with the max ind
                        tuneCateg = [tuneCateg; cellNum] ; % categorical (>= 90 or <= 90)
                    end

                    if isempty(find(diff(sign(diff(statMeans))),1)) % everything is going up or down 
                        tuneRampStrict = [tuneRampStrict; cellNum]; % ramping up or down
                    end
                    
                elseif responseP(cellid) < thresholdResponseTotal
                    if meanTotal(cellid) > tempThreshold(1) || meanTotal(cellid) < tempThreshold(2)
                        spkAll = cell2mat(cellfun(@(x) x, spk, 'uniformoutput', false));
                        cellsNTResponse = [cellsNTResponse; cellNum];
                        tempTimeseries = nanmean(spkAll(:,baseFrameNum+frames));
                        NTRAmplitude = [NTRAmplitude; abs(mean(timeAverageSpk))];
                        if meanTotal(cellid) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif meanTotal(cellid) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries < min(tempTimeseries)*0.9,1);                            
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; 0]; % error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(spkAll,1)
                            tempTS = spkAll(ri,baseFrameNum + frames);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(timeAverageSpk) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(timeAverageSpk) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(spkAll(:,baseFrameNum:end));
                        rho = zeros(size(spkAll,1),1);
                        for ri = 1 : size(spkAll,1)
                            tempTS = spkAll(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cellNum];
                    end
                else
                    cellsNTNR = [cellsNTNR; cellNum];
                end
            else
                if responseP(cellid) < thresholdResponseTotal
                    if meanTotal(cellid) > tempThreshold(1) || meanTotal(cellid) < tempThreshold(2)
                        spkAll = cell2mat(cellfun(@(x) x, spk, 'uniformoutput', false));
                        cellsNTResponse = [cellsNTResponse; cellNum];
                        tempTimeseries = nanmean(spkAll(:,baseFrameNum+frames));
                        NTRAmplitude = [NTRAmplitude; abs(mean(timeAverageSpk))];
                        if meanTotal(cellid) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif meanTotal(cellid) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries < min(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; -1]; % another kind of error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(spkAll,1)
                            tempTS = spkAll(ri,baseFrameNum + frames);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(timeAverageSpk) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(timeAverageSpk) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(spkAll(:,baseFrameNum:end));
                        rho = zeros(size(spkAll,1),1);
                        for ri = 1 : size(spkAll,1)
                            tempTS = spkAll(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cellNum];
                    end
                else
                    cellsNTNR = [cellsNTNR; cellNum];
                end
            end
        end
        
        noise = u.noise;
        celly = u.celly;
        cellx = u.cellx;
        c2ypoints = u.c2ypoints;
        c2xpoints = u.c2xpoints;
        fovsize = u.fovsize;
        fovxrange = u.fovxrange;
        fovyrange = u.fovyrange;
        fovdepth = u.fovdepth;

        save(savefn, 'cell*','tune*','NTR*', 'spkTotal', '*ctype', 'angles', 'baseFrameNum', 'afterFrameNum', 'threshold*', ...
            'excludeDrinkingTime', 'onlyBeforeDecision', 'onlyAfterDecision', 'allowOverlap', 'onlyFirstTouch','frameRate', 'oneSampleH', '*P', 'noise', 'c2*points', 'fov*')
%         end
    end
end