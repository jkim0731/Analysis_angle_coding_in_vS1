% Angle-tuning from touch model. 
% Remove all the other behavior categories.
% It will be touch-only reconstructed model. 
% It will be used to calculate impact on angle-tuning by other behaviors,
% and also for angle-selectivity across training.
% 2020/06/21 JK

clear

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

angles = 45:15:135;
thresholdAnovaP = 0.05; 
thresholdPermAnovaP = 0.05;
anovactype = 'hsd';
numResampling = 10000;

% for mi = 1 : length(mice)
for mi = 1
    for si = 1
%     for si = 1 : length(sessions{mi})
        mouse = mice(mi);
        session = sessions{mi}(si);
        
        fprintf('Processing JK%03d S%02d.\n', mouse, session)
        
        savefn = sprintf('touchModel_angle_tuning_model_v9_JK%03dS%02d', mouse, session);
        
        tuning = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC_PTC', baseDir, mouse, mouse, session), 'spk');
        
        touchID = tuning.spk.touchID;
        tuned = tuning.spk.tuned;
        
        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session), 'u')
        load(sprintf('%s%03d\\glmResponseType_JK%03dS%02d_lasso_NC_R01',baseDir, mouse, mouse, session), 'cIDAll', 'fitCoeffs', 'allPredictors', 'indPartial');
        
        if ~isempty(find(cIDAll - u.cellNums))
            error('u.cellNums and cIDAll mismatch.')
        end
        tempInd = ismember(cIDAll, touchID);
        cIDAll = cIDAll(tempInd);
        if ~isempty(find(cIDAll - touchID'))
            error('touch ID and cell ID All mismatch.')
        end
        depth = u.cellDepths(tempInd);
        fitCoeffs = fitCoeffs(tempInd);
        
        coeffs = zeros(size(fitCoeffs,1),length(fitCoeffs{1}),10);
        for ri = 1 : 10
            glmfn = sprintf('%s%03d\\glmResponseType_JK%03dS%02d_lasso_NC_R%02d', baseDir, mouse, mouse, session, ri); % these are only from touch response cells
            load(glmfn, 'fitCoeffs')
            fitCoeffs = fitCoeffs(tempInd);
            emptyCoeffsInd = find(cellfun(@isempty, fitCoeffs));
            firstNonemptyInd = find(1-cellfun(@isempty, fitCoeffs),1,'first');
            if ~isempty(emptyCoeffsInd)
                for emptyi = 1 : length(emptyCoeffsInd)
                    fitCoeffs{emptyCoeffsInd(emptyi)} = nan(size(fitCoeffs{firstNonemptyInd}));
                end
            end
            coeffs(:,:,ri) = cell2mat(fitCoeffs')';
        end
        meanCoeffsTouch = nanmean(coeffs,3);

        % convert touch glm predictors back to each trial
        % there are 8 NaNs between trials. first trial starts with 4 NaNs, and the last ends with 4 NaNs.
        allPredictorsTouch = cell(8,1);
        for i = 1 : 8
            if ~isempty(allPredictors{i}) % it can be empty in case of JK027 S09 and S10
                naninds = find(isnan(allPredictors{i}(:,1)));
                nanindChanges = find(diff(naninds)>1);
                numTrials = length(nanindChanges);
                allPredictorsTouch{i} = cell(numTrials,1);
                for j = 1 : numTrials
                    allPredictorsTouch{i}{j} = allPredictors{i}(naninds(nanindChanges(j))+1:naninds(nanindChanges(j)+1)-1,:);
                end
            end
        end
        
        upperPlaneInd = find(cellfun(@(x) length(find(x.planes==1)), u.trials));
        lowerPlaneInd = find(cellfun(@(x) length(find(x.planes==5)), u.trials));
        upperPlaneSpikes = cell2mat(cellfun(@(x) x.spk, u.trials(upperPlaneInd)', 'un', 0));
        lowerPlaneSpikes = cell2mat(cellfun(@(x) x.spk, u.trials(lowerPlaneInd)', 'un', 0));
        upperTouchID = cIDAll(cIDAll < 5000);
        lowerTouchID = cIDAll(cIDAll > 5000);
        upperID = u.cellNums(u.cellNums < 5000);
        lowerID = u.cellNums(u.cellNums > 5000);
        upperTouchInd = find(ismember(upperID, upperTouchID)); % sorted
        lowerTouchInd = find(ismember(lowerID, lowerTouchID)); % sorted
        upperTouchSpikes = upperPlaneSpikes(upperTouchInd,:);
        lowerTouchSpikes = lowerPlaneSpikes(lowerTouchInd,:);
        touchSpikes = [mat2cell(upperTouchSpikes, ones(1,size(upperTouchSpikes,1)), size(upperTouchSpikes,2)); ...
            mat2cell(lowerTouchSpikes, ones(1,size(lowerTouchSpikes,1)), size(lowerTouchSpikes,2))];

        % then, all models (from touch neurons only)
        touchFullModels = cell(length(cIDAll),1);
        touchOnlyModels = cell(length(cIDAll),1);
        deFull = zeros(length(cIDAll),1);
        deTouchOnly = zeros(length(cIDAll),1);
        parfor ci = 1 : length(cIDAll)
            pi = floor(cIDAll(ci)/1000);
            tempPredictors = cell2mat(allPredictorsTouch{pi});
            tempInput = [ones(size(tempPredictors,1),1), tempPredictors];
            touchFullModels{ci} = exp(tempInput*meanCoeffsTouch(ci,:)');

            tempInputWhiskerOnly = tempInput;
            tempInputWhiskerOnly(:,indPartial{1}(end)+2:end) = 0;
            touchOnlyModels{ci} = exp(tempInputWhiskerOnly*meanCoeffsTouch(ci,:)');

            spkTest = touchSpikes{ci};
            mu = nanmean(spkTest);
            nullLogLikelihood = nansum(log(poisspdf(spkTest,mu)));
            saturatedLogLikelihood = nansum(log(poisspdf(spkTest,spkTest)));
            fullLogLikelihood = nansum(log(poisspdf(spkTest',touchFullModels{ci})));
            deFull(ci) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            touchOnlyLogLikelihood = nansum(log(poisspdf(spkTest',touchOnlyModels{ci})));
            deTouchOnly(ci) = (touchOnlyLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
        end
       
        % set up touch frames, baseline frames, indices, etc.
        answerTime = cell(length(u.trials),1);
        for di = 1 : length(answerTime)
            if isempty(u.trials{di}.answerLickTime)
                answerTime{di} = u.trials{di}.poleDownOnsetTime;
            else
                answerTime{di} = u.trials{di}.answerLickTime;
            end
        end
        tempTouchTrialInd = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        paTouchInd = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y, u.trials(tempTouchTrialInd), answerTime(tempTouchTrialInd)));
        touchTrialInd = tempTouchTrialInd(paTouchInd); % index of u
        numPlane = length(u.mimg);
        planeTouchTrialsInd = cell(numPlane,1); % pre-answer touch trials index (of u) in each plane
        planeTrialsInd = cell(numPlane,1); % all trials index (of u) in each plane, to match with predictors later
        poleUpFrames = cell(numPlane,1);
        beforePoleUpFrames = cell(numPlane,1);
        touchFrames = cell(numPlane,1);
        numTouchPreAnswer = cell(numPlane,1);
        angleTouchTrialInds = cell(numPlane,1);
        for pi = 1 : numPlane
            planeTrialsInd{pi} = find(cellfun(@(x) ismember(pi, x.planes), u.trials)); % having exactly same # as allPredictorsWKV
            planeTouchTrialsInd{pi} = intersect(planeTrialsInd{pi}, touchTrialInd); % index of u
            tempInd = find(u.trials{planeTouchTrialsInd{pi}(1)}.planes == pi);
            poleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} >= x.poleUpTime(1) & x.tpmTime{tempInd} <= x.poleUpTime(end)), u.trials(planeTouchTrialsInd{pi}), 'uniformoutput', false);
            beforePoleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} < x.poleUpOnsetTime), u.trials(planeTouchTrialsInd{pi}), 'uniformoutput', false);
            touchFrames{pi} = cell(length(planeTouchTrialsInd{pi}),1);
            numTouchPreAnswer{pi} = zeros(length(planeTouchTrialsInd{pi}),1);
            for ti = 1 : length(planeTouchTrialsInd{pi})
                tempTrial = u.trials{planeTouchTrialsInd{pi}(ti)};
                if isempty(tempTrial.answerLickTime)
                    tempAnswerTime = tempTrial.poleDownOnsetTime;
                else
                    tempAnswerTime = tempTrial.answerLickTime;
                end
                preAnswerInd = find(cellfun(@(x) tempTrial.whiskerTime(x(1)) < tempAnswerTime, tempTrial.protractionTouchChunksByWhisking));
                tempFrames = cell(1, length(preAnswerInd));
                for ptci = 1 : length(tempFrames)
                    tempFrames{ptci} = [0:1] + find(tempTrial.tpmTime{tempInd} >= tempTrial.whiskerTime(tempTrial.protractionTouchChunksByWhisking{ptci}(1)), 1, 'first');
                end
                touchFrames{pi}{ti} = unique(cell2mat(tempFrames));
                numTouchPreAnswer{pi}(ti) = length(preAnswerInd);
            end
            angleTouchTrialInds{pi} = cell(length(angles),1); % index of planeTouchTrialsInd{pi}
            for ai = 1 : length(angles)
                angleTouchTrialInds{pi}{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(planeTouchTrialsInd{pi})));
            end
        end
        
        spikeAngleAll = cell(length(cIDAll), 1);
        fullModelAngleAll = cell(length(cIDAll), 1);
        touchOnlyAngleAll = cell(length(cIDAll), 1);
        atSpikeAll = cell(length(cIDAll),1);
        atFullAll = cell(length(cIDAll),1);
        atTouchAll = cell(length(cIDAll),1);
        
        spkAnovaPAll = zeros(length(cIDAll),1);
        fullAnovaPAll = zeros(length(cIDAll),1);
        touchAnovaPAll = zeros(length(cIDAll),1);
        
        tunedSpk = zeros(length(cIDAll),1);
        tunedFull = zeros(length(cIDAll),1);
        tunedTouch = zeros(length(cIDAll),1);
        parfor ci = 1 : length(cIDAll)
%         for ci = 1 : length(cIDAll)
%             fprintf('cell # %d\n', ci)
            cID = cIDAll(ci);
            pi = floor(cID/1000);
            tempInd = find(cellfun(@(x) ismember(pi, x.planes), u.trials), 1, 'first');
            tempcIDlist = u.trials{tempInd}.neuindSession;
            cind = find(tempcIDlist == cID);
            spikeAngle = cell(1, length(angles));
            fullModelAngle = cell(1, length(angles));
            touchOnlyAngle = cell(1, length(angles));

            touchFramesAngle = cell(1, length(angles));
            baselineFramesAngle = cell(1, length(angles));
            numTouchAngle = cell(1, length(angles));
            for ai = 1 : length(angles)
                tempTouchTrialInd = angleTouchTrialInds{pi}{ai}; % index of planeTouchTrialsInd{pi}
                tempTrialInd = planeTouchTrialsInd{pi}(tempTouchTrialInd); % index of u
                spikeAngle{ai} = cell(length(tempTrialInd),1);
                % feature removal
                fullModelAngle{ai} = cell(length(tempTrialInd),1);
                touchOnlyAngle{ai} = cell(length(tempTrialInd),1);

                touchFramesAngle{ai} = cell(length(tempTrialInd),1);
                baselineFramesAngle{ai} = cell(length(tempTrialInd),1);
                numTouchAngle{ai} = cell(length(tempTrialInd),1);
                for ti = 1 : length(tempTrialInd)
                    spikeAngle{ai}{ti} = u.trials{tempTrialInd(ti)}.spk(cind,:);

                    tempPredictorInd = find(planeTrialsInd{pi} == tempTrialInd(ti)); % index matched with the order in allPredictorsWKV
                    tempPredictor = allPredictorsTouch{pi}{tempPredictorInd};
                    tempInput = [ones(size(tempPredictor,1),1), tempPredictor];

                    tempCoeff = meanCoeffsTouch(ci,:);
                    fullModelAngle{ai}{ti} = exp(tempInput * tempCoeff');

                    tempCoeff(1,indPartial{1}(end)+2:end) = 0;
                    touchOnlyAngle{ai}{ti} = exp(tempInput * tempCoeff');


                    touchFramesAngle{ai}{ti} = touchFrames{pi}{tempTouchTrialInd(ti)};
                    baselineFramesAngle{ai}{ti} = beforePoleUpFrames{pi}{tempTouchTrialInd(ti)};
                    numTouchAngle{ai}{ti} = numTouchPreAnswer{pi}(tempTouchTrialInd(ti));
                end
            end
            
            % trial-by-trial model reconstruction
            touchResponseSpike = cell(1, length(angles));
            touchResponseFull = cell(1, length(angles));
            touchResponseTouch = cell(1, length(angles));

            for ai = 1 : length(angles)
                touchResponseSpike{ai} = zeros(length(spikeAngle{ai}),1);
                
                touchResponseFull{ai} = zeros(length(spikeAngle{ai}),1);
                touchResponseTouch{ai} = zeros(length(spikeAngle{ai}),1);
                
                for ti = 1 : length(spikeAngle{ai})
                    touchResponseSpike{ai}(ti) = sum( spikeAngle{ai}{ti}(touchFramesAngle{ai}{ti})  - mean( spikeAngle{ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti};
                    
                    touchResponseFull{ai}(ti) = sum( fullModelAngle{ai}{ti}(touchFramesAngle{ai}{ti})  - nanmean( fullModelAngle{ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti}; % nanmean because models start with NaN values
                    touchResponseTouch{ai}(ti) = sum( touchOnlyAngle{ai}{ti}(touchFramesAngle{ai}{ti})  - nanmean( touchOnlyAngle{ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti};

                end
            end

            % save some values
            atSpike = cellfun(@mean, touchResponseSpike);
            atFull = cellfun(@mean, touchResponseFull);
            atWhisker = cellfun(@mean, touchResponseTouch);
            
            atSpikeAll{ci} = atSpike;
            atFullAll{ci} = atFull;
            atTouchAll{ci} = atWhisker;
            
            spikeAngleAll{ci} = touchResponseSpike;
            fullModelAngleAll{ci} = touchResponseFull;
            touchOnlyAngleAll{ci} = touchResponseTouch;
            
            % run anova and shuffling to determine angle-tuning
            % spk first, for confirmation (there could be a coding error)
            
            groupAnova = zeros(sum(cellfun(@length, touchResponseSpike)),1);
            angleLengths = [0,cumsum(cellfun(@length, touchResponseSpike))];
            for ai = 1 : length(angles)
                groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
            end
            
            if ~isempty(find(cellfun(@length, touchResponseSpike) - cellfun(@length, touchResponseFull)))
                error('# of trials mismatch between spikes and full model')
            elseif ~isempty(find(cellfun(@length, touchResponseSpike) - cellfun(@length, touchResponseTouch)))
                error('# of trials mismatch between spikes and full model')
            end
            
            spkAnovaVal = zeros(sum(cellfun(@length, touchResponseSpike)),1);
            for ai = 1 : length(angles)
                spkAnovaVal(angleLengths(ai)+1:angleLengths(ai+1)) = touchResponseSpike{ai};
            end
            if ~isempty(find(isnan(spkAnovaVal)))
                error('nan values')
            end
            
            spkAnovaP = anova1(spkAnovaVal, groupAnova, 'off');
            spkAnovaPAll(ci) = spkAnovaP;

            tempH = cellfun(@(x) ttest(x), touchResponseSpike);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            if spkAnovaP < thresholdAnovaP && ~isempty(sigInd)
                % permutation test
                permAnovaP = zeros(numResampling,1);
                for ri = 1 : numResampling
                    tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                    permAnovaP(ri) = anova1(spkAnovaVal, tempG, 'off');
                end
                
                if length(find(permAnovaP < spkAnovaP)) < 0.05 * numResampling % passing the test
                    tunedSpk(ci) = 1;
                end
            end
            
            
            
            % then full model 
            fullAnovaVal = zeros(sum(cellfun(@length, touchResponseFull)),1);
            for ai = 1 : length(angles)
                fullAnovaVal(angleLengths(ai)+1:angleLengths(ai+1)) = touchResponseFull{ai};
            end
            if ~isempty(find(isnan(fullAnovaVal)))
                error('nan values')
            end
            
            fullAnovaP = anova1(fullAnovaVal, groupAnova, 'off');
            fullAnovaPAll(ci) = fullAnovaP;

            tempH = cellfun(@(x) ttest(x), touchResponseFull);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            if fullAnovaP < thresholdAnovaP && ~isempty(sigInd)
                % permutation test
                permAnovaP = zeros(numResampling,1);
                for ri = 1 : numResampling
                    tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                    permAnovaP(ri) = anova1(fullAnovaVal, tempG, 'off');
                end
                
                if length(find(permAnovaP < fullAnovaP)) < 0.05 * numResampling % passing the test
                    tunedFull(ci) = 1;
                end
            end
            
            
            
            % lastly, touch-only model
            touchAnovaVal = zeros(sum(cellfun(@length, touchResponseTouch)),1);
            for ai = 1 : length(angles)
                touchAnovaVal(angleLengths(ai)+1:angleLengths(ai+1)) = touchResponseTouch{ai};
            end
            if ~isempty(find(isnan(touchAnovaVal)))
                error('nan values')
            end
            
            touchAnovaP = anova1(touchAnovaVal, groupAnova, 'off');
            touchAnovaPAll(ci) = touchAnovaP;

            tempH = cellfun(@(x) ttest(x), touchResponseTouch);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            if touchAnovaP < thresholdAnovaP && ~isempty(sigInd)
                % permutation test
                permAnovaP = zeros(numResampling,1);
                for ri = 1 : numResampling
                    tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                    permAnovaP(ri) = anova1(touchAnovaVal, tempG, 'off');
                end
                
                if length(find(permAnovaP < touchAnovaP)) < 0.05 * numResampling % passing the test
                    tunedTouch(ci) = 1;
                end
            end
            
            
        end
        
        save(sprintf('%s%03d\\%s', baseDir, mouse, savefn), '*All', 'beforePoleUpFrames', 'touchFrames', 'numTouchPreAnswer', ...
            'deFull', 'deTouchOnly', 'touchSpikes', 'touchFullModels', 'touchOnlyModels', 'depth', 'tuned', 'tunedSpk', 'tunedFull', 'tunedTouch', 'fullAnovaPAll', 'touchAnovaPAll')
        
    end
end






%% Making the summarized file
clear

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

saveFn = 'touchModel_angle_tuning_model_v9';

% naive first
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    loadfn = sprintf('%s%03d\\touchModel_angle_tuning_model_v9_JK%03dS%02d', baseDir, mouse, mouse, session);
    data = load(loadfn, '*All', 'tuned*','*AnovaPAll','de*','at*','*AngleAll');
    % remove touchSpikes, touchFullModels, and touchOnlyModels from data
%     data.cIDAll = cIDAll;
%     data.tuned = tuned;
%     data.tunedSpk = tunedSpk;
%     data.tunedFull = tunedFull;
%     data.tunedTouch = tunedTouch;
%     data.spkAnovaPAll = spkAnovaPAll;
%     data.fullAnovaPAll = fullAnovaPAll;
%     data.touchAnovaPAll = touchAnovaPAll;
%     data.deFull = deFull;
%     data.deTouchOnly = deTouchOnly;
%     data.depth = depth;
%     data.atSpikeAll = atSpikeAll;
%     data.atFullAll = atFullAll;
%     data.atTouchAll = atTouchAll;
%     data.spikeAngleAll = spikeAngleAll;
%     data.fullModelAngleAll = fullModelAngleAll;
%     data.touchOnlyAngleAll = touchOnlyAngleAll;
    naive(mi) = data;
end

% naive first
expertInd = find(cellfun(@length, sessions) == 2);
for mi = 1 : length(expertInd)
    mouse = mice(expertInd(mi));
    session = sessions{expertInd(mi)}(2);
    loadfn = sprintf('%s%03d\\touchModel_angle_tuning_model_v9_JK%03dS%02d', baseDir, mouse, mouse, session);
    data = load(loadfn, '*All', 'tuned*','*AnovaPAll','de*','at*','*AngleAll');
    expert(mi) = data;
end

save([baseDir, saveFn], 'naive', 'expert')



