% Impact on angle tuning.
% Compare angle-tuning in sensorimotor model.
% Remove all the other behavior categories.
% From full model and this whisker-only model, remove each feature and
% calculate angle-tuning correlation with spikes.
% Call this v9 because I found an error from previous calculations during
% Object Angle Tuning in S1 revision v9.
% This method is confirmed from d200618_closer_look_at_wkv_models.m
% 2020/06/18 JK

% fix: add tuned info during compiling 2020/06/21 JK

clear

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

angles = 45:15:135;

for mi = 1 : length(mice)
    for si = 1 : length(sessions{mi})
        mouse = mice(mi);
        session = sessions{mi}(si);
        
        fprintf('Processing JK%03d S%02d.\n', mouse, session)
        
        savefn = sprintf('wkv_angle_tuning_model_v9_JK%03dS%02d', mouse, session);
        
        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session), 'u')

        load(sprintf('%s%03d\\glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R01',baseDir, mouse, mouse, session), 'cIDAll', 'fitCoeffs', 'allPredictors', 'indPartial');
        coeffs = zeros(size(fitCoeffs,1),length(fitCoeffs{1}),10);
        for ri = 1 : 10
            glmfn = sprintf('%s%03d\\glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R%02d', baseDir, mouse, mouse, session, ri); % these are only from touch response cells
            load(glmfn, 'fitCoeffs')
            emptyCoeffsInd = find(cellfun(@isempty, fitCoeffs));
            firstNonemptyInd = find(1-cellfun(@isempty, fitCoeffs),1,'first');
            if ~isempty(emptyCoeffsInd)
                for emptyi = 1 : length(emptyCoeffsInd)
                    fitCoeffs{emptyCoeffsInd(emptyi)} = nan(size(fitCoeffs{firstNonemptyInd}));
                end
            end
            coeffs(:,:,ri) = cell2mat(fitCoeffs')';
        end
        meanCoeffsWKV = nanmean(coeffs,3);

        % convert wkv glm predictors back to each trial
        % there are 8 NaNs between trials. first trial starts with 4 NaNs, and the last ends with 4 NaNs.
        allPredictorsWKV = cell(8,1);
        for i = 1 : 8
            if ~isempty(allPredictors{i}) % it can be empty in case of JK027 S09 and S10
                naninds = find(isnan(allPredictors{i}(:,1)));
                nanindChanges = find(diff(naninds)>1);
                numTrials = length(nanindChanges);
                allPredictorsWKV{i} = cell(numTrials,1);
                for j = 1 : numTrials
                    allPredictorsWKV{i}{j} = allPredictors{i}(naninds(nanindChanges(j))+1:naninds(nanindChanges(j)+1)-1,:);
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
        whiskerFullModels = cell(length(cIDAll),1);
        whiskerOnlyModels = cell(length(cIDAll),1);
        deFull = zeros(length(cIDAll),1);
        deWhiskerOnly = zeros(length(cIDAll),1);
        parfor ci = 1 : length(cIDAll)
            pi = floor(cIDAll(ci)/1000);
            tempPredictors = cell2mat(allPredictorsWKV{pi});
            tempInput = [ones(size(tempPredictors,1),1), tempPredictors];
            whiskerFullModels{ci} = exp(tempInput*meanCoeffsWKV(ci,:)');

            tempInputWhiskerOnly = tempInput;
            tempInputWhiskerOnly(:,indPartial{1}(end)+2:end) = 0;
            whiskerOnlyModels{ci} = exp(tempInputWhiskerOnly*meanCoeffsWKV(ci,:)');

            spkTest = touchSpikes{ci};
            mu = nanmean(spkTest);
            nullLogLikelihood = nansum(log(poisspdf(spkTest,mu)));
            saturatedLogLikelihood = nansum(log(poisspdf(spkTest,spkTest)));
            fullLogLikelihood = nansum(log(poisspdf(spkTest',whiskerFullModels{ci})));
            deFull(ci) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            whiskerOnlyLogLikelihodd = nansum(log(poisspdf(spkTest',whiskerOnlyModels{ci})));
            deWhiskerOnly(ci) = (whiskerOnlyLogLikelihodd - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
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
        
        
        atCorrFull = cell(length(cIDAll), 1); %(:,1) full model, (:,2:13) removing each feature
        atCorrWhisker = cell(length(cIDAll), 1); %(:,1) whisker-only model, (:,2:13) removing each feature
        removeInds = cell(1,13);
        removeInds{1} = [];
        for i = 2 : 13
            removeInds{i} = (i-2)*3 + 1 + [1:3];
        end
        
        
        spikeAngleAll = cell(length(cIDAll), 1);
        fullModelAngleAll = cell(length(cIDAll), 1);
        whiskerOnlyAngleAll = cell(length(cIDAll), 1);
        atSpikeAll = cell(length(cIDAll),1);
        atFullAll = cell(length(cIDAll),1);
        atWhiskerAll = cell(length(cIDAll),1);
        parfor ci = 1 : length(cIDAll)
%         for ci = 1 : length(cIDAll)
%             fprintf('cell # %d\n', ci)
            cID = cIDAll(ci);
            pi = floor(cID/1000);
            tempInd = find(cellfun(@(x) ismember(pi, x.planes), u.trials), 1, 'first');
            tempcIDlist = u.trials{tempInd}.neuindSession;
            cind = find(tempcIDlist == cID);
            spikeAngle = cell(1, length(angles));
            fullModelAngle = cell(13, length(angles));
            whiskerOnlyAngle = cell(13, length(angles));

            touchFramesAngle = cell(1, length(angles));
            baselineFramesAngle = cell(1, length(angles));
            numTouchAngle = cell(1, length(angles));
            for ai = 1 : length(angles)
                tempTouchTrialInd = angleTouchTrialInds{pi}{ai}; % index of planeTouchTrialsInd{pi}
                tempTrialInd = planeTouchTrialsInd{pi}(tempTouchTrialInd); % index of u
                spikeAngle{ai} = cell(length(tempTrialInd),1);
                % feature removal
                for fi = 1 : 13
                    fullModelAngle{fi,ai} = cell(length(tempTrialInd),1);
                    whiskerOnlyAngle{fi,ai} = cell(length(tempTrialInd),1);
                end

                touchFramesAngle{ai} = cell(length(tempTrialInd),1);
                baselineFramesAngle{ai} = cell(length(tempTrialInd),1);
                numTouchAngle{ai} = cell(length(tempTrialInd),1);
                for ti = 1 : length(tempTrialInd)
                    spikeAngle{ai}{ti} = u.trials{tempTrialInd(ti)}.spk(cind,:);

                    tempPredictorInd = find(planeTrialsInd{pi} == tempTrialInd(ti)); % index matched with the order in allPredictorsWKV
                    tempPredictor = allPredictorsWKV{pi}{tempPredictorInd};
                    tempInput = [ones(size(tempPredictor,1),1), tempPredictor];
                    for fi = 1 : 13
                        tempCoeff = meanCoeffsWKV(ci,:);
                        tempCoeff(1,removeInds{fi}) = 0;
                        fullModelAngle{fi,ai}{ti} = exp(tempInput * tempCoeff');

                        tempCoeff(1,indPartial{1}(end)+2:end) = 0;
                        whiskerOnlyAngle{fi,ai}{ti} = exp(tempInput * tempCoeff');
                    end

                    touchFramesAngle{ai}{ti} = touchFrames{pi}{tempTouchTrialInd(ti)};
                    baselineFramesAngle{ai}{ti} = beforePoleUpFrames{pi}{tempTouchTrialInd(ti)};
                    numTouchAngle{ai}{ti} = numTouchPreAnswer{pi}(tempTouchTrialInd(ti));
                end
            end
            
            % trial-by-trial model reconstruction
            touchResponseSpike = cell(1, length(angles));
            touchResponseFull = cell(13, length(angles));
            touchResponseWhisker = cell(13, length(angles));

            for ai = 1 : length(angles)
                touchResponseSpike{ai} = zeros(length(spikeAngle{ai}),1);
                for fi = 1 : 13
                    touchResponseFull{fi,ai} = zeros(length(spikeAngle{ai}),1);
                    touchResponseWhisker{fi,ai} = zeros(length(spikeAngle{ai}),1);
                end
                for ti = 1 : length(spikeAngle{ai})
                    touchResponseSpike{ai}(ti) = (  sum( spikeAngle{ai}{ti}(touchFramesAngle{ai}{ti}) ) - mean( spikeAngle{ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti};
                    for fi = 1 : 13
                        touchResponseFull{fi,ai}(ti) = (  sum( fullModelAngle{fi,ai}{ti}(touchFramesAngle{ai}{ti}) ) - nanmean( fullModelAngle{fi,ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti}; % nanmean because models start with NaN values
                        touchResponseWhisker{fi,ai}(ti) = (  sum( whiskerOnlyAngle{fi,ai}{ti}(touchFramesAngle{ai}{ti}) ) - nanmean( whiskerOnlyAngle{fi,ai}{ti}(baselineFramesAngle{ai}{ti}) )  ) / numTouchAngle{ai}{ti};
                    end
                end
            end

            % calculate angle-tuning curve correlation
            atSpike = cellfun(@mean, touchResponseSpike);
            atFull = cellfun(@mean, touchResponseFull);
            atWhisker = cellfun(@mean, touchResponseWhisker);
            
            tempCorrFull = zeros(1,13);
            tempCorrWhisker = zeros(1,13);
            for fi = 1 : 13
                tempCorrFull(1,fi) = corr(atSpike', atFull(fi,:)');
                tempCorrWhisker(1,fi) = corr(atSpike', atWhisker(fi,:)');
            end
            
            atSpikeAll{ci} = atSpike;
            atFullAll{ci} = atFull;
            atWhiskerAll{ci} = atWhisker;
            
            atCorrFull{ci} = tempCorrFull;
            atCorrWhisker{ci} = tempCorrWhisker;
            
            spikeAngleAll{ci} = touchResponseSpike;
            fullModelAngleAll{ci} = touchResponseFull;
            whiskerOnlyAngleAll{ci} = touchResponseWhisker;
        end
        
        save(sprintf('%s%03d\\%s', baseDir, mouse, savefn), '*All', 'atCorr*', 'beforePoleUpFrames', 'touchFrames', 'numTouchPreAnswer', 'deFull', 'deWhiskerOnly')
        
    end
end






%% Making the summarized file
clear

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

saveFn = 'wkv_angle_tuning_model_v9';

tune = load([baseDir, 'angle_tuning_summary_preAnswer_perTouch_NC']);
% naive first
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    loadfn = sprintf('%s%03d\\wkv_angle_tuning_model_v9_JK%03dS%02d', baseDir, mouse, mouse, session);
    dat = load(loadfn);
    if isempty(find(tune.naive(mi).touchID - dat.cIDAll))
        dat.tuned = tune.naive(mi).tuned;
        dat.depth = tune.naive(mi).depth;
    end
    naive(mi) = dat;
end

% naive first
expertInd = find(cellfun(@length, sessions) == 2);
for mi = 1 : length(expertInd)
    mouse = mice(expertInd(mi));
    session = sessions{expertInd(mi)}(2);
    loadfn = sprintf('%s%03d\\wkv_angle_tuning_model_v9_JK%03dS%02d', baseDir, mouse, mouse, session);
    dat = load(loadfn);
    if isempty(find(tune.expert(mi).touchID - dat.cIDAll))
        dat.tuned = tune.expert(mi).tuned;
        dat.depth = tune.expert(mi).depth;
    end
    expert(mi) = dat;
end

save([baseDir, saveFn], 'naive', 'expert')



