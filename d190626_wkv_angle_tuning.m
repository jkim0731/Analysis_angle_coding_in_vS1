% Making spike matrix as in u, but from a model with removing selected whisker kinematic variables.
% This is to show which variables are important in contructing angle tuning.
% Only consider previously angle-tuned cells.
% In different settings, run angle tuning, and record if the cell is deemed tuned or not.

% First, build spike matrix from the full model.
% Then, build a full model from touch glm.
% After that, build a full model from wkv glm.
% List of variables to remove: test set - dKv, slide distance, touch quality control - touch num, and arc length, negative control - dKh and dPhi.
% Combination of variables to remove: dKv + slide distance (test). Touch num + arc length (control). dKh + dPhi
% List of variables to include: dKv, slide distance, dKv + slide distance, touch num + arc length, dKh + dPhi

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
sessions = {[],[3,10],[3,21],[1,17],[7],[],[1,23],[],[],[],[],[]};
sessions = {[],[3],[],[1],[7],[],[],[],[],[],[],[]};
cd(baseDir)
load('cellFunctionRidgeDE010_JK027_S10.mat')

naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];
angles = 45:15:135;
numResampling = 10000;
thresholdAnovaP = 0.05;
thresholdPermutation = 0.05;

% naiveModelTune = struct;
% expertModelTune = struct;
for mi = 1 : length(mice)
% for mi = 2
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    for si = 1 : length(sessions{mi})
%     for si = 1
        session = sessions{mi}(si);        
        savefn = sprintf('angle_tuning_model_JK%03dS%02d', mouse, session);
        
        % load uber
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load(ufn)
        
        % load wkv glm results and get coefficients (for all touch cells)
        glmfn = sprintf('glmWhiskerTouchVariablesONLYlasso_JK%03dS%02d_R10',mouse, session); % these are only from touch response cells
        load(glmfn, 'cIDAll', 'fitCoeffs', 'allPredictors')
        ap1 = allPredictors;
        cIDAllWKV = cIDAll;        
        coeffs = zeros(length(fitCoeffs),length(fitCoeffs{1}),10);
        for gi = 1 : 10
            glmfn = sprintf('glmWhiskerTouchVariablesONLYlasso_JK%03dS%02d_R%02d',mouse, session, gi); % these are only from touch response cells
            load(glmfn, 'fitCoeffs')
            emptyCoeffsInd = find(cellfun(@isempty, fitCoeffs));
            firstNonemptyInd = find(1-cellfun(@isempty, fitCoeffs),1,'first');
            if ~isempty(emptyCoeffsInd)
                for emptyi = 1 : length(emptyCoeffsInd)                    
                    fitCoeffs{emptyCoeffsInd(emptyi)} = nan(size(fitCoeffs{firstNonemptyInd}));
                end
            end
            coeffs(:,:,gi) = cell2mat(fitCoeffs')';
        end
        meanCoeffsWKV = nanmean(coeffs,3);
        % convert wkv glm predictors back to each trial
        % there are 8 NaNs between trials. first trial starts with 4 NaNs, and the last ends with 4 NaNs.
        allPredictorsWKV= cell(8,1);
%         for i = 1 : 8
%             if ~isempty(allPredictors{i}) % it can be empty in case of JK027 S09 and S10
%                 
%                 naninds = find(isnan(allPredictors{i}(:,1)));
%                 nanindChanges = find(diff(naninds)>1);
%                 numTrials = length(nanindChanges);
%                 allPredictorsWKV{i} = cell(numTrials,1);
%                 for j = 1 : numTrials
%                     allPredictorsWKV{i}{j} = allPredictors{i}(naninds(nanindChanges(j))+1:naninds(nanindChanges(j)+1)-1,:);
%                 end
%             end
%         end
        
        % load touch glm results and get coefficients (for all active cells)
        glmfn = sprintf('glmResponseType_JK%03dS%02d_m44_R10',mouse, session); % these are only from all active cells
        load(glmfn, 'cIDAll', 'fitCoeffs', 'allPredictors')
        ap2 = allPredictors;
        cIDAllTouch = cIDAll;
        coeffs = zeros(length(fitCoeffs),length(fitCoeffs{1}),10);
        for gi = 1 : 10
            glmfn = sprintf('glmResponseType_JK%03dS%02d_m44_R%02d',mouse, session, gi); % these are only from touch response cells
            load(glmfn, 'fitCoeffs')
            emptyCoeffsInd = find(cellfun(@isempty, fitCoeffs));
            firstNonemptyInd = find(1-cellfun(@isempty, fitCoeffs),1,'first');
            if ~isempty(emptyCoeffsInd)
                for emptyi = 1 : length(emptyCoeffsInd)                    
                    fitCoeffs{emptyCoeffsInd(emptyi)} = nan(size(fitCoeffs{firstNonemptyInd}));
                end
            end
            coeffs(:,:,gi) = cell2mat(fitCoeffs')';
        end
        meanCoeffsTouch = nanmean(coeffs,3);
        % convert touch glm predictors back to each trial
        allPredictorsTouch = cell(8,1);
        for i = 1 : 8
            if ~isempty(allPredictors{i}) % it can be empty in case of JK027 S09 and S10
                naninds = find(isnan(allPredictors{i}(:,1)));
                nanindChanges = find(diff(naninds)>1);
                numTrials = length(nanindChanges);
                allPredictorsTouch{i} = cell(numTrials,1);
                allPredictorsWKV{i} = cell(numTrials,1);
                for j = 1 : numTrials
                    allPredictorsTouch{i}{j} = ap2{i}(naninds(nanindChanges(j))+1:naninds(nanindChanges(j)+1)-1,:);
                    allPredictorsWKV{i}{j} = ap1{i}(naninds(nanindChanges(j))+1:naninds(nanindChanges(j)+1)-1,:);
                end
            end
        end
        
        
        
        % find glm results
        if si == 2
            glmi = find(expertMi == mi);
            glm = expert(glmi);
        else
            glm = naive(mi);
        end        
        touchID = glm.touchID;
        if size(touchID,1) < size(touchID,2)
            touchID = touchID';
        end
        
        % dKv: [5:7]
        % slide distance: [14:16]
        % touch num: [35:37]
        % arc length: [32:34]
        % dKh: [2:4]
        % dPhi: [11:13]
        
%         total 17 sets for testing angle tuning
%         (1) raw dF
%         (2) full model with touch angle
%         (3) full model with wkv
%         (4) from wkv, remove dKv (remove 5:7)
%         (5) from wkv, remove slide distance (remove 14:16)
%         (6) from wkv, remove touch num (remove 35:37)
%         (7) from wkv, remove arc length (remove 32:34)
%         (8) from wkv, remove dKh (remove 2:4)
%         (9) from wkv, remove dPhi (remove 11:13)
%         (10) from wkv, remove dKv + slide distance (remove [5:7, 14:16])
%         (11) from wkv, remove touch num + arc length (remove 32:37)
%         (12) from wkv, remove dKh + dPhi (remove [2:4, 11:13])
%         (13) from wkv, use only dKv (use [1, 5:7])
%         (14) from wkv, use only slide distance (use [1, 14:16])
%         (15) from wkv, use only dKv + slide distance (use [1, 5:7, 14:16])
%         (16) from wkv, use only touch num + arc length (use [1, 32:37])
%         (17) from wkv, use only dKh + dPhi (use [1:4, 11:13])
%         
        
%         First, align them all (model and dF), and align with baseline frames and touch frames for each corresponding trial
%         Each in different planes
        
        % making templates
        % find predecision touch trials
        decisionTime = cell(length(u.trials),1);
        for di = 1 : length(decisionTime)
            if isempty(u.trials{di}.answerLickTime)
                decisionTime{di} = u.trials{di}.poleDownOnsetTime;
            else
                decisionTime{di} = u.trials{di}.answerLickTime;
            end
        end
        tempTouchTrialInd = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        pdTouchInd = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y, u.trials(tempTouchTrialInd), decisionTime(tempTouchTrialInd)));
        touchTrialInd = tempTouchTrialInd(pdTouchInd); % index of u
        numPlane = length(u.mimg);
        planeTrialsInd = cell(numPlane,1);
        planeTrialsNum = cell(numPlane,1);
        poleUpFrames = cell(numPlane,1);
        beforePoleUpFrames = cell(numPlane,1);
        touchFrames = cell(numPlane,1);
        angleTrialInds = cell(numPlane,1);
        for pi = 1 : numPlane
            planeTrialsInd{pi} = intersect(find(cellfun(@(x) ismember(pi, x.planes), u.trials)), touchTrialInd); % index of u
            tempInd = find(u.trials{planeTrialsInd{pi}(1)}.planes == pi);
            poleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} >= x.poleUpTime(1) & x.tpmTime{tempInd} <= x.poleUpTime(end)), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            beforePoleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} < x.poleUpOnsetTime), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            touchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            for ti = 1 : length(planeTrialsInd{pi})
                tempTrial = u.trials{planeTrialsInd{pi}(ti)};
                if isempty(tempTrial.answerLickTime)
                    tempDecisionTime = tempTrial.poleDownOnsetTime;
                else
                    tempDecisionTime = tempTrial.answerLickTime;
                end
                preDecisionInd = find(cellfun(@(x) tempTrial.whiskerTime(x(1)) < tempDecisionTime, tempTrial.protractionTouchChunksByWhisking));
                tempFrames = cell(1, length(preDecisionInd));
                for ptci = 1 : length(tempFrames)
                    tempFrames{ptci} = [0:1] + find(tempTrial.tpmTime{tempInd} >= tempTrial.whiskerTime(tempTrial.protractionTouchChunksByWhisking{ptci}(1)), 1, 'first');
                end
                touchFrames{pi}{ti} = unique(cell2mat(tempFrames));
            end
            angleTrialInds{pi} = cell(length(angles),1); % index of planeTrialsInd{pi}
            for ai = 1 : length(angles)
                angleTrialInds{pi}{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(planeTrialsInd{pi})));
            end
        end
        
        
        % Each allPredictors should have corresponding trial index of uber
        % array, considering existence of predecision touches and angles
        % It should be the same between WKV and Touch GLM
        allPredictorsTouchAngleInds = cell(8,1); % index of predictors{i}
        for i = 1 : 8
            tempInd = find(cellfun(@(x) ismember(i, x.planes), u.trials));
            allPredictorsTouchAngleInds{i} = cell(length(angles),1);
            for ai = 1 : length(angles)
                allPredictorsTouchAngleInds{i}{ai} = find(ismember( tempInd, intersect(intersect(tempInd, touchTrialInd), ...
                    find(cellfun(@(x) x.angle == angles(ai), u.trials))) ));
            end
        end
        
        
        % settings for variables to be saved later
%         spkValAllCell = cell(length(touchID),17);
        spkValAllCell = cell(length(touchID),17);
        anovaPAllCell = zeros(length(touchID),17);
        tunedAllCell = zeros(length(touchID),17);
        tuneAngleAllCell = zeros(length(touchID),17);
        
        % angle tuning in each cell
        parfor ci = 1:length(touchID)
%         for ci = 1:length(touchID)
%         for ci = 27
            fprintf('Processing JK%03d S%02d touch cell %d / %d\n', mouse, session, ci, length(touchID))
            cellNum = touchID(ci);
            plane = floor(cellNum/1000);
            trialInds = planeTrialsInd{plane}; % index of u
            calciumPoleUpFrames = poleUpFrames{plane};
            spkTouchFrames = touchFrames{plane};
            baselineFrames = beforePoleUpFrames{plane};
            angleInds = angleTrialInds{plane}; % index of trialInds            
            modelTouchAngleInds = allPredictorsTouchAngleInds{plane}; % index of allPredictorsTouch & allPredictorsWKV
            
            spkValAll = cell(1,17);
            
            % all spikes
            cind = find(u.trials{trialInds(1)}.neuindSession == cellNum);
            tempSpk = cellfun(@(x) x.spk(cind,:), u.trials(trialInds), 'uniformoutput', false);
            spkValAll{1} = cell(length(angles),1);
            for ai = 1 : length(angles)
                trialAngleInd = angleInds{ai};               
                spkValAll{1}{ai} = zeros(length(trialAngleInd),1);
                for ti = 1 : length(trialAngleInd)
                    tempInd = trialAngleInd(ti);
                    spkValAll{1}{ai}(ti) = mean(tempSpk{tempInd}(spkTouchFrames{tempInd})) - mean(tempSpk{tempInd}(baselineFrames{tempInd}));
                end
            end
            
            % from full touch glm (#2)
            coeffInd = find(cIDAllTouch == cellNum);
            spkValAll{2} = cell(length(angles),1);
            for ai = 1 : length(angles)
                trialAngleInd = modelTouchAngleInds{ai};
                spkValAll{2}{ai} = zeros(length(trialAngleInd),1);
                for ti = 1 : length(trialAngleInd)
                    tempInd = trialAngleInd(ti); % index of allPredictorsTouch
                    tempPredictor = allPredictorsTouch{plane}{tempInd};
                    tempLength = size(tempPredictor,1);
                    tempInput = [ones(tempLength,1), tempPredictor];
                    tempCoeff = meanCoeffsTouch(coeffInd,:);
                    model = exp(tempInput*tempCoeff');
                    
                    tempUInd = find(cellfun(@(x) ismember(plane, x.planes), u.trials));
                    uInd = tempUInd(tempInd);
                    matchingInd = find(trialInds == uInd);
                    spkValAll{2}{ai}(ti) = mean(model(spkTouchFrames{matchingInd})) - nanmean(model(baselineFrames{matchingInd}));
                end
            end
            
            
            % from full wkv glm (#3 ~ #17)
            coeffInd = find(cIDAllWKV == cellNum);
            %         (4) from wkv, remove dKv (remove 5:7)
            %         (5) from wkv, remove slide distance (remove 14:16)
            %         (6) from wkv, remove touch num (remove 35:37)
            %         (7) from wkv, remove arc length (remove 32:34)
            %         (8) from wkv, remove dKh (remove 2:4)
            %         (9) from wkv, remove dPhi (remove 11:13)
            %         (10) from wkv, remove dKv + slide distance (remove [5:7, 14:16])
            %         (11) from wkv, remove touch num + arc length (remove 32:37)
            %         (12) from wkv, remove dKh + dPhi (remove [2:4, 11:13])
            %         (13) from wkv, use only dKv (use [1, 5:7])
            %         (14) from wkv, use only slide distance (use [1, 14:16])
            %         (15) from wkv, use only dKv + slide distance (use [1, 5:7, 14:16])
            %         (16) from wkv, use only touch num + arc length (use [1, 32:37])
            %         (17) from wkv, use only dKh + dPhi (use [1:4, 11:13])
            inds = cell(15,1);
            coeffLength = size(meanCoeffsWKV,2);
            inds{1} = 1:coeffLength;
            inds{2} = setdiff(1:coeffLength, 5:7);
            inds{3} = setdiff(1:coeffLength, 14:16);
            inds{4} = setdiff(1:coeffLength, 35:37);
            inds{5} = setdiff(1:coeffLength, 32:34);
            inds{6} = setdiff(1:coeffLength, 2:4);
            inds{7} = setdiff(1:coeffLength, 11:13);
            inds{8} = setdiff(1:coeffLength, [5:7, 14:16]);
            inds{9} = setdiff(1:coeffLength, 32:37);
            inds{10} = setdiff(1:coeffLength, [2:4, 11:13]);
            inds{11} = [1,5:7];
            inds{12} = [1,14:16];
            inds{13} = [1,5:7, 14:16];
            inds{14} = [1, 32:37];
            inds{15} = [1:4, 11:13];
            for i = 3 : 17
                spkValAll{i} = cell(length(angles),1);
            end
            
            anovaPcell = zeros(1,17);
            tunedCell = zeros(1,17);
            tuneAngleCell = nan(1,17);
            for ai = 1 : length(angles)
                trialAngleInd = modelTouchAngleInds{ai};
                for i = 3 : 17
                    spkValAll{i}{ai} = zeros(length(trialAngleInd),1);
                end
                for ti = 1 : length(trialAngleInd)
                    tempInd = trialAngleInd(ti); % index of allPredictorsTouch & WKV
                    tempPredictor = allPredictorsWKV{plane}{tempInd};
                    tempLength = size(tempPredictor,1);
                    tempInput = [ones(tempLength,1), tempPredictor];
                    tempCoeff = meanCoeffsWKV(coeffInd,:);

                    tempUInd = find(cellfun(@(x) ismember(plane, x.planes), u.trials));
                    uInd = tempUInd(tempInd);
                    matchingInd = find(trialInds == uInd);
                    
                    for i = 1 : 15
                        model = exp(tempInput(:,inds{i})*tempCoeff(inds{i})');
                        spkValAll{i+2}{ai}(ti) = nanmean(model(spkTouchFrames{matchingInd})) - nanmean(model(baselineFrames{matchingInd}));
                    end
                end
            end
            
            % ANOVA in each configuration            
            anovaVal = cell2mat(spkValAll{1});
            groupAnova = zeros(size(anovaVal));
            angleLengths = [0;cumsum(cellfun(@length, spkValAll{1}))];
            for ai = 1 : length(angles)
                groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
            end
            for i = 1 : size(spkValAll,2)
                anovaVal = cell2mat(spkValAll{i});
                anovaP = anova1(anovaVal, groupAnova, 'off');
                anovaPcell(i) = anovaP;
                meanVal = cellfun(@nanmean, spkValAll{i});
                tempH = cellfun(@(x) ttest(x), spkValAll{i});
                tempH(isnan(tempH)) = deal(0);
                sigInd = find(tempH); % significant indices
                if anovaP < thresholdAnovaP && ~isempty(sigInd)
                    % permutation test
                    permAnovaP = zeros(numResampling,1);
                    for ri = 1 : numResampling
                        tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                        permAnovaP(ri) = anova1(anovaVal, tempG, 'off');                        
                    end
                    if length(find(permAnovaP < anovaP)) < thresholdPermutation * numResampling % passed permutation test
                        tunedCell(i) = 1;                        
                        [~, maxind] = max(abs(meanVal(sigInd)));
                        tunedAngleInd = sigInd(maxind);
                        tuneAngleCell(i) = angles(tunedAngleInd);                    
                    end
                end
            end
            spkValAllCell(ci,:) = spkValAll;
            anovaPAllCell(ci,:) = anovaPcell;
            tunedAllCell(ci,:) = tunedCell;
            tuneAngleAllCell(ci,:) = tuneAngleCell;            
        end
        
        save(savefn, '*AllCell')
%         if si == 2
%             glmi = find(expertMi == mi);
%             expertModelTune(glmi).spkValAllCell = spkValAllCell;
%             expertModelTune(glmi).anovaPAllCell = spkValAllCell;
%             expertModelTune(glmi).angleTunedAllCell = spkValAllCell;
%         else
%             naiveModelTune(glmi).spkValAllCell = spkValAllCell;
%             naiveModelTune(glmi).anovaPAllCell = spkValAllCell;
%             naiveModelTune(glmi).angleTunedAllCell = spkValAllCell;
%         end
    end
end
        
% save(sprintf('%sangle_tuning_model',baseDir), '*ModelTune')
        
% %%
% ts1 = cellfun(@mean, spkValAll{1});
% ts1 = (ts1 - min(ts1)) / (max(ts1) - min(ts1));
% ts2 = cellfun(@mean, spkValAll{2});
% ts2 = (ts2 - min(ts2)) / (max(ts2) - min(ts2));
% ts3 = cellfun(@nanmean, spkValAll{3});
% ts3 = (ts3 - min(ts3)) / (max(ts3) - min(ts3));
% 
% figure, 
% subplot(121), hold on, plot(cellfun(@mean, spkValAll{1})), plot(cellfun(@mean, spkValAll{2})), plot(cellfun(@nanmean, spkValAll{3}))
% subplot(122), hold on, plot(ts1), plot(ts2), plot(ts3)
% 
% legend({'spikes', 'touch model', 'wkv model'})
% 
% %% test. Compare DE between WKV and touch GLM (both lasso)
% 
% 
% 
% 
