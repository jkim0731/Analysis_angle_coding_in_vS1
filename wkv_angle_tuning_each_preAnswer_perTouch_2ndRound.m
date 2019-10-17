% Second round to wkv model angle tuning reconstruction.
% 3 features, dPhi, dKv, slide distance, selected from both removing and singe correlation
% Including different combinations, based on the observation in the first pass 2019/10/04
% Also considering other variables than the whisker variables
% 2019/10/05 JK
%%

% Making spike matrix as in u, but from a model with removing selected whisker kinematic variables.
% This is to show which variables are important in contructing angle tuning.
% Only consider previously angle-tuned cells.
% In different settings, run angle tuning, and record if the cell is deemed tuned or not.

% First, build spike matrix from the full model.
% Then, build a full model from touch glm.
% After that, build a full model from wkv glm.

% 2019/09/27 JK
% Copied from wkv_angl_tuning_each
% Major change: calculate angle tuning based on per-touch spikes,
% same as in angle_tuning_preAnswer_perTouch_spkOnly
% Added some properties to be calculated. (type, modulation, sharpness)

%% Basic settings
clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
% sessions = {[],[3,10],[3,21],[1,17],[7],[],[1,23],[],[],[],[],[]};
% sessions = {[],[3],[],[1],[7],[],[],[],[],[],[],[]};
cd(baseDir)
load('cellFunctionLasso_NC.mat')

naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];
angles = 45:15:135;
numResampling = 10000;
thresholdAnovaP = 0.05;
thresholdPermutation = 0.05;
anovactype = 'hsd';
thresholdCategory = 0.05;

%% target feature setting
%         total 25 sets for testing angle tuning
%             1: without any whisker variable
%             2: maxDtheta + other variables
%             3: maxDphi + other variables
%             4: maxDkH + other variables
%             5: maxDkV + other variables
%             6: max(Slide distance) + other variables
%             7: max(duration) + other variables
%             8: thetaAtTouch + other variables
%             9: phiAtTouch + other variables
%             10: kHAtTouch + other variables
%             11: kVAtTouch + other variables
%             12: arc length + other variables
%             13: touch count + other variables
%             14: -(maxDphi + maxDkV)
%             15: -(maxDphi + max(Slide distance))
%             16: -(maxDkV + max(Slide distance))
%             17: -(maxDphi + maxDkV + max(Slide distance))
%             18: dPhi + dKv
%             19: dPhi + slide distance
%             20: dKv + slide distance
%             21: dPhi + dKv + slide distance
%             22: (dPhi + dKv) + other
%             23: (dPhi + slide distance) + other
%             24: (dKv + slide distance) + other
%             25: (maxDphi + maxDkV + max(Slide distance)) + other
numFeature = 25;
featureNames = {'noWhiskerVariable', ...
    'maxDtheta+other', ...
    'maxDphi+other', ...
    'maxDkH+other', ...
    'maxDkV+other', ...
    'max(SlideDistance)+other', ...
    'max(duration)+other', ...
    'thetaAtTouch+other', ...
    'phiAtTouch+other', ...
    'kHAtTouch+other', ...
    'kVAtTouch+other', ...
    'arcLength+other', ...
    'touchCount+other', ...
    'all-(maxDphi+maxDkV)', ...
    'all-(maxDphi+max(SlideDistance)', ...
    'all-(maxDphi+max(SlideDistance)', ...
    'all-(maxDphi+maxDkV+max(SlideDistance))', ...
    'dPhi + dKv', ...
    'dPhi + slide distance',...
    'dKv + slide distance',...
    'dPhi + dKv + slide distance',...
    '(dPhi + dKv) + other',...
    '(dPhi + slide distance) + other',...
    '(dKv + slide distance) + other',...
    '(maxDphi + maxDkV + max(Slide distance)) + other'};


% naiveModelTune = struct;
% expertModelTune = struct;
for mi = 1 : length(mice)
% for emi = 1:length(expertMi)
%     mi = expertMi(emi);
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    for si = 1 : length(sessions{mi})
%     for si = 2
        session = sessions{mi}(si);        
        savefn = sprintf('angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK%03dS%02d_2ndRound', mouse, session);
        
        % load uber
        ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
        load(ufn)
        
        % load wkv glm results and get coefficients (for all touch cells)
        glmfn = sprintf('glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R10',mouse, session); % these are only from touch response cells
        load(glmfn, 'cIDAll', 'fitCoeffs', 'allPredictors')
        ap1 = allPredictors;
        cIDAllWKV = cIDAll;        
        coeffs = zeros(length(fitCoeffs),length(fitCoeffs{1}),10);
        for gi = 1 : 10
            glmfn = sprintf('glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R%02d',mouse, session, gi); % these are only from touch response cells
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
        allPredictorsWKV = cell(8,1);
        
        for i = 1 : 8
            if ~isempty(allPredictors{i}) % it can be empty in case of JK027 S09 and S10
                naninds = find(isnan(allPredictors{i}(:,1)));
                nanindChanges = find(diff(naninds)>1);
                numTrials = length(nanindChanges);
                allPredictorsWKV{i} = cell(numTrials,1);
                for j = 1 : numTrials
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
        

        
        % making templates
        % find preAnswer touch trials
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
        planeTrialsInd = cell(numPlane,1);
        planeTrialsNum = cell(numPlane,1);
        poleUpFrames = cell(numPlane,1);
        beforePoleUpFrames = cell(numPlane,1);
        touchFrames = cell(numPlane,1);
        numTouchPreAnswer = cell(numPlane,1);
        angleTrialInds = cell(numPlane,1);
        for pi = 1 : numPlane
            planeTrialsInd{pi} = intersect(find(cellfun(@(x) ismember(pi, x.planes), u.trials)), touchTrialInd); % index of u
            tempInd = find(u.trials{planeTrialsInd{pi}(1)}.planes == pi);
            poleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} >= x.poleUpTime(1) & x.tpmTime{tempInd} <= x.poleUpTime(end)), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            beforePoleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} < x.poleUpOnsetTime), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            touchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            numTouchPreAnswer{pi} = zeros(length(planeTrialsInd{pi}),1);
            for ti = 1 : length(planeTrialsInd{pi})
                tempTrial = u.trials{planeTrialsInd{pi}(ti)};
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
        spkValAllCell = cell(length(touchID),numFeature);
        anovaPAllCell = zeros(length(touchID),numFeature);
        tunedAllCell = zeros(length(touchID),numFeature);
        tuneAngleAllCell = zeros(length(touchID),numFeature);
        tuneModulationAllCell = zeros(length(touchID), numFeature);
        tuneSharpnessAllCell = zeros(length(touchID),numFeature);
        unimodalSingleAllCell = zeros(length(touchID),numFeature);
        unimodalBroadAllCell = zeros(length(touchID),numFeature);
        multimodalAllCell = zeros(length(touchID),numFeature);
            
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
            numTouch = numTouchPreAnswer{plane};
            
            spkValAll = cell(1,numFeature);
            
            % whiskerTouchMat = [maxDthetaMat, maxDphiMat, maxDkappaHMat, maxDkappaVMat, maxSlideDistanceMat, maxDurationMat, ...    
%                             thetaAtTouchMat, phiAtTouchMat, kappaHAtTouchMat, kappaVAtTouchMat, arcLengthAtTouchMat, touchCountMat];
            
%             1: without any whisker variable
%             2: maxDtheta + other variables
%             3: maxDphi + other variables
%             4: maxDkH + other variables
%             5: maxDkV + other variables
%             6: max(Slide distance) + other variables
%             7: max(duration) + other variables
%             8: thetaAtTouch + other variables
%             9: phiAtTouch + other variables
%             10: kHAtTouch + other variables
%             11: kVAtTouch + other variables
%             12: arc length + other variables
%             13: touch count + other variables
%             14: -(maxDphi + maxDkV)
%             15: -(maxDphi + max(Slide distance))
%             16: -(maxDkV + max(Slide distance))
%             17: -(maxDphi + maxDkV + max(Slide distance))
%             18: dPhi + dKv
%             19: dPhi + slide distance
%             20: dKv + slide distance
%             21: dPhi + dKv + slide distance
%             22: (dPhi + dKv) + other
%             23: (dPhi + slide distance) + other
%             24: (dKv + slide distance) + other
%             25: (maxDphi + maxDkV + max(Slide distance)) + other

            coeffInd = find(cIDAllWKV == cellNum);
            inds = cell(numFeature,1);
            coeffLength = size(meanCoeffsWKV,2);
            
            inds{1} = setdiff(1:coeffLength, 2:37); % without any whisker variable
            inds{2} = [1,2:4, 38:coeffLength]; % maxDtheta + other variables
            inds{3} = [1,5:7, 38:coeffLength]; % maxDphi + other variables
            inds{4} = [1,8:10, 38:coeffLength]; % maxDkH + other variables
            inds{5} = [1,11:13, 38:coeffLength]; % maxDkV + other variables
            inds{6} = [1,14:16, 38:coeffLength]; % max(Slide distance) + other variables
            inds{7} = [1,17:19, 38:coeffLength]; % max(duration) + other variables
            inds{8} = [1,20:22, 38:coeffLength]; % thetaAtTouch + other variables
            inds{9} = [1,23:25, 38:coeffLength]; % phiAtTouch + other variables
            inds{10} = [1,26:28, 38:coeffLength]; % kHAtTouch + other variables
            inds{11} = [1,29:31, 38:coeffLength]; % kVAtTouch + other variables
            inds{12} = [1,32:34, 38:coeffLength]; % arc length + other variables
            inds{13} = [1,35:37, 38:coeffLength]; % touch count + other variables
            
            inds{14} = setdiff(1:coeffLength, [5:7, 11:13]); % -(dPhi + dKv)
            inds{15} = setdiff(1:coeffLength, [5:7, 14:16]); % -(dPhi + slide distance)
            inds{16} = setdiff(1:coeffLength, [11:13, 14:16]); % -(dKv + slide distance)
            inds{17} = setdiff(1:coeffLength, [5:7, 11:13, 14:16]); % - (maxDphi + maxDkV + max(Slide distance))
            
            inds{18} = [1, 5:7, 11:13]; % (maxDphi + maxDkV)
            inds{19} = [1, 5:7, 14:16]; % (maxDphi + max(Slide distance))
            inds{20} = [1, 11:13, 14:16]; % (maxDkV + max(Slide distance))
            inds{21} = [1, 5:7, 11:13, 14:16]; % (maxDphi + maxDkV + max(Slide distance))
            
            inds{22} = [1, 5:7, 11:13, 38:coeffLength]; % (dPhi + dKv) + other
            inds{23} = [1, 5:7, 14:16, 38:coeffLength]; % (dPhi + max(Slide distance)) + other
            inds{24} = [1, 11:13, 14:16, 38:coeffLength]; % (dKv + max(Slide distance)) + other
            inds{25} = [1, 5:7, 11:13, 14:16, 38:coeffLength]; % (maxDphi + maxDkV + max(Slide distance)) + other
            
            for i = 1 : numFeature
                spkValAll{i} = cell(length(angles),1);
            end
            
            for ai = 1 : length(angles)
                trialAngleInd = modelTouchAngleInds{ai};
                for i = 1 : numFeature
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
                    
                    for i = 1 : numFeature
                        model = exp(tempInput(:,inds{i})*tempCoeff(inds{i})');
                        spkValAll{i}{ai}(ti) = nansum(model(spkTouchFrames{matchingInd}) - nanmean(model(baselineFrames{matchingInd}))) / numTouch(matchingInd);                        
                    end
                end
            end
            
            
            % ANOVA in each configuration            
            anovaPcell = zeros(1,numFeature);
            tunedCell = zeros(1,numFeature);
            tuneAngleCell = nan(1,numFeature);
            tuneModulationCell = zeros(1,numFeature);
            tuneSharpnessCell = zeros(1,numFeature);
            unimodalSingleCell = zeros(1,numFeature);
            unimodalBroadCell = zeros(1,numFeature);
            multimodalCell = zeros(1,numFeature);
            
            anovaVal = cell2mat(spkValAll{1});
            groupAnova = zeros(size(anovaVal));
            angleLengths = [0;cumsum(cellfun(@length, spkValAll{1}))];
            for ai = 1 : length(angles)
                groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
            end
            for i = 1 : size(spkValAll,2)
                anovaVal = cell2mat(spkValAll{i});
                
                [anovaP, ~, anovaStat] = anova1(anovaVal, groupAnova, 'off');
                spkPairComp = multcompare(anovaStat, 'Ctype', anovactype, 'Display', 'off');
                spkMeans = anovaStat.means;
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
                        tuneModulationCell(i) = max(spkMeans) - min(spkMeans);
                        tuneSharpnessCell(i) = spkMeans(tunedAngleInd) - mean(spkMeans(setdiff(1:length(angles), tunedAngleInd)));
                        
                        % Categorization
                        ind__1 = find(spkPairComp(:,1) == tunedAngleInd);
                        ind__2 = find(spkPairComp(:,2) == tunedAngleInd);
                        testInd = union(ind__1, ind__2);
                        insigDiffInd = find(spkPairComp(testInd,6) >= thresholdCategory);
                        sigDiffInd = find(spkPairComp(testInd,6) < thresholdCategory);
                        temp = spkPairComp(testInd(insigDiffInd),1:2);
                        insigDiffIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd, except when there's nothing

                        if isempty(insigDiffIndGroup)
                            unimodalSingleCell(i) = 1;
                        else
                            temp = spkPairComp(testInd(sigDiffInd),1:2);
                            sigDiffIndGroup = setdiff(unique(temp(:)), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.

                            broadInd = intersect(sigInd,insigDiffIndGroup);
                            if length(broadInd) < 2
                                unimodalSingleCell(i) = 1;
                            else                            
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
                                    else
                                        break
                                    end
                                end
                                if broadNum == length(broadInd)
                                    unimodalBroadCell(i) = 1;                                    
                                else
                                    multimodalCell(i) = 1;
                                end
                            end
                        end
                        
                    end
                end
            end
            spkValAllCell(ci,:) = spkValAll;
            anovaPAllCell(ci,:) = anovaPcell;
            tunedAllCell(ci,:) = tunedCell;
            tuneAngleAllCell(ci,:) = tuneAngleCell;
            tuneModulationAllCell(ci,:) = tuneModulationCell;
            tuneSharpnessAllCell(ci,:) = tuneSharpnessCell;
            unimodalSingleAllCell(ci,:) = unimodalSingleCell;
            unimodalBroadAllCell(ci,:) = unimodalBroadCell;
            multimodalAllCell(ci,:) = multimodalCell;
        end
        
        save(savefn, '*AllCell', 'featureNames')
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
