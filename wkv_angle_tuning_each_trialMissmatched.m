%%
% Temporary fix for JK038 S02, where touch glm was run and then one trial was removed because of w3 fix (about dome-like mask and kappaV projection plane).
% WKV glm was run after that fix again, but not touch GLM (just because it takes too long, and also changing it will change very minute results too.)
% There is one less trial in WKV model (#55, and uber too), and treat each model differently to see if that makes it work.


%%
% Making spike matrix as in u, but from a model with removing selected whisker kinematic variables.
% This is to show which variables are important in contructing angle tuning.
% Only consider previously angle-tuned cells.
% In different settings, run angle tuning, and record if the cell is deemed tuned or not.

% First, build spike matrix from the full model.
% Then, build a full model from touch glm.
% After that, build a full model from wkv glm.


%%

baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% mice = [38];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
% sessions = {[2]};
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

% naiveModelTune = struct;
% expertModelTune = struct;
% for mi = 1 : length(mice)
%% this has to match with the structure of 'cellFunctionLasso_NC'
for mi = 6 
%%
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    %% this has to match with the structure of 'cellFunctionLasso_NC'
    for si = 1
    %%
        session = sessions{mi}(si);        
        savefn = sprintf('angle_tuning_model_lasso_NC_JK%03dS%02d', mouse, session);
        
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
        
        % load touch glm results and get coefficients (for all active cells)
        % remove the corresponding predictors (trial # 55, index 54)
        glmfn = sprintf('glmResponseType_JK%03dS%02d_lasso_NC_R10',mouse, session); % these are only from all active cells
        load(glmfn, 'cIDAll', 'fitCoeffs', 'allPredictors')
        ap2 = allPredictors;
        
%% Here comes the treatment for JK038 S02 error trial (trial # 55 of lower plane is included in ap2, but not in ap1 and u)
        lowerplaneInd = find(cellfun(@(x) ismember(5,x.planes), u.trials));
        lowerplaneTrialNums = cellfun(@(x) x.trialNum, u.trials(lowerplaneInd));
        ind2removeFromAp2 = find(lowerplaneTrialNums > 55, 1, 'first');
        % now, remove this index from ap2{5:8};
        % first, divide ap2{5:8} into trials using nans. calculate from 5th
        % plane and then apply for planes 5-8 
        naninds = find(isnan(ap2{5}(:,1)));
        nanindChanges = find(diff(naninds)>1);
        numTrials = length(nanindChanges);
        ap2temp = ap2;
        % confirmed that this numTrials is 1 larger than
        % lowerplaneTrialNums
        for i = 5 : 8
            ap2{i}(naninds(nanindChanges(ind2removeFromAp2))-4+1 : naninds(nanindChanges(ind2removeFromAp2+1))-4,:) = [];
        end
%% test if I removed the right one
%     nanindsAfter = find(isnan(ap2{5}(:,1)));
%     nanindChangesAfter = find(diff(nanindsAfter)>1);
% 
%     compare(ap2temp{5}(naninds(nanindChanges(ind2removeFromAp2-1))-4+1:naninds(nanindChanges(ind2removeFromAp2))-4,:), ...
%         ap2{5}(nanindsAfter(nanindChangesAfter(ind2removeFromAp2-1))-4+1:nanindsAfter(nanindChangesAfter(ind2removeFromAp2))-4,:))
% 
%     compare(ap2temp{5}(naninds(nanindChanges(ind2removeFromAp2+1))-4+1:naninds(nanindChanges(ind2removeFromAp2+2))-4,:), ...
%         ap2{5}(nanindsAfter(nanindChangesAfter(ind2removeFromAp2))-4+1:nanindsAfter(nanindChangesAfter(ind2removeFromAp2+1))-4,:))
%     
% naninds1 = find(isnan(ap1{5}(:,1)));
% naninds2 = find(isnan(ap2{5}(:,1)));
% compare(naninds1, naninds2)
%% The rest continues the same         
        cIDAllTouch = cIDAll;
        coeffs = zeros(length(fitCoeffs),length(fitCoeffs{1}),10);
        for gi = 1 : 10
            glmfn = sprintf('glmResponseType_JK%03dS%02d_lasso_NC_R%02d',mouse, session, gi); % these are only from touch response cells
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
            if ~isempty(ap2{i}) % it can be empty in case of JK027 S09 and S10
                naninds = find(isnan(ap2{i}(:,1)));
                nanindChanges = find(diff(naninds)>1);
                numTrials = length(nanindChanges);
                allPredictorsTouch{i} = cell(numTrials,1);
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
        
%         total 15 sets for testing angle tuning
%         (1) raw dF
%         (2) full model with touch angle
%         (3) full model with wkv
%         (4-15) from wkv, remove each WKV        
%         (16-27) from wkv model without any wkv, add each WKV
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
        % Treat WKV and Touch GLM the same, since the corresponding error
        % trial in touch GLM model is gone.
        allPredictorsTouchAngleInds = cell(8,1); % index of ap1{i} & ap2{i}
        for i = 1 : 8
            tempIndTouch = find(cellfun(@(x) ismember(i, x.planes), u.trials));
            allPredictorsTouchAngleInds{i} = cell(length(angles),1);
            for ai = 1 : length(angles)
                allPredictorsTouchAngleInds{i}{ai} = find(ismember( tempIndTouch, intersect(intersect(tempIndTouch, touchTrialInd), ...
                    find(cellfun(@(x) x.angle == angles(ai), u.trials))) ));                
            end
        end
        
        % settings for variables to be saved later
%         spkValAllCell = cell(length(touchID),17);
        spkValAllCell = cell(length(touchID),27);
        anovaPAllCell = zeros(length(touchID),27);
        tunedAllCell = zeros(length(touchID),27);
        tuneAngleAllCell = zeros(length(touchID),27);
        
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
            modelTouchAngleInds = allPredictorsTouchAngleInds{plane}; % index of allPredictorsTouch & WKV
            
            spkValAll = cell(1,27);
            
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
            
            % whiskerTouchMat = [maxDkappaHMat, maxDkappaVMat, maxDthetaMat, maxDphiMat, maxSlideDistanceMat, maxDurationMat, ...    
%                             thetaAtTouchMat, phiAtTouchMat, kappaHAtTouchMat, kappaVAtTouchMat, arcLengthAtTouchMat, touchCountMat];

            % from full wkv glm (#3 ~ #15)
            coeffInd = find(cIDAllWKV == cellNum);
            inds = cell(25,1);
            coeffLength = size(meanCoeffsWKV,2);
            inds{1} = 1:coeffLength;
            inds{2} = setdiff(1:coeffLength, 8:10); % maxDthetaMat
            inds{3} = setdiff(1:coeffLength, 11:13); % maxDphiMat
            inds{4} = setdiff(1:coeffLength, 2:4); % maxDkappaHMat
            inds{5} = setdiff(1:coeffLength, 5:7); % maxDkappaVMat
            inds{6} = setdiff(1:coeffLength, 14:16); % maxSlideDistanceMat
            inds{7} = setdiff(1:coeffLength, 17:19); % maxDurationMat
            inds{8} = setdiff(1:coeffLength, 20:22); % thetaAtTouchMat
            inds{9} = setdiff(1:coeffLength, 23:25); % phiAtTouchMat
            inds{10} = setdiff(1:coeffLength, 26:28); % kappaHAtTouchMat
            inds{11} = setdiff(1:coeffLength, 29:31); % kappaVAtTouchMat
            inds{12} = setdiff(1:coeffLength, 32:34); % arcLengthAtTouchMat
            inds{13} = setdiff(1:coeffLength, 35:37); % touchCountMat
            
            inds{14} = [1,8:10];
            inds{15} = [1,11:13];
            inds{16} = [1,2:4];
            inds{17} = [1,5:7];
            inds{18} = [1,14:16];
            inds{19} = [1,17:19];
            inds{20} = [1,20:22];
            inds{21} = [1,23:25];
            inds{22} = [1,26:28];
            inds{23} = [1,29:31];
            inds{24} = [1,32:34];
            inds{25} = [1,35:37];

            for i = 3 : 27
                spkValAll{i} = cell(length(angles),1);
            end
            
            anovaPcell = zeros(1,27);
            tunedCell = zeros(1,27);
            tuneAngleCell = nan(1,27);
            for ai = 1 : length(angles)
                trialAngleInd = modelTouchAngleInds{ai};
                for i = 3 : 27
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
                    
                    for i = 1 : 25
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
