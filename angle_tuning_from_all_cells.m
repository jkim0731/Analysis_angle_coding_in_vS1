% Only within cells responding to touch (from glmFunctionRidgeDE010.mat)

% 1) z-score of calcium response. pole up frames / before pole up frames (calcium)
% 2) # of spikes during touch frames + [0:2] / before pole up (spkPole)
% from all the cells
% run ANOVA
% permutation tests for tuned cells
% final decision of tuned cells
% 
% From tuned cells calculate:
% 1) tuned angle (max abs)
% 2) tuning direction (excited, inhibited, bipolar)
% 3) unimodal-single, unimodal-broad, multimodal, leave-one-out, categorical, ramp
% 4) modultaion: max-min.     sharpness: response of the max - mean response of the rest
% From not-tuned cells calculate:
% 1) response direction (excited, inhibited)
% 2) response amplitude

% 2019/04/08 JK

% Updates:
% include not tuned response cells 2019/05/28 JK

% settings
clear
baseDir = 'D:\TPM\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
% naiveMi = 1:12;
% expertMi = [1,2,3,4,7,9];
% L4Mi = 13:16;

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,9,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];

angles = 45:15:135;
thresholdAnovaP = 0.05; 
thresholdTtestNeighbors = 0.05;
thresholdTtestResponse = 0.05;
thresholdCategory = 0.05;
anovactype = 'hsd';
numResampling = 10000; % permutation test

% Load ridge results file
% It should be at the base directory
cd(baseDir)


for mi = 1 : length(mice)    
% for mi = 1
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
%     for si = 1 : length(sessions{mi})
    for si = 2 : length(sessions{mi})

        session = sessions{mi}(si);
        
        % load uber
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load(ufn)
        
        % still some settings
        savefn = [u.mouseName,u.sessionName,'angle_tuning_all_cell.mat']; %

        % making templates
        % find touch trials 
        touchTrialInd = find(cellfun(@(x) ~isempty(x.protractionTouchChunks), u.trials));
        touchTrialNum = cellfun(@(x) x.trialNum, u.trials(touchTrialInd));
        numPlane = length(u.mimg);
        planeTrialsInd = cell(numPlane,1);
        planeTrialsNum = cell(numPlane,1);
        poleUpFrames = cell(numPlane,1);
        beforePoleUpFrames = cell(numPlane,1);
        touchFrames = cell(numPlane,1);
        nonTouchFrames = cell(numPlane,1); % within pole up frames        
        angleTrialInds = cell(numPlane,1);
        for pi = 1 : numPlane
            planeTrialsInd{pi} = intersect(find(cellfun(@(x) ismember(pi, x.planes), u.trials)), touchTrialInd);
            tempInd = find(u.trials{planeTrialsInd{pi}(1)}.planes == pi);
            planeTrialsNum{pi} = cellfun(@(x) x.trialNum, u.trials(planeTrialsInd{pi}));
            
            poleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} >= x.poleUpTime(1) & x.tpmTime{tempInd} <= x.poleUpTime(end)), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            beforePoleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} < x.poleUpOnsetTime), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);            
            touchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            nonTouchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            for ti = 1 : length(planeTrialsInd{pi})
                tempTrial = u.trials{planeTrialsInd{pi}(ti)};
                tempFrames = cell(1, length(tempTrial.protractionTouchChunks));
                for ptci = 1 : length(tempFrames)
                    tempFrames{ptci} = [0:2] + find(tempTrial.tpmTime{tempInd} >= tempTrial.whiskerTime(tempTrial.protractionTouchChunks{ptci}(1)), 1, 'first');
                end
                tempTouchFrames = unique(cell2mat(tempFrames));
                touchFrames{pi}{ti} = tempTouchFrames(tempTouchFrames <= length(tempTrial.tpmTime{tempInd}));
                nonTouchFrames{pi}{ti} = setdiff(poleUpFrames{pi}{ti}, touchFrames{pi}{ti});
            end
            angleTrialInds{pi} = cell(length(angles),1); % index of planeTrialsInd{pi}
            for ai = 1 : length(angles)
                angleTrialInds{pi}{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(planeTrialsInd{pi})));
            end
        end
        
        touchID = u.cellNums;
        
        ca.touchID = touchID;
        catuned = zeros(length(touchID),1);
        caAnovaPAll = zeros(length(touchID),1);
        catunedAngle = zeros(length(touchID),1);
        catuneDirection = zeros(length(touchID),1);
        caunimodalSingle = zeros(length(touchID),1);
        caunimodalBroad = zeros(length(touchID),1);
        camultimodal = zeros(length(touchID),1);
        caleaveOneOut = zeros(length(touchID),1);
        cacategorical = zeros(length(touchID),1);
        caramp = zeros(length(touchID),1);
        camodulation = zeros(length(touchID),1);
        casharpness = zeros(length(touchID),1);
        caNTR = zeros(length(touchID),1);
        caNTRamplitude = zeros(length(touchID),1); % for not-tuned cells
        caNTRdirection = zeros(length(touchID),1); % for not-tuned cells
        caValAll = cell(length(touchID),1);
        
        spk.touchID = touchID;
        spktuned = zeros(length(touchID),1);
        spkAnovaPAll = zeros(length(touchID),1);
        spktunedAngle = zeros(length(touchID),1);
        spktuneDirection = zeros(length(touchID),1);
        spkunimodalSingle = zeros(length(touchID),1);
        spkunimodalBroad = zeros(length(touchID),1);
        spkmultimodal = zeros(length(touchID),1);
        spkleaveOneOut = zeros(length(touchID),1);
        spkcategorical = zeros(length(touchID),1);
        spkramp = zeros(length(touchID),1);
        spkmodulation = zeros(length(touchID),1);
        spksharpness = zeros(length(touchID),1);
        spkNTR = zeros(length(touchID),1);
        spkNTRamplitude = zeros(length(touchID),1); % for not-tuned cells
        spkNTRdirection = zeros(length(touchID),1); % for not-tuned cells        
        spkValAll = cell(length(touchID),1);
        
        parfor ci = 1:length(touchID)
            fprintf('Processing JK%03d S%02d touch cell %d / %d\n', mouse, session, ci, length(touchID))
            cellNum = touchID(ci);
            plane = floor(cellNum/1000);
            trialInds = planeTrialsInd{plane};
            calciumPoleUpFrames = poleUpFrames{plane};
            spkTouchFrames = touchFrames{plane};
            baselineFrames = beforePoleUpFrames{plane};
            angleInds = angleTrialInds{plane}; % index of trialInds
            
            % z-score calcium signal in each cell
            cind = find(u.trials{trialInds(1)}.neuindSession == cellNum);
            tempF = cellfun(@(x) x.dF(cind,:), u.trials(trialInds), 'uniformoutput', false);
            meanF = mean(cell2mat(tempF'));
            stdF = std(cell2mat(tempF'));
            caZscore = cellfun(@(x) (x-meanF)/stdF, tempF, 'uniformoutput', false);
            
            % all spikes
            tempSpk = cellfun(@(x) x.spk(cind,:), u.trials(trialInds), 'uniformoutput', false);
            
            caValAll{ci} = cell(length(angles),1);
            spkValAll{ci} = cell(length(angles),1);
            for ai = 1 : length(angles)
                trialAngleInd = angleInds{ai};               
                caValAll{ci}{ai} = zeros(length(trialAngleInd),1);
                spkValAll{ci}{ai} = zeros(length(trialAngleInd),1);
                for ti = 1 : length(trialAngleInd)
                    tempInd = trialAngleInd(ti);
                    caValAll{ci}{ai}(ti) = mean(tempF{tempInd}(calciumPoleUpFrames{tempInd})) - mean(tempF{tempInd}(baselineFrames{tempInd}));
                    spkValAll{ci}{ai}(ti) = mean(tempSpk{tempInd}(spkTouchFrames{tempInd})) - mean(tempSpk{tempInd}(baselineFrames{tempInd}));
                end
            end
            
            caVal = caValAll{ci};
            spkVal = spkValAll{ci};
            %% ANOVA
            caAnovaVal = cell2mat(caVal);
            spkAnovaVal = cell2mat(spkVal);
            if length(caAnovaVal) ~= length(spkAnovaVal)
                error('calcium and spike has different lengths')
            end
            if ~isempty(find(isnan(caAnovaVal))) || ~isempty(find(isnan(spkAnovaVal)))
                error('nan values')
            end
            groupAnova = zeros(size(caAnovaVal));
            angleLengths = [0;cumsum(cellfun(@length, caVal))];
            for ai = 1 : length(angles)
                groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
            end
                        
            [caAnovaP, ~, caAnovaStat] = anova1(caAnovaVal, groupAnova, 'off');
            caAnovaPAll(ci) = caAnovaP;
            [spkAnovaP, ~, spkAnovaStat] = anova1(spkAnovaVal, groupAnova, 'off');
            spkAnovaPAll(ci) = spkAnovaP;
            caPairComp = multcompare(caAnovaStat, 'Ctype', anovactype, 'Display', 'off');
            spkPairComp = multcompare(spkAnovaStat, 'Ctype', anovactype, 'Display', 'off');
            caMeans = caAnovaStat.means;
            spkMeans = spkAnovaStat.means;
            
            %% Start with calcium first
            tempH = cellfun(@(x) ttest(x), caVal);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            tunedFlag = 0;
            if ~isnan(caAnovaP)
                if caAnovaP < thresholdAnovaP && ~isempty(sigInd)
                    % permutation test
                    maxmod = max(caMeans) - min(caMeans);
                    permAnovaP = zeros(numResampling,1);
                    permmaxmod = zeros(numResampling,1);
    %                 parfor ri = 1 : numResampling
                    for ri = 1 : numResampling
                        tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                        [permAnovaP(ri), ~, permStats] = anova1(caAnovaVal, tempG, 'off');
                        permmaxmod(ri) = max(permStats.means) - min(permStats.means);
                    end

                    if length(find(permAnovaP < caAnovaP)) < 0.05 * numResampling % passed permutation test. Tuned.
                        tunedFlag = 1; % necessary to have this for parfor loop
                        catuned(ci) = 1;
                        [~, maxind] = max(abs(caMeans(sigInd)));
                        tunedAngleInd = sigInd(maxind);
                        catunedAngle(ci) = angles(tunedAngleInd);

                        maxVal = max(caMeans(sigInd));
                        minVal = min(caMeans(sigInd));
                        if minVal > 0
                            catuneDirection(ci) = 1;
                        elseif maxVal < 0 
                            catuneDirection(ci) = 2;
                        elseif maxVal > 0 && minVal < 0
                            catuneDirection(ci) = 3;
                        else
                            catuneDirection(ci) = -1; % error
                        end
                        camodulation(ci) = max(caMeans) - min(caMeans);
                        casharpness(ci) = caMeans(tunedAngleInd) - mean(caMeans(setdiff(1:length(angles), tunedAngleInd)));

                        % Categorization
                        ind__1 = find(caPairComp(:,1) == tunedAngleInd);
                        ind__2 = find(caPairComp(:,2) == tunedAngleInd);
                        testInd = union(ind__1, ind__2);
                        insigDiffInd = find(caPairComp(testInd,6) >= thresholdCategory);
                        sigDiffInd = find(caPairComp(testInd,6) < thresholdCategory);
                        temp = caPairComp(testInd(insigDiffInd),1:2);
                        insigDiffIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd, except when there's nothing

                        if isempty(insigDiffIndGroup)
                            caunimodalSingle(ci) = 1;
                        else
                            broadInd = intersect(sigInd,insigDiffIndGroup);
                            if length(broadInd) < 2
                                caunimodalSingle(ci) = 1;
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
                                    end
                                end
                                if broadNum == length(broadInd)
                                    caunimodalBroad(ci) = 1;
                                    % if broad, then it can be a categorical
                                    center = (length(angles)+1) / 2;
                                    compInd = union(find(caPairComp(:,1) == tunedAngleInd), find(caPairComp(:,2) == tunedAngleInd));
                                    indMat = caPairComp(compInd,1:2);
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
                                    if isempty(find(caPairComp(compInd(withinInd),6) < thresholdCategory, 1)) && ... % nothing within the same half is different from the max ind
                                            isempty(find(caPairComp(compInd(betweenInd),6) >= thresholdCategory, 1)) % nothing between different half is same with the max ind
                                        cacategorical(ci) = 1; % categorical (>= 90 or <= 90)
                                    end
                                else
                                    camultimodal(ci) = 1;
                                end
                            end
                            temp = caPairComp(testInd(sigDiffInd),1:2);
                            sigIndGroup = setdiff(temp(:), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.
                            if ~isempty(find(diff(insigDiffIndGroup)>1,1))
                                if sum(tempH(sigIndGroup))
                                    camultimodal(ci) = 1; % multimodal. Including bipolar.
                                end
                                if length(sigIndGroup) == 1 && ... % only one bin is significantly different from the tuned bin. (can't be larger in response because of the way tuned bin is defined)
                                        all(tempH(insigDiffIndGroup)) % and all insignicant indices are different from 0
                                    caleaveOneOut(ci) = 1 ; % leave-one-out. Part of multimodal in definition.
                                end
                            end
                            if isempty(find(diff(sign(diff(caMeans))),1)) % everything is going up or down 
                                caramp(ci) = 1; % ramping up or down
                            end
                        end
                    end
                end
            end
            if tunedFlag == 0 % not tuned.
                tval = ttest(caAnovaVal)
                % NTR: not tuned response
                if ~isnan(tval)
                    if tval % if the response is different from 0
                        caNTR(ci) = 1;
                        if mean(caAnovaVal) > 0
                            caNTRdirection(ci) = 1;
                        else
                            caNTRdirection(ci) = 2;
                        end
                        caNTRamplitude(ci) = mean(caAnovaVal);
                    end
                end
            end
            
            %% Then with spikes
            tempH = cellfun(@(x) ttest(x), spkVal);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            tunedFlag = 0;
            if ~isnan(spkAnovaP)
                if spkAnovaP < thresholdAnovaP && ~isempty(sigInd)
                    % permutation test
                    maxmod = max(spkMeans) - min(spkMeans);
                    permAnovaP = zeros(numResampling,1);
                    permmaxmod = zeros(numResampling,1);
    %                 parfor ri = 1 : numResampling
                    for ri = 1 : numResampling
                        tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                        [permAnovaP(ri), ~, permStats] = anova1(spkAnovaVal, tempG, 'off');
                        permmaxmod(ri) = max(permStats.means) - min(permStats.means);
                    end

                    if length(find(permAnovaP < spkAnovaP)) < 0.05 * numResampling % passed permutation test. Tuned.
                        tunedFlag = 1; % necessary to have this for parfor loop
                        spktuned(ci) = 1;
                        [~, maxind] = max(abs(spkMeans(sigInd)));
                        tunedAngleInd = sigInd(maxind);
                        spktunedAngle(ci) = angles(tunedAngleInd);

                        maxVal = max(spkMeans(sigInd));
                        minVal = min(spkMeans(sigInd));
                        if minVal > 0
                            spktuneDirection(ci) = 1;
                        elseif maxVal < 0 
                            spktuneDirection(ci) = 2;
                        elseif maxVal > 0 && minVal < 0
                            spktuneDirection(ci) = 3;
                        else
                            spktuneDirection(ci) = -1; % error
                        end
                        spkmodulation(ci) = max(spkMeans) - min(spkMeans);
                        spksharpness(ci) = spkMeans(tunedAngleInd) - mean(spkMeans(setdiff(1:length(angles), tunedAngleInd)));

                        % Categorization
                        ind__1 = find(spkPairComp(:,1) == tunedAngleInd);
                        ind__2 = find(spkPairComp(:,2) == tunedAngleInd);
                        testInd = union(ind__1, ind__2);
                        insigDiffInd = find(spkPairComp(testInd,6) >= thresholdCategory);
                        sigDiffInd = find(spkPairComp(testInd,6) < thresholdCategory);
                        temp = spkPairComp(testInd(insigDiffInd),1:2);
                        insigDiffIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd, except when there's nothing

                        if isempty(insigDiffIndGroup)
                            spkunimodalSingle(ci) = 1;
                        else
                            broadInd = intersect(sigInd,insigDiffIndGroup);
                            if length(broadInd) < 2
                                spkunimodalSingle(ci) = 1;
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
                                    spkunimodalBroad(ci) = 1;
                                    % if broad, then it can be a categorical
                                    center = (length(angles)+1) / 2;
                                    compInd = union(find(spkPairComp(:,1) == tunedAngleInd), find(spkPairComp(:,2) == tunedAngleInd));
                                    indMat = spkPairComp(compInd,1:2);
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
                                    if isempty(find(spkPairComp(compInd(withinInd),6) < thresholdCategory, 1)) && ... % nothing within the same half is different from the max ind
                                            isempty(find(spkPairComp(compInd(betweenInd),6) >= thresholdCategory, 1)) % nothing between different half is same with the max ind
                                        spkcategorical(ci) = 1; % categorical (>= 90 or <= 90)
                                    end
                                else
                                    spkmultimodal(ci) = 1;
                                end
                            end
                            temp = spkPairComp(testInd(sigDiffInd),1:2);
                            sigIndGroup = setdiff(temp(:), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.
                            if ~isempty(find(diff(insigDiffIndGroup)>1,1))
                                if sum(tempH(insigDiffIndGroup))>1 % to exclude tuned angle
                                    spkmultimodal(ci) = 1; % multimodal. Including bipolar.
                                end
                                if length(sigIndGroup) == 1 && ... % only one bin is significantly different from the tuned bin. (can't be larger in response because of the way tuned bin is defined)
                                        all(tempH(insigDiffIndGroup)) % and all insignicant indices are different from 0
                                    spkleaveOneOut(ci) = 1 ; % leave-one-out. Part of multimodal in definition.
                                end
                            end
                            if isempty(find(diff(sign(diff(spkMeans))),1)) % everything is going up or down 
                                spkramp(ci) = 1; % ramping up or down
                            end
                        end
                    end
                end
            end
            
            if tunedFlag == 0 % not tuned
                % NTR: not tuned response
                tval = ttest(spkAnovaVal);
                if ~isnan(tval)
                    if tval % if the response is different from 0
                        spkNTR(ci) = 1;
                        if mean(spkAnovaVal) > 0
                            spkNTRdirection(ci) = 1;
                        else
                            spkNTRdirection(ci) = 2;
                        end
                        spkNTRamplitude(ci) = mean(spkAnovaVal);
                    end
                end
            end
        end
        
        ca.tuned = catuned;
        ca.tunedAngle = catunedAngle;
        ca.anovaP = caAnovaPAll;
        ca.tuneDirection = catuneDirection;
        ca.unimodalSingle = caunimodalSingle;
        ca.unimodalBroad = caunimodalBroad;
        ca.multimodal = camultimodal;
        ca.leaveOneOut = caleaveOneOut;
        ca.categorical = cacategorical;
        ca.ramp = caramp;
        ca.modulation = camodulation;
        ca.sharpness = casharpness;
        ca.NTR = caNTR;
        ca.NTRamplitude = caNTRamplitude;
        ca.NTRdirection = caNTRdirection;
        ca.val = caValAll;
        
        spk.tuned = spktuned;
        spk.anovaP = spkAnovaPAll;
        spk.tunedAngle = spktunedAngle;
        spk.tuneDirection = spktuneDirection;
        spk.unimodalSingle = spkunimodalSingle;
        spk.unimodalBroad = spkunimodalBroad;
        spk.multimodal = spkmultimodal;
        spk.leaveOneOut = spkleaveOneOut;
        spk.categorical = spkcategorical;
        spk.ramp = spkramp;
        spk.modulation = spkmodulation;
        spk.sharpness = spksharpness;
        spk.NTR = spkNTR;
        spk.NTRamplitude = spkNTRamplitude;
        spk.NTRdirection = spkNTRdirection;
        spk.val = spkValAll;
        
        info.cellID = u.cellNums;
        info.celly = u.celly;
        info.cellx = u.cellx;
        info.c2ypoints = u.c2ypoints;
        info.c2xpoints = u.c2xpoints;
        info.fovsize = u.fovsize;
        info.fovxrange = u.fovxrange;
        info.fovyrange = u.fovyrange;
        info.fovdepth = u.fovdepth;

        save(savefn, 'ca','spk','info')

    end
end