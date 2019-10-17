% Modified from angle_tuning_preAnswer_perTouch_spkOnly
% Consider pre-first lick, instead of pre-answer
% It is because of high correlation between lick and angle in expert mice
% 2019/10/16 JK

% Modified from angle_tuning_predecision
% Only within cells responding to touch (from glmFunctionLasso_NC.mat)

% ( # of spikes during touch frames + [0:1] - before pole up (spkPole) ) / # of touches
% Only from trials with touch
% run ANOVA
% permutation tests from tuned cells
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
% 2019/06/25 - Only pre-decision periods. addition of frames were just 1,
% instead of 2. Just compensating for a possible spike detection frame
% error (at most single is assumed). Consider touch chunks "by whisking".

% 2019/09/27 JK
% Copied and modified from angle_tuning_predecision
% Major change: 
% 1) Only consider inferred spikes
% 2) Divide (Sum diff Spikes #) by (# of touches). Average per-touch response in each trial.
% 3) naming: predecision -> preAnswer. we don't know when mice decided.
% 4) Removed one of multimodal selection. Before, it was sorted as multimodal even when the insigDiff angle response is insignificant from 0. 
% Possibly, this led to some (though very minor) overlap between multimodal and unimodalBroad. Before, in this case, it was sorted as unimodalBroad, not multimodal.
% 5) Confirmed that tpmTime calculation (frame assignment) was correct.

%%
% settings
clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];

angles = 45:15:135;
thresholdAnovaP = 0.05; 
thresholdPermAnovaP = 0.05;
thresholdTtestNeighbors = 0.05;
thresholdTtestResponse = 0.05;
thresholdCategory = 0.05;
anovactype = 'hsd';
numResampling = 10000; % permutation test

% Load lasso results file
% It should be at the base directory
cd(baseDir)
load('cellFunctionLasso_NC.mat')

for mi = 4 : length(mice)    
% for mi = 2
    mouse = mice(mi);
    cd(sprintf('%s%03d',baseDir,mouse))
    for si = 1 : length(sessions{mi})
%     for si = 2
        session = sessions{mi}(si);
        
        % load uber
        ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
        load(ufn)
        
        % still some settings
        savefn = [u.mouseName,u.sessionName,'angle_tuning_lasso_preLick_perTouch_spkOnly_NC.mat']; %

        % making templates
        % find preLick touch trials
        firstLickTime = cell(length(u.trials),1);
        allLicks = cellfun(@(x) union(union(union(x.leftLickTime, x.rightLickTime), x.answerLickTime), x.poleDownOnsetTime), u.trials, 'un', 0);
        lickIndsAfterPoleIn = cellfun(@(x,y) find(x>y.poleUpOnsetTime,1), allLicks, u.trials);
        % poleUpTime(1) should be better, but there is almost no
        % difference in the result, and it's hard to deal with catch trials
        % 2019/10/16 JK
        for licki = 1 : length(firstLickTime)
            firstLickTime{licki} = allLicks{licki}(lickIndsAfterPoleIn(licki));
        end
        tempTouchTrialInd = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        plTouchInd = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y, u.trials(tempTouchTrialInd), firstLickTime(tempTouchTrialInd)));
        touchTrialInd = tempTouchTrialInd(plTouchInd);
        numPlane = length(u.mimg);
        planeTrialsInd = cell(numPlane,1);
        planeTrialsNum = cell(numPlane,1);
        poleUpFrames = cell(numPlane,1);
        beforePoleUpFrames = cell(numPlane,1);
        touchFrames = cell(numPlane,1);        
        nonTouchFrames = cell(numPlane,1); % within pole up frames. Need for confirmation whether the cell is really touch-responsive, compared to general task-responsive.    
        numTouchPreFirstLick = cell(numPlane,1);
        angleTrialInds = cell(numPlane,1);
        for pi = 1 : numPlane
            planeTrialsInd{pi} = intersect(find(cellfun(@(x) ismember(pi, x.planes), u.trials)), touchTrialInd);
            tempInd = find(u.trials{planeTrialsInd{pi}(1)}.planes == pi);
            planeTrialsNum{pi} = cellfun(@(x) x.trialNum, u.trials(planeTrialsInd{pi}));
            
            poleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} >= x.poleUpTime(1) & x.tpmTime{tempInd} <= x.poleUpTime(end)), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);
            beforePoleUpFrames{pi} = cellfun(@(x) find(x.tpmTime{tempInd} < x.poleUpOnsetTime), u.trials(planeTrialsInd{pi}), 'uniformoutput', false);            
            touchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            nonTouchFrames{pi} = cell(length(planeTrialsInd{pi}),1);
            numTouchPreFirstLick{pi} = zeros(length(planeTrialsInd{pi}),1);
            for ti = 1 : length(planeTrialsInd{pi})
                tempTrial = u.trials{planeTrialsInd{pi}(ti)};
                tempFirstLickTime = firstLickTime{planeTrialsInd{pi}(ti)};
                preFirstLickInd = find(cellfun(@(x) tempTrial.whiskerTime(x(1)) < tempFirstLickTime, tempTrial.protractionTouchChunksByWhisking));
                
                tempFrames = cell(1, length(preFirstLickInd));
                
                for ptci = 1 : length(tempFrames)
                    tempFrames{ptci} = [0:1] + find(tempTrial.tpmTime{tempInd} >= tempTrial.whiskerTime(tempTrial.protractionTouchChunksByWhisking{ptci}(1)), 1, 'first');
                    % tpmTime is the beginning timepoint of each frame from the trial start point (TTL1 signal). 
                    % Considering Ca2+ and GCaMP signal rise time, any event should be assigned to frames having the time point >= to that event timing.
                    % (ref, C:\Users\shires\Documents\GitHub\jksbx\+Calcium\@CalciumDataArray\CalciumDataArray.m)
                    % +1 frame for delayed response (~ > 120 ms)
                end
                touchFrames{pi}{ti} = unique(cell2mat(tempFrames));
                nonTouchFrames{pi}{ti} = setdiff(poleUpFrames{pi}{ti}:find(tempTrial.tpmTime{tempInd} < tempFirstLickTime), touchFrames{pi}{ti});
                numTouchPreFirstLick{pi}(ti) = length(preFirstLickInd);
            end
            angleTrialInds{pi} = cell(length(angles),1); % index of planeTrialsInd{pi}
            for ai = 1 : length(angles)
                angleTrialInds{pi}{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(planeTrialsInd{pi})));
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
        spkNTamplitude = zeros(length(touchID),1); % for not-tuned cells
        spkNTdirection = zeros(length(touchID),1); % for not-tuned cells
        spkValAll = cell(length(touchID),1);
        spkValAllNontouch = cell(length(touchID),1);
        spkNumTouch = cell(length(touchID),1);
        spkNumTouchFrames = cell(length(touchID),1);
        spkNumNontouchFrames = cell(length(touchID),1);
        
        parfor ci = 1:length(touchID)
%         for ci = 195
            fprintf('Processing JK%03d S%02d touch cell %d / %d\n', mouse, session, ci, length(touchID))
            cellNum = touchID(ci);
%             if mouse == 27 && (session == 9 || session == 10)
%                 plane = floor(cellNum/1000) - 4;
%             else
                plane = floor(cellNum/1000);
%             end
            trialInds = planeTrialsInd{plane};
            spkTouchFrames = touchFrames{plane};
            spkNontouchFrames = nonTouchFrames{plane};
            baselineFrames = beforePoleUpFrames{plane};
            angleInds = angleTrialInds{plane}; % index of trialInds
            numTouch = numTouchPreFirstLick{plane};
            
            % all spikes
            cind = find(u.trials{trialInds(1)}.neuindSession == cellNum);
            tempSpk = cellfun(@(x) x.spk(cind,:), u.trials(trialInds), 'uniformoutput', false);

            spkValAll{ci} = cell(length(angles),1);
            spkValAllNontouch{ci} = cell(length(angles),1);
            spkNumTouch{ci} = cell(length(angles),1);
            spkNumTouchFrames{ci} = cell(length(angles),1);
            spkNumNontouchFrames{ci} = cell(length(angles),1);
            
            for ai = 1 : length(angles)
                trialAngleInd = angleInds{ai};                               
                spkValAll{ci}{ai} = zeros(length(trialAngleInd),1);
                spkValAllNontouch{ci}{ai} = nan(length(trialAngleInd),1);
                spkNumTouch{ci}{ai} = zeros(length(trialAngleInd),1);
                spkNumTouchFrames{ci}{ai} = zeros(length(trialAngleInd),1);
                spkNumNontouchFrames{ci}{ai} = zeros(length(trialAngleInd),1);
                for ti = 1 : length(trialAngleInd)
                    tempInd = trialAngleInd(ti);
                    spkValAll{ci}{ai}(ti) = sum( tempSpk{tempInd}(spkTouchFrames{tempInd}) - mean(tempSpk{tempInd}(baselineFrames{tempInd})) ) / numTouch(tempInd);
                    % Delta inferred spike per touch 2019/09/27
                    spkNumTouch{ci}{ai}(ti) = numTouch(tempInd);
                    spkNumTouchFrames{ci}{ai}(ti) = length(spkTouchFrames{tempInd});
                    spkNumNontouchFrames{ci}{ai}(ti) = length(spkNontouchFrames{tempInd});
                    if ~isempty(spkNontouchFrames{tempInd})
                        spkValAllNontouch{ci}{ai}(ti) = mean(tempSpk{tempInd}(spkNontouchFrames{tempInd})) - mean(tempSpk{tempInd}(baselineFrames{tempInd})) ;
                    end
                end
            end

            spkVal = spkValAll{ci};
            
            %% ANOVA
            spkAnovaVal = cell2mat(spkVal);
            if ~isempty(find(isnan(spkAnovaVal)))
                error('nan values')
            end
            groupAnova = zeros(size(spkAnovaVal));
            angleLengths = [0;cumsum(cellfun(@length, spkVal))];
            for ai = 1 : length(angles)
                groupAnova(angleLengths(ai)+1:angleLengths(ai+1)) = deal(ai);
            end
                        
            [spkAnovaP, ~, spkAnovaStat] = anova1(spkAnovaVal, groupAnova, 'off');
            spkAnovaPAll(ci) = spkAnovaP;
            spkPairComp = multcompare(spkAnovaStat, 'Ctype', anovactype, 'Display', 'off');
            spkMeans = spkAnovaStat.means;

            %% Then with spikes
            tempH = cellfun(@(x) ttest(x), spkVal);
            tempH(isnan(tempH)) = deal(0);
            sigInd = find(tempH); % significant indices
            if spkAnovaP < thresholdAnovaP && ~isempty(sigInd)
                % permutation test
                maxmod = max(spkMeans) - min(spkMeans);
                permAnovaP = zeros(numResampling,1);
%                 permmaxmod = zeros(numResampling,1);
%                 parfor ri = 1 : numResampling
                for ri = 1 : numResampling
                    tempG = groupAnova(randperm(length(groupAnova),length(groupAnova)));
                    permAnovaP(ri) = anova1(spkAnovaVal, tempG, 'off');
%                     [permAnovaP(ri), ~, permStats] = anova1(spkAnovaVal, tempG, 'off');
%                     permmaxmod(ri) = max(permStats.means) - min(permStats.means);
                end

                if length(find(permAnovaP < spkAnovaP)) >= thresholdPermAnovaP * numResampling % failed to pass permutation test
                    % NT: not tuned
                    if mean(spkAnovaVal) > 0
                        spkNTdirection(ci) = 1;
                    else
                        spkNTdirection(ci) = 2;
                    end
                    spkNTamplitude(ci) = mean(spkAnovaVal);
                else % passed permutation test. Tuned.
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
                        temp = spkPairComp(testInd(sigDiffInd),1:2);
                        sigDiffIndGroup = setdiff(unique(temp(:)), tunedAngleInd); % exclude tunedAngleInd. Any index that is significantly different from the tuned angle index.
                        
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
%                                 categorical is possible only in unimodalBroad
                                if isempty(find(spkPairComp(compInd(withinInd),6) < thresholdCategory, 1)) && ... % nothing within the same half is different from the max ind
                                        isempty(find(spkPairComp(compInd(betweenInd),6) >= thresholdCategory, 1)) % nothing between different half is same with the max ind
                                    spkcategorical(ci) = 1; % categorical (>= 90 or <= 90)
                                end
                            else
                                spkmultimodal(ci) = 1;
%                                 leaveOneOut is possible only in multimodal
                                if length(sigDiffIndGroup) == 1 && ... % only one bin is significantly different from the tuned bin. (can't be larger in response because of the way tuned bin is defined)
                                    all(tempH(insigDiffIndGroup)) % and all insignicant indices are different from 0
                                    spkleaveOneOut(ci) = 1 ; % leave-one-out. Part of multimodal in definition.
                                end
                            end
                        end
                        
                        if isempty(find(diff(sign(diff(spkMeans))),1)) % everything is going up or down 
                            spkramp(ci) = 1; % ramping up or down
                        end
                    end
                end
            else % NT: not tuned
                if mean(spkAnovaVal) > 0
                    spkNTdirection(ci) = 1;
                else
                    spkNTdirection(ci) = 2;
                end
                spkNTamplitude(ci) = mean(spkAnovaVal);
            end            
        end
        
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
        spk.NTamplitude = spkNTamplitude;
        spk.NTdirection = spkNTdirection;
        spk.val = spkValAll;
        spk.valNontouchFrames = spkValAllNontouch;
        spk.numTouch = spkNumTouch;
        spk.numTouchFrames = spkNumTouchFrames;
        spk.numNontouchFrames = spkNumNontouchFrames;
        
        info.cellID = u.cellNums;
        info.celly = u.celly;
        info.cellx = u.cellx;
        info.c2ypoints = u.c2ypoints;
        info.c2xpoints = u.c2xpoints;
        info.fovsize = u.fovsize;
        info.fovxrange = u.fovxrange;
        info.fovyrange = u.fovyrange;
        info.fovdepth = u.fovdepth;

        save(savefn, 'spk','info')
    end
end