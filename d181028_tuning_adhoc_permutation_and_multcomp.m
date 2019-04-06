clear
% mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};  
mice = [38,41];
sessions = {[2],[3]};

baseDir = 'Y:\Whiskernas\JK\suite2p\';

%% for calcium
percentFailPerm = zeros(sum(cellfun(@(x) length(x), sessions)),1);
touchResponse = zeros(sum(cellfun(@(x) length(x), sessions)),1);
failNum = zeros(sum(cellfun(@(x) length(x), sessions)),1);
percentFailSharpness = zeros(sum(cellfun(@(x) length(x), sessions)),1);

for mi = 1 : length(mice)
% for mi = 6:length(mice)  
    cd([baseDir, sprintf('%03d',mice(mi))])
    for si = 1 : length(sessions{mi})
        mouse = mice(mi);
        session = sessions{mi}(si);
        fprintf('Processing JK%03d S%02d \n', mice(mi), sessions{mi}(si));
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load(ufn)
%         u = Uber.buildUberArray(mice(mi), sessions{mi}(si));        
        frameRate = u.frameRate;
        fn = [u.mouseName,u.sessionName,'singleCell_anova_calcium.mat']; %
        savefn = [u.mouseName,u.sessionName,'singleCell_anova_calcium_final.mat']; 
        permfn = [u.mouseName, u.sessionName, 'permtest_calcium.mat'];
        load(fn);
        load(permfn)
        
        if length(fail) == length(cellsTuned) % meaning that permutation results are not applied yet
            currentInd = sum(cellfun(@(x) length(x), sessions(1:(mi-1))))+1;
            failNum(currentInd) = sum(fail);
            percentFailPerm(currentInd) = sum(fail)/length(fail);

            permfailInd = find(fail);
            cellsFailPerm = cellsTuned(permfailInd);
            changedToNTR = [];
            changedToNTNR = [];
            for ci = 1 : sum(fail)
                cell = cellsTuned(permfailInd(ci));
                cellind = find(cellsTotal == cell);
                dF = cell2mat(dFtotal{cellind});
                response = nanmean(dF(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2);
                if ttest(response)
                    tempThreshold = [max(thresholdResponsePositive, u.noise(cellid)), min(thresholdResponseNegative, -u.noise(cellid)/2)];
                    if mean(response) > tempThreshold(1) || mean(response) < tempThreshold(2)
                        changedToNTR = [changedToNTR; cell];
                        
                        cellsNTResponse = [cellsNTResponse;cell];
                        tempTimeseries = nanmean(dF(:,baseFrameNum+1:baseFrameNum+afterFrameNum),1);
                        NTRAmplitude = [NTRAmplitude; abs(mean(response))];
                        if mean(response) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif mean(response) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries > min(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; 0]; % error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(dF,1)
                            tempTS = dF(ri,baseFrameNum + 1 : baseFrameNum + afterFrameNum);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(response) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(response) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(dF(:,baseFrameNum:end));
                        rho = zeros(size(dF,1),1);
                        for ri = 1 : size(dF,1)
                            tempTS = dF(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cell];
                        changedToNTNR = [changedToNTNR; cell];
                    end
                else
                    cellsNTNR = [cellsNTNR; cell];
                    changedToNTNR = [changedToNTNR; cell];
                end                
            end
            cellsTuned(permfailInd) = [];
            tuneAmplitude(permfailInd) = [];
            tuneAngle(permfailInd) = [];
            tuneDirection(permfailInd) = [];
            tuneMaxResponseTimepoint(permfailInd) = [];
            tuneModulationMaxmin(permfailInd) = [];
            tuneReliability(permfailInd) = [];
            tuneResponseProb(permfailInd) = [];
            tuneSharpness(permfailInd) = [];
            tuneCateg = setdiff(tuneCateg, cellsFailPerm);
            tuneLOO = setdiff(tuneLOO, cellsFailPerm);
            tuneMM = setdiff(tuneMM, cellsFailPerm);
            tuneRampStrict = setdiff(tuneRampStrict, cellsFailPerm);
            touchResponse(currentInd) = length(changedToNTR);
            cellsNTNR = sort(cellsNTNR);
            [cellsNTResponse, sortind] = sort(cellsNTResponse);
            NTRAmplitude = NTRAmplitude(sortind);
            NTRdirection = NTRdirection(sortind);
            NTRMaxResponseTimepoint = NTRMaxResponseTimepoint(sortind);
            NTRresponseProb = NTRresponseProb(sortind);
            NTRreliability = NTRreliability(sortind);
        end
        
        sharpfailInd = find(tuneSharpness == length(angles));
        if ~isempty(sharpfailInd)
            cellsFailSharp = cellsTuned(sharpfailInd);
            percentFailPerm(currentInd) = length(sharpfailInd)/length(cellsTuned);
            for ci = 1 : length(sharpfailInd)
                cell = cellsTuned(sharpfailInd(ci));
                cellind = find(cellsTotal == cell);
                dF = cell2mat(dFtotal{cellind});
                response = nanmean(dF(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2);
                if ttest(response)
                    tempThreshold = [max(thresholdResponsePositive, u.noise(cellid)), min(thresholdResponseNegative, -u.noise(cellid)/2)];
                    if mean(response) > tempThreshold(1) || mean(response) < tempThreshold(2)
                        changedToNTR = [changedToNTR; cell];
                        
                        cellsNTResponse = [cellsNTResponse;cell];
                        tempTimeseries = nanmean(dF(:,baseFrameNum+1:baseFrameNum+afterFrameNum),1);
                        NTRAmplitude = [NTRAmplitude; abs(mean(response))];
                        if mean(response) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif mean(response) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries > min(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; 0]; % error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(dF,1)
                            tempTS = dF(ri,baseFrameNum + 1 : baseFrameNum + afterFrameNum);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(response) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(response) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(dF(:,baseFrameNum:end));
                        rho = zeros(size(dF,1),1);
                        for ri = 1 : size(dF,1)
                            tempTS = dF(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cell];
                        changedToNTNR = [changedToNTNR; cell];
                    end
                else
                    cellsNTNR = [cellsNTNR; cell];
                    changedToNTNR = [changedToNTNR; cell];
                end 
            end
            cellsTuned(sharpfailInd) = [];
            tuneAmplitude(sharpfailInd) = [];
            tuneAngle(sharpfailInd) = [];
            tuneDirection(sharpfailInd) = [];
            tuneMaxResponseTimepoint(sharpfailInd) = [];
            tuneModulationMaxmin(sharpfailInd) = [];
            tuneReliability(sharpfailInd) = [];
            tuneResponseProb(sharpfailInd) = [];
            tuneSharpness(sharpfailInd) = [];
            tuneCateg = setdiff(tuneCateg, cellsFailSharp);
            tuneLOO = setdiff(tuneLOO, cellsFailSharp);
            tuneMM = setdiff(tuneMM, cellsFailSharp);
            tuneRampStrict = setdiff(tuneRampStrict, cellsFailSharp);
            touchResponse(currentInd) = length(changedToNTR);
            cellsNTNR = sort(cellsNTNR);
            [cellsNTResponse, sortind] = sort(cellsNTResponse);
            NTRAmplitude = NTRAmplitude(sortind);
            NTRdirection = NTRdirection(sortind);
            NTRMaxResponseTimepoint = NTRMaxResponseTimepoint(sortind);
            NTRresponseProb = NTRresponseProb(sortind);
            NTRreliability = NTRreliability(sortind);
        end
            
        
        save(savefn, 'cell*','tune*','NTR*', 'dFtotal', '*ctype', 'angles', 'baseFrameNum', 'afterFrameNum', 'threshold*', ...
            'excludeDrinkingTime', 'onlyBeforeDecision', 'onlyAfterDecision', 'allowOverlap', 'onlyFirstTouch', 'ksh*', 'ksp*', 'leveneTest', 'frameRate', 'oneSampleH', '*P', 'noise', 'c2*points', 'fov*',...
            'changed*', '*fail*', 'anovaP', 'permmaxmod', 'touchResponse', 'percentFail*')
    end
end

%% for spikes
percentFailPerm = zeros(sum(cellfun(@(x) length(x), sessions)),1);
touchResponse = zeros(sum(cellfun(@(x) length(x), sessions)),1);
failNum = zeros(sum(cellfun(@(x) length(x), sessions)),1);
percentFailSharpness = zeros(sum(cellfun(@(x) length(x), sessions)),1);

for mi = 1 : length(mice)
% for mi = 6:length(mice)  
    cd([baseDir, sprintf('%03d',mice(mi))])
    for si = 1 : length(sessions{mi})        
        mouse = mice(mi);
        session = sessions{mi}(si);
        fprintf('Processing JK%03d S%02d \n', mice(mi), sessions{mi}(si));
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load(ufn)
%         u = Uber.buildUberArray(mice(mi), sessions{mi}(si));        
        frameRate = u.frameRate;
        fn = [u.mouseName,u.sessionName,'singleCell_anova_spk.mat']; %
        savefn = [u.mouseName,u.sessionName,'singleCell_anova_spk_final.mat']; 
        permfn = [u.mouseName, u.sessionName, 'permtest_spk.mat'];
        load(fn);
        load(permfn)
        
        if length(fail) == length(cellsTuned) % meaning that permutation results are not applied yet
            currentInd = sum(cellfun(@(x) length(x), sessions(1:(mi-1))))+1;
            failNum(currentInd) = sum(fail);
            percentFailPerm(currentInd) = sum(fail)/length(fail);

            permfailInd = find(fail);
            cellsFailPerm = cellsTuned(permfailInd);
            changedToNTR = [];
            changedToNTNR = [];
            for ci = 1 : sum(fail)
                cell = cellsTuned(permfailInd(ci));
                cellind = find(cellsTotal == cell);
                infspk = cell2mat(spkTotal{cellind});
                response = nanmean(infspk(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2);
                if ttest(response)
                    tempThreshold = [max(thresholdResponsePositive, u.noise(cellid)), min(thresholdResponseNegative, -u.noise(cellid)/2)];
                    if mean(response) > tempThreshold(1) || mean(response) < tempThreshold(2)
                        changedToNTR = [changedToNTR; cell];
                        
                        cellsNTResponse = [cellsNTResponse;cell];
                        tempTimeseries = nanmean(infspk(:,baseFrameNum+1:baseFrameNum+afterFrameNum),1);
                        NTRAmplitude = [NTRAmplitude; abs(mean(response))];
                        if mean(response) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif mean(response) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries > min(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; 0]; % error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(infspk,1)
                            tempTS = infspk(ri,baseFrameNum + 1 : baseFrameNum + afterFrameNum);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(response) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(response) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(infspk(:,baseFrameNum:end));
                        rho = zeros(size(infspk,1),1);
                        for ri = 1 : size(infspk,1)
                            tempTS = infspk(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cell];
                        changedToNTNR = [changedToNTNR; cell];
                    end
                else
                    cellsNTNR = [cellsNTNR; cell];
                    changedToNTNR = [changedToNTNR; cell];
                end                
            end
            cellsTuned(permfailInd) = [];
            tuneAmplitude(permfailInd) = [];
            tuneAngle(permfailInd) = [];
            tuneDirection(permfailInd) = [];
            tuneMaxResponseTimepoint(permfailInd) = [];
            tuneModulationMaxmin(permfailInd) = [];
            tuneReliability(permfailInd) = [];
            tuneResponseProb(permfailInd) = [];
            tuneSharpness(permfailInd) = [];
            tuneCateg = setdiff(tuneCateg, cellsFailPerm);
            tuneLOO = setdiff(tuneLOO, cellsFailPerm);
            tuneMM = setdiff(tuneMM, cellsFailPerm);
            tuneRampStrict = setdiff(tuneRampStrict, cellsFailPerm);
            touchResponse(currentInd) = length(changedToNTR);
            cellsNTNR = sort(cellsNTNR);
            [cellsNTResponse, sortind] = sort(cellsNTResponse);
            NTRAmplitude = NTRAmplitude(sortind);
            NTRdirection = NTRdirection(sortind);
            NTRMaxResponseTimepoint = NTRMaxResponseTimepoint(sortind);
            NTRresponseProb = NTRresponseProb(sortind);
            NTRreliability = NTRreliability(sortind);
        end
        
        sharpfailInd = find(tuneSharpness == length(angles));
        if ~isempty(sharpfailInd)
            cellsFailSharp = cellsTuned(sharpfailInd);
            percentFailPerm(currentInd) = length(sharpfailInd)/length(cellsTuned);
            for ci = 1 : length(sharpfailInd)
                cell = cellsTuned(sharpfailInd(ci));
                cellind = find(cellsTotal == cell);
                infspk = cell2mat(spkTotal{cellind});
                response = nanmean(infspk(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2);
                if ttest(response)
                    tempThreshold = [max(thresholdResponsePositive, u.noise(cellid)), min(thresholdResponseNegative, -u.noise(cellid)/2)];
                    if mean(response) > tempThreshold(1) || mean(response) < tempThreshold(2)
                        changedToNTR = [changedToNTR; cell];
                        
                        cellsNTResponse = [cellsNTResponse;cell];
                        tempTimeseries = nanmean(infspk(:,baseFrameNum+1:baseFrameNum+afterFrameNum),1);
                        NTRAmplitude = [NTRAmplitude; abs(mean(response))];
                        if mean(response) > tempThreshold(1)
                            NTRdirection = [NTRdirection; 1]; % increase
                            mt = find(tempTimeseries > max(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        elseif mean(response) < tempThreshold(2)
                            NTRdirection = [NTRdirection; 2]; % decrease
                            mt = find(tempTimeseries > min(tempTimeseries)*0.9,1);
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; mt / frameRate];
                        else
                            NTRdirection = [NTRdirection; 0]; % error
                            NTRMaxResponseTimepoint = [NTRMaxResponseTimepoint; NaN];
                        end

                        % Response probability
                        responseNum = 0;
                        for ri = 1 : size(infspk,1)
                            tempTS = infspk(ri,baseFrameNum + 1 : baseFrameNum + afterFrameNum);
                            tempTS = tempTS(~isnan(tempTS));
                            [~, p] = ttest(tempTS');
                            if p < thresholdTtestResponse
                                if (mean(response) > tempThreshold(1) && mean(tempTS) > tempThreshold(1)) || ...
                                        (mean(response) < tempThreshold(2) && mean(tempTS) < tempThreshold(2))
                                    responseNum = responseNum + 1;
                                end
                            end
                        end
                        NTRresponseProb = [NTRresponseProb; responseNum/ri*100];

                        % Response reliability
                        template = nanmean(infspk(:,baseFrameNum:end));
                        rho = zeros(size(infspk,1),1);
                        for ri = 1 : size(infspk,1)
                            tempTS = infspk(ri,baseFrameNum:end);
                            inds = find(~isnan(tempTS));
                            rho(ri) = corr(template(inds)', tempTS(inds)');
                        end
                        NTRreliability = [NTRreliability; mean(rho)];
                    else
                        cellsNTNR = [cellsNTNR; cell];
                        changedToNTNR = [changedToNTNR; cell];
                    end
                else
                    cellsNTNR = [cellsNTNR; cell];
                    changedToNTNR = [changedToNTNR; cell];
                end 
            end
            cellsTuned(sharpfailInd) = [];
            tuneAmplitude(sharpfailInd) = [];
            tuneAngle(sharpfailInd) = [];
            tuneDirection(sharpfailInd) = [];
            tuneMaxResponseTimepoint(sharpfailInd) = [];
            tuneModulationMaxmin(sharpfailInd) = [];
            tuneReliability(sharpfailInd) = [];
            tuneResponseProb(sharpfailInd) = [];
            tuneSharpness(sharpfailInd) = [];
            tuneCateg = setdiff(tuneCateg, cellsFailSharp);
            tuneLOO = setdiff(tuneLOO, cellsFailSharp);
            tuneMM = setdiff(tuneMM, cellsFailSharp);
            tuneRampStrict = setdiff(tuneRampStrict, cellsFailSharp);
            touchResponse(currentInd) = length(changedToNTR);
            cellsNTNR = sort(cellsNTNR);
            [cellsNTResponse, sortind] = sort(cellsNTResponse);
            NTRAmplitude = NTRAmplitude(sortind);
            NTRdirection = NTRdirection(sortind);
            NTRMaxResponseTimepoint = NTRMaxResponseTimepoint(sortind);
            NTRresponseProb = NTRresponseProb(sortind);
            NTRreliability = NTRreliability(sortind);
        end
            
        
        save(savefn, 'cell*','tune*','NTR*', 'spkTotal', '*ctype', 'angles', 'baseFrameNum', 'afterFrameNum', 'threshold*', ...
            'excludeDrinkingTime', 'onlyBeforeDecision', 'onlyAfterDecision', 'allowOverlap', 'onlyFirstTouch', 'ksh*', 'ksp*', 'leveneTest', 'frameRate', 'oneSampleH', '*P', 'noise', 'c2*points', 'fov*',...
            'changed*', '*fail*', 'anovaP', 'permmaxmod', 'touchResponse', 'percentFail*')
    end
end