function glm = glm_results_cell_function_shuffling_WTV_ONLY(mouse, session, baseDir)

%% For figure 3. GLM result plots fro function assignment in each cell
%%
% same as in "glm_results_cell_function_shuffling", except that
% touch angles are replaced by whisker touch variables

% Notes:
%     Assume u is fixed. No change of u.cellNums (very important in indexing)
%     After running 10 repeats of glmnet_cell_tupe.m
%     Based on d190322_glm_repetition.m (and former checking codes)
%
% Input: 
%     mouse
%     session
%     baseDir
% 
% Files:
%     10 repeated ridge results
%     Uber file
%     ca anova results
%     spk anova results
% 
% Output:
%     Proportion of touch cells, whisking cells, touch & whisking cells
%     Proportion of angle-tuned touch cells.
%     Comparison with anova test results
%     These in C2, non-C2, L2/3, L4, L2/3 C2, L2/3 non-C2, (To compare with Peron et al), L4 C2, and L4 non-C2
% 
%     Proportion of sound cells, reward cells, and licking cells
% 
%     Example map of cells in different function assignment (information for these)
% 
%     glm.cellFitID
%     glm.cellFitIndC2
%     glm.cellFitIndL23
%     glm.cellFitDepths
%     glm.cellFitxpoint
%     glm.cellFitypoint
%     glm.cellNums = u.cellNums
%     glm.cellDepths = u.cellDepths
%     glm.isC2 = u.isC2
%     glm.cellFunction : 1 - touch, 2 - sound, 3 - reward, 4 - whisking, 5 - licking
%     glm.tunedID : cell ID from u.cellNum, that has tuning (either excited or inhibited, whichever is larger in absolute value)
%     glm.tunedAngle : best tuned angle (45-135). Same length as tuned. Calculated by area under curve of each angle (same as summation of the coefficients)
%     glm.tuneDirection : 1 - excited, 2 - inhibited
%     glm.touchID : cell ID of touch response cells, EXCLUDING tuned cells
%     ca.tunedID
%     ca.tunedAngle
%     ca.tuneDireaction
%     ca.touchID
%     spk.tunedID
%     spk.tunedAngle
%     spk.tuneDirection
%     spk.touchID
% 

% 2019/04/02 JK


% Updates:
% Forget about tuning here. It can't be well inferred from glm. 2019/04/08 JK
%
% Instead of excluding one by one, shuffle the predictors. 
% This is because when excluding one group of coefficients, the intercept does not explain the partial model well
% (veryfied by the fact that inclusion method gives negative deviance explained)
% When shuffling, maintain trial sequences. Shuffle only within trials.
% 
% 2109/04/11 JK

%% basic settings
% chi2pvalThreshold = 0.001; % less than 0.001 for fitting
deThreshold = 0.05; % include 0.1 as fit
coeffThreshold = 0; % include 0.01 as a coefficient
repeat = 10;
glm = struct;
% ca = struct;
% spk = struct;
L4depth = 350; % include 350 um as L4 (350 is the starting point)
% angles = 45:15:135;

%% dependent settings
ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
glmfnBase = sprintf('glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R', mouse, session);
% cafn = sprintf('JK%03dS%02dsingleCell_anova_calcium_final', mouse, session);
% spkfn = sprintf('JK%03dS%02dsingleCell_anova_spk_final', mouse, session);

%% load uber
cd(sprintf('%s%03d',baseDir, mouse))
load(ufn, 'u') % loading u
u = u;
%% select cells with average DE > deThreshold (0.1) and average coefficients
load(sprintf('%s%02d',glmfnBase, repeat), 'cIDAll', 'allPredictors', 'indPartial', 'posShift')
cIDAll = cIDAll;
allPredictors = allPredictors;
indPartial = indPartial;
posShift = posShift;
numCells = length(cIDAll);

averageDE = zeros(numCells,1);
allCoeff = cell(numCells,repeat);
for ri = 1 : repeat
    load(sprintf('%s%02d',glmfnBase, ri), 'fitCoeffs', 'fitDevExplained')
    averageDE = averageDE + fitDevExplained/repeat;
    allCoeff(:,ri) = fitCoeffs;
end
averageCoeff = cell(numCells,1);
for ci = 1 : numCells
    averageCoeff{ci} = mean(cell2mat(allCoeff(ci,:)),2);
end
% deFitInd = find(averageDE >= deThreshold);

%% assigning functions to each cell
deviance = zeros(numCells,1);
errorRatio = cell(numCells,1);
devExp = zeros(numCells,1);
DEdiff = cell(numCells,1);
exclusionER = cell(numCells,1);

whiskerVariableER = cell(numCells,1);
whiskerVariableExclusionER = zeros(numCells,13);
whiskerVariableDEdiff = zeros(numCells,13);

parfor ci = 1 : numCells
    if ~isempty(averageCoeff{ci})
        fprintf('Processing %d/%d\n', ci, numCells)
        cID = cIDAll(ci);
        tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
        cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
        planeInd = floor(cID/1000);

        testInput = allPredictors{planeInd};
        spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

        if length(testInput) ~= length(spkTest)
            error('input matrix and spike length mismatch')
        end

        finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
        spkTest = spkTest(finiteIndTest);

        coeff = averageCoeff{ci};
        model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*coeff);
        mu = mean(spkTest); % null poisson parameter
        nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
        saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));                    
        fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
        tempDeviance = 2 * (fullLogLikelihood - nullLogLikelihood);
        devExplained = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);

    %     tempFit = zeros(1,length(indPartial) + 1);
    %     if devExplained >= deThreshold  % for ridge
    %         tempFit(1) = 1;

            numPermute = 100;
            tempPartialDEsub = zeros(1,length(indPartial));
            tempExclusionER = zeros(1,length(indPartial));

            permER = zeros(length(indPartial),numPermute);
            permWTV = zeros(13,numPermute);
            tempWVDEdiff = zeros(1,13);
            tempWTVexclusionER = zeros(1,13);

            permER = zeros(length(indPartial),numPermute);
            permWTV = zeros(6,numPermute);
            tempWVDEdiff = zeros(1,11);
            tempWTVexclusionER = zeros(1,11);

            for pi = 1 : length(indPartial)
                %% exclusion method
                partialInds = setdiff(1:length(coeff), indPartial{pi}+1); % including intercept
                partialCoeffs = coeff(partialInds);
                partialModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,partialInds(2:end)-1)]*partialCoeffs);
                partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
                partialDevExp = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                tempPartialDEsub(pi) = devExplained - partialDevExp;            
                tempExclusionER(pi) = (saturatedLogLikelihood - partialLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);

                %% permutation method
                switch pi
                    case 1 % whisker touch variables. shifts 3. 13 groups in sequential shift                    
                        aptest = testInput(:,indPartial{pi}(1));

                        nonanind = find(~isnan(aptest));
                        numGroups = length(find(diff(nonanind)>1));
                        indIntervals = [0;find(diff(nonanind)>1)]; % (i)+1:(i+1)

                        indGroups = cell(numGroups,1);

                        randGroups = cell(numGroups,numPermute);
                        for gi = 1 : numGroups
                            indGroups{gi} = nonanind(indIntervals(gi)+1:indIntervals(gi+1));    
                            for ri = 1 : numPermute
                                randGroups{gi,ri} = indGroups{gi}(randperm(length(indGroups{gi})));
                            end
                        end
                        randGroups = cell2mat(randGroups);
                        indGroups = cell2mat(indGroups);

                        for ri = 1 : numPermute
                            tempPartialInputNodelay = testInput(:,indPartial{pi}(1:3:37));
                            tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                            tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), size(tempPartialInputNodelay,2)*3);
                            for di = 1 : 3
                                tempPartialInputAll(:, di : 3 : di+36 ) = ...
                                    circshift(tempPartialInputNodelay, [0 di-1]);
                            end
                            tempInput = testInput;
                            tempInput(:,indPartial{pi}) = tempPartialInputAll;
                            permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                            permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                            permER(pi,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end

                        % comparing between whisker touch variables. There are
                        % 6 of them currently
                        % 11 of them 2019/04/16 JK
                        % 13 of them 2019/04/18 JK

                        for ri = 1 : numPermute
                            for j = 1 : 13
                                tempPartialInputNodelay = testInput(:,indPartial{pi}((j-1)*3+1));
                                tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                                tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), 3);
                                for di = 1 : 3
                                    tempPartialInputAll(:, di ) = ...
                                        circshift(tempPartialInputNodelay, [0 di-1]);
                                end
                                tempInput = testInput;
                                tempInput(:,indPartial{pi}((j-1)*3+1:j*3)) = tempPartialInputAll;
                                permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                                permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                                permWTV(j,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                            end
                        end
                        
                        for j = 1 : 13
                            partialInds = setdiff(1:length(coeff), indPartial{pi}((j-1)*3+1:j*3) + 1); % including intercept
                            partialCoeffs = coeff(partialInds);
                            partialModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,partialInds(2:end)-1)]*partialCoeffs);
                            partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
                            partialDevExp = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                            tempWVDEdiff(j) = devExplained - partialDevExp;            
                            tempWTVexclusionER(j) = (saturatedLogLikelihood - partialLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end
                        
                    case 2 % in other cases, it repeats in every segment
                        % sound. shifts 3
                        aptest = testInput(:,indPartial{pi}(1));

                        nonanind = find(~isnan(aptest));
                        numGroups = length(find(diff(nonanind)>1));
                        indIntervals = [0;find(diff(nonanind)>1)]; % (i)+1:(i+1)

                        indGroups = cell(numGroups,1);

                        randGroups = cell(numGroups,numPermute);
                        for gi = 1 : numGroups
                            indGroups{gi} = nonanind(indIntervals(gi)+1:indIntervals(gi+1));    
                            for ri = 1 : numPermute
                                randGroups{gi,ri} = indGroups{gi}(randperm(length(indGroups{gi})));
                            end
                        end
                        randGroups = cell2mat(randGroups);
                        indGroups = cell2mat(indGroups);

                        for ri = 1 : numPermute
                            tempPartialInputNodelay = testInput(:,indPartial{pi}(1));                    
                            tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                            tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), size(tempPartialInputNodelay,2)*4);
                            for di = 1 : 4
                                tempPartialInputAll(:,(di-1)*size(tempPartialInputNodelay,2)+1 : di*size(tempPartialInputNodelay,2)) = ...
                                    circshift(tempPartialInputNodelay, [0 di-1]);
                            end
                            tempInput = testInput;
                            tempInput(:,indPartial{pi}) = tempPartialInputAll;
                            permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                            permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                            permER(pi,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end
                    case 3
                        % reward. shifts 3. grouped in angles, so same as in touch
                        aptest = testInput(:,indPartial{pi}(1));

                        nonanind = find(~isnan(aptest));
                        numGroups = length(find(diff(nonanind)>1));
                        indIntervals = [0;find(diff(nonanind)>1)]; % (i)+1:(i+1)

                        indGroups = cell(numGroups,1);

                        randGroups = cell(numGroups,numPermute);
                        for gi = 1 : numGroups
                            indGroups{gi} = nonanind(indIntervals(gi)+1:indIntervals(gi+1));    
                            for ri = 1 : numPermute
                                randGroups{gi,ri} = indGroups{gi}(randperm(length(indGroups{gi})));
                            end
                        end
                        randGroups = cell2mat(randGroups);
                        indGroups = cell2mat(indGroups);

                        for ri = 1 : numPermute
                            tempPartialInputNodelay = testInput(:,indPartial{pi}(1:8));                    
                            tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                            tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), size(tempPartialInputNodelay,2)*4);
                            for di = 1 : 4
                                tempPartialInputAll(:,(di-1)*size(tempPartialInputNodelay,2)+1 : di*size(tempPartialInputNodelay,2)) = ...
                                    circshift(tempPartialInputNodelay, [0 di-1]);
                            end
                            tempInput = testInput;
                            tempInput(:,indPartial{pi}) = tempPartialInputAll;
                            permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                            permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                            permER(pi,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end
                    case 4
                        % whisking. shifts 7. 3 groups in sequential shift
                        aptest = testInput(:,indPartial{pi}(1));

                        nonanind = find(~isnan(aptest));
                        numGroups = length(find(diff(nonanind)>1));
                        indIntervals = [0;find(diff(nonanind)>1)]; % (i)+1:(i+1)

                        indGroups = cell(numGroups,1);

                        randGroups = cell(numGroups,numPermute);
                        for gi = 1 : numGroups
                            indGroups{gi} = nonanind(indIntervals(gi)+1:indIntervals(gi+1));    
                            for ri = 1 : numPermute
                                randGroups{gi,ri} = indGroups{gi}(randperm(length(indGroups{gi})));
                            end
                        end
                        randGroups = cell2mat(randGroups);
                        indGroups = cell2mat(indGroups);

                        for ri = 1 : numPermute
                            tempPartialInputNodelay = testInput(:,indPartial{pi}(1:7:21));
                            tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                            tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), size(tempPartialInputNodelay,2)*7);
                            for di = 1 : 7
                                tempPartialInputAll(:, di : 7 : di+14 ) = ...
                                    circshift(tempPartialInputNodelay, [0 di-1]);
                            end
                            tempInput = testInput;
                            tempInput(:,indPartial{pi}) = tempPartialInputAll;
                            permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                            permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                            permER(pi,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end
                    case 5
                        % licking. shifts 4. 2 groups in sequential shift
                        aptest = testInput(:,indPartial{pi}(1));

                        nonanind = find(~isnan(aptest));
                        numGroups = length(find(diff(nonanind)>1));
                        indIntervals = [0;find(diff(nonanind)>1)]; % (i)+1:(i+1)

                        indGroups = cell(numGroups,1);

                        randGroups = cell(numGroups,numPermute);
                        for gi = 1 : numGroups
                            indGroups{gi} = nonanind(indIntervals(gi)+1:indIntervals(gi+1));    
                            for ri = 1 : numPermute
                                randGroups{gi,ri} = indGroups{gi}(randperm(length(indGroups{gi})));
                            end
                        end
                        randGroups = cell2mat(randGroups);
                        indGroups = cell2mat(indGroups);

                        for ri = 1 : numPermute
                            tempPartialInputNodelay = testInput(:,indPartial{pi}(1:4:5));
                            tempPartialInputNodelay(indGroups,:) = tempPartialInputNodelay(randGroups(:,ri),:);
                            tempPartialInputAll = zeros(size(tempPartialInputNodelay,1), size(tempPartialInputNodelay,2)*4);
                            for di = 1 : 4
                                tempPartialInputAll(:, di : 4 : di+4 ) = ...
                                    circshift(tempPartialInputNodelay, [0 di-1]);
                            end
                            tempInput = testInput;
                            tempInput(:,indPartial{pi}) = tempPartialInputAll;
                            permModel = exp([ones(length(finiteIndTest),1),tempInput(finiteIndTest,:)]*coeff);
                            permLogLikelihood = sum(log(poisspdf(spkTest',permModel)));
                            permER(pi,ri) = (saturatedLogLikelihood - permLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                        end
                end

            end
    %         if sum(tempFit(2:end)) == 0 % when there is no function fit, the cell is regarded not fit
    %             tempFit(1) = 0;
    %         end
    %     end
    %     cellFunction{ci} = tempFit;
        whiskerVariableDEdiff(ci,:) = tempWVDEdiff;
        whiskerVariableExclusionER(ci,:) = tempWTVexclusionER;
        deviance(ci) = tempDeviance;
        errorRatio{ci} = permER;
        devExp(ci) = devExplained;
        DEdiff{ci} = tempPartialDEsub;
        exclusionER{ci} = tempExclusionER;
        whiskerVariableER{ci} = permWTV;
    end
end

glm.cID = cIDAll;
glm.deviance = deviance;
glm.errorRatio = errorRatio;
glm.devExp = devExp;
glm.DEdiff = DEdiff;
glm.exclusionER = exclusionER;
glm.whiskerVariableER = whiskerVariableER;
glm.whiskerVariableExclusionER = whiskerVariableExclusionER;
glm.whiskerVariableDEdiff = whiskerVariableDEdiff;


