function glm = glm_results_cell_function_touch(mouse, session, baseDir)

%% For figure 3. GLM result plots fro function assignment in each cell
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
%     10 repeated lasso results
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
% 2019/04/02 JK


% Updates:
% Forget about tuning here. It can't be well inferred from glm. 2019/04/08 JK

%% basic settings
% chi2pvalThreshold = 0.001; % less than 0.001 for fitting
deThreshold = 0.1; % include 0.1 as fit
coeffThreshold = 0; % include 0.01 as a coefficient
repeat = 10;
glm = struct;
% ca = struct;
% spk = struct;
L4depth = 350; % include 350 um as L4 (350 is the starting point)
% angles = 45:15:135;

%% dependent settings
ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_lasso_NC_R', mouse, session);
% cafn = sprintf('JK%03dS%02dsingleCell_anova_calcium_final', mouse, session);
% spkfn = sprintf('JK%03dS%02dsingleCell_anova_spk_final', mouse, session);

%% load uber
cd(sprintf('%s%03d',baseDir, mouse))
load(ufn, 'u') % loading u
u = u;
%% select cells with average DE > deThreshold (0.1) and average coefficients
averageDE = zeros(length(u.cellNums),1);
allCoeff = cell(length(u.cellNums),repeat);
for ri = 1 : repeat
    load(sprintf('%s%02d',glmfnBase, ri), 'fitCoeffs', 'fitDevExplained')
    averageDE = averageDE + fitDevExplained/repeat;
    allCoeff(:,ri) = fitCoeffs;
end
load(sprintf('%s%02d',glmfnBase, repeat), 'allPredictors', 'indPartial', 'posShift')
allPredictors = allPredictors;
indPartial = indPartial;
posShift = posShift;
averageCoeff = cell(length(u.cellNums),1);
for ci = 1 : length(u.cellNums)
    tempCoeff = mean(cell2mat(allCoeff(ci,:)),2);
    tempCoeff(abs(tempCoeff) < coeffThreshold) = deal(0);
    averageCoeff{ci} = tempCoeff;
end
deFitInd = find(averageDE >= deThreshold);

%% assigning functions to each cell
cellFunction = cell(length(deFitInd), 1);
parfor ci = 1 : length(deFitInd)
% for ci = 620
    cID = u.cellNums(deFitInd(ci));
    tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
    planeInd = floor(cID/1000);
    
    testInput = allPredictors{planeInd};
    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));
    
    if length(testInput) ~= length(spkTest)
        error(sprintf('input matrix and spike length mismatch in ci %d', ci))
    end

    finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
    spkTest = spkTest(finiteIndTest);

    coeff = averageCoeff{deFitInd(ci)};
    coeffInds = find(coeff(2:end));
    model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*coeff);
    mu = mean(spkTest); % null poisson parameter
    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));                    
    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
    devExplained = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    dfFull = length(coeffInds);
    
    tempFit = zeros(1,length(indPartial) + 1);
%     if devExplained >= deThreshold && devianceFullNull > chi2inv(1-chi2pvalThreshold, dfFull) % re-evaluation led to full model fit
    if devExplained >= deThreshold  % for ridge
        tempFit(1) = 1;
        for pi = 1 : length(indPartial)
            if sum(ismember(coeffInds, indPartial{pi})) > 0
                resInd = setdiff(1:size(testInput,2), indPartial{pi});
                temptest = testInput(finiteIndTest,resInd);
                partialModel = exp([ones(length(finiteIndTest),1),temptest]*coeff([1,resInd+1]));
                partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
                devExpPartial = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                DEdiff = devExplained - devExpPartial;
                devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
                dfPartial = length(find(coeff(resInd+1)));
%                 if DEdiff >= deThreshold && devianceFullPartial > chi2inv(1-chi2pvalThreshold, dfFull - dfPartial)
                if DEdiff >= deThreshold % for ridge    
                    tempFit(pi+1) = 1;
                end
            end
        end
        if sum(tempFit(2:end)) == 0 % when there is no function fit, the cell is regarded not fit
            tempFit(1) = 0;
        end
    end
    cellFunction{ci} = tempFit;
end

%% re-assigning based on the calculation above, and assigning 
deFit = cell2mat(cellFunction);
finalFitInd = find(deFit(:,1)); % indices from average DE>0.1 cells
fitInd = deFitInd(finalFitInd); % indices from all cells
glm.cellFitID = u.cellNums(fitInd);
glm.cellFunction = cellfun(@(x) x(2:end), cellFunction(finalFitInd), 'uniformoutput', false);
% averageCoeff = averageCoeff(fitInd);

%% transfer information from u to glm
glm.cellFitIndC2 = find(u.isC2(fitInd));
glm.cellFitIndL23 = find(u.cellDepths(fitInd) < L4depth);
glm.cellFitDepths = u.cellDepths(fitInd);
glm.cellFitxpoint = u.cellx(fitInd);
glm.cellFitypoint = u.celly(fitInd);
glm.cellNums = u.cellNums(fitInd);
glm.cellDepths = u.cellDepths(fitInd);
glm.isC2 = u.isC2(fitInd);
glm.allCell.cellNums = u.cellNums;
glm.allCell.cellDepths = u.cellDepths;
glm.allCell.isC2 = u.isC2;

%% function assignment
glm.touchID = glm.cellFitID(find(cellfun(@(x) x(1), glm.cellFunction)));
glm.soundID = glm.cellFitID(find(cellfun(@(x) x(2), glm.cellFunction)));
glm.rewardID = glm.cellFitID(find(cellfun(@(x) x(3), glm.cellFunction)));
glm.whiskingID = glm.cellFitID(find(cellfun(@(x) x(4), glm.cellFunction)));
glm.lickingID = glm.cellFitID(find(cellfun(@(x) x(5), glm.cellFunction)));

