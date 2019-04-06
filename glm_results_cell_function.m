function [glm, ca, spk] = glm_results_cell_function(mouse, session, baseDir)

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

%% basic settings
chi2pvalThreshold = 0.001; % less than 0.001 for fitting
deThreshold = 0.05; % include 0.1 as fit
coeffThreshold = 0.01; % include 0.01 as a coefficient
repeat = 10;
glm = struct;
ca = struct;
spk = struct;
L4depth = 350; % include 350 um as L4 (350 is the starting point)
angles = 45:15:135;

%% dependent settings
ufn = sprintf('UberJK%03dS%02d',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_m44_R', mouse, session);
cafn = sprintf('JK%03dS%02dsingleCell_anova_calcium_final', mouse, session);
spkfn = sprintf('JK%03dS%02dsingleCell_anova_spk_final', mouse, session);

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
cellFitInd = find(averageDE >= deThreshold);

%% assigning functions to each cell
cellFunction = cell(length(cellFitInd), 1);
parfor ci = 1 : length(cellFitInd)
    cID = u.cellNums(cellFitInd(ci));
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

    coeff = averageCoeff{cellFitInd(ci)};
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
    if devExplained >= deThreshold && devianceFullNull > chi2inv(1-chi2pvalThreshold, dfFull) % re-evaluation led to full model fit
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
                if DEdiff >= deThreshold && devianceFullPartial > chi2inv(1-chi2pvalThreshold, dfFull - dfPartial)
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
allFit = cell2mat(cellFunction);
allFitInd = find(allFit(:,1));
fitInd = cellFitInd(allFitInd);
glm.cellFitID = u.cellNums(fitInd);
glm.cellFunction = cellfun(@(x) x(2:end), cellFunction(allFitInd), 'uniformoutput', false);
averageCoeff = averageCoeff(allFitInd);

%% transfer information from u to glm
glm.cellFitIndC2 = find(u.isC2(fitInd));
glm.cellFitIndL23 = find(u.cellDepths(fitInd) < L4depth);
glm.cellFitDepths = u.cellDepths(fitInd);
glm.cellFitxpoint = u.cellx(fitInd);
glm.cellFitypoint = u.celly(fitInd);
glm.cellNums = u.cellNums;
glm.cellDepths = u.cellDepths;
glm.isC2 = u.isC2;

%% finding tuning
% here, simply assume anything that has coefficient in indPartial{1}, other
% than all touch coefficients (mod(x,8) ~=0)
allFit = cell2mat(glm.cellFunction);
touchFitInd = find(allFit(:,1));
tunedInd = zeros(size(allFit,1),1);
tunedAngle = zeros(size(allFit,1),1);
tuneDirection = zeros(size(allFit,1),1);
tuneCoeffInd = setdiff(indPartial{1}, find(mod(indPartial{1},length(angles)+1)==0));
angleCoeffInd = cell(length(angles),1);
for ai = 1 : length(angles)
    angleCoeffInd{ai} = find(mod(indPartial{1},length(angles)+1) == ai);
end
for ci = 1 : length(touchFitInd)
    ind = touchFitInd(ci);    
    coeff = averageCoeff{ind};
    if ~isempty(find(coeff(tuneCoeffInd+1)))
        tunedInd(ind) = 1;
        coeffSum = zeros(length(angles),1);
        for ai = 1 : length(angles)            
            coeffSum(ai) = sum(coeff(angleCoeffInd{ai}+1)); % either excited or inhibited
        end
        [~, maxInd] = max(abs(coeffSum));
        tunedAngle(ind) = angles(maxInd);
        if coeffSum(maxInd) > 0
            tuneDirection(ind) = 1;
        else
            tuneDirection(ind) = 2;
        end
    end
end
glm.tunedID = glm.cellFitID(find(tunedInd));
glm.tunedAngle = tunedAngle(find(tunedAngle));
glm.tuneDirection = tuneDirection(find(tuneDirection));
glm.touchID = glm.cellFitID(setdiff(touchFitInd, find(tunedInd)));

%% transfer information fron anova results
cadat = load(cafn);
ca.tunedID = cadat.cellsTuned;
ca.tunedAngle = cadat.tuneAngle;
ca.tuneDirection = cadat.tuneDirection;
ca.touch = cadat.cellsNTResponse;

spkdat = load(spkfn);
spk.tunedID = spkdat.cellsTuned;
spk.tunedAngle = spkdat.tuneAngle;
spk.tuneDirection = spkdat.tuneDirection;
spk.touch = spkdat.cellsNTResponse;

