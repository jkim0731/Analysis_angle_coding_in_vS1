function devExp = glm_results_dev_exp_power(mouse, session, baseDir)

%% For figure 4. Explanatory power of pole angle and whisker touch variables
%% Also suits for other general purposes
% Notes:
%     Assume u is fixed. No change of u.cellNums (very important in indexing)
%     After running 10 repeats of glm.
%
% Input: 
%     mouse
%     session
%     baseDir
% 
% Files:
%     10 repeated glm results, with both touch and 
%     Uber file
% 
% Output:
%     Deviance explained, from all coefficients, and from each group of coefficients.
%     Other information: calcium, spikes, full model, predictors, cellID, average
%     coefficients, correlation between full model and spikes.
%
%     devExp.allDE
%     devExp.averageDE
%     devExp.partialSub{}
%     devExp.cellID
%     devExp.coeffs{} (including the intercept)
%     devExp.corrVal

% 2019/04/11 JK

%% basic settings
repeat = 10;
devExp = struct;
%% dependent settings
ufn = sprintf('UberJK%03dS%02d',mouse, session);
% glmfnBase = sprintf('glmResponseType_JK%03dS%02d_lasso_R', mouse, session);
glmfnBase = sprintf('glmWhisker_lasso_touchCell_JK%03dS%02d_R', mouse, session);
%% load uber
cd(sprintf('%s%03d',baseDir, mouse))
load(ufn, 'u') % loading u
u = u;
%% gather information from 10 repeat glm results
load(sprintf('%s%02d',glmfnBase, 1), 'allPredictors', 'indPartial', 'posShift', 'cIDAll')
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

%% calculating deviance explained
% before going into parfor...
allDE = nan(numCells,1);
partialDEsub = cell(numCells,1);
corrVal = nan(numCells,1); 
parfor ci = 1 : numCells
    fprintf('Processing JK%03d S%02d %d/%d\n', mouse, session, ci, numCells)
    coeff = averageCoeff{ci};
    if ~isempty(coeff)
        cID = cIDAll(ci);
        tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
        cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
        planeInd = floor(cID/1000);

        testInput = allPredictors{planeInd};
        spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));
        dF = cell2mat(cellfun(@(x) [nan(1,posShift), x.dF(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

        if length(testInput) ~= length(spkTest)
            error('input matrix and spike length mismatch')
        end

        finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
        spkTest = spkTest(finiteIndTest);
        dF = dF(finiteIndTest);

        coeff = averageCoeff{ci};
        model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*coeff);
        mu = mean(spkTest); % null poisson parameter
        nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
        saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
        fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
        devExplained = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);

        tempPartialDEsub = zeros(1,length(indPartial));
        for pi = 1 : length(indPartial)
            partialInds = setdiff(1:length(coeff), indPartial{pi}+1); % including intercept
            partialCoeffs = coeff(partialInds);
            partialModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,partialInds(2:end)-1)]*partialCoeffs);
            partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
            partialDevExp = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            tempPartialDEsub(pi) = devExplained - partialDevExp;
        end

        % assigning values within parfor loop
        allDE(ci) = devExplained;
        partialDEsub{ci} = tempPartialDEsub;
        corrVal(ci) = corr(spkTest', model);
    end
end

%% assigning parfor loop (and also some before) values to output 
devExp.allDE = allDE;
devExp.averageDE = averageDE;
devExp.partialSub = partialDEsub;
devExp.cellID = cIDAll;
devExp.coeffs = averageCoeff;
devExp.corrVal = corrVal;