function glmWhisker = glm_results_WKV_exclusion(mouse, session, baseDir)

%% For figure 3. GLM result plots fro function assignment in each cell
%%
% same as in "glm_results_cell_function_shuffling", except that
% touch angles are replaced by whisker touch variables

% Notes:
%     Assume u is fixed. No change of u.cellNums (very important in indexing)
%     After running 10 repeats of glmnet_cell_type_lasso.m
%     After running 10 repeats of glmnet_whisker_lasso_touchCells.m
%     Based on d190322_glm_repetition.m (and former checking codes)
%     Based on glm_results_cell_function_shuffling_WTV_ONLY.m
%         Only using the exclsion method, not shuffling (it takes long and
%         the result is similar)
%
% Input: 
%     mouse
%     session
%     baseDir
% 
% Files:
%     10 repeated lasso results
%     Uber file
% 
% Output:
%     glmWhisker (struct including cID, deviance, deviance explained, DE difference, whisker variables DE difference, exclusion error ratio, and whisker variables error ratio)
% 2019/07/16 JK

%% basic settings
repeat = 10;
glmWhisker = struct;

                            wkvPosShift =  3;
                            wkvNumVar = 12;
                            
                            
%% dependent settings
ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
glmfnBase = sprintf('glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R', mouse, session);

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

%% assigning functions to each cell
deviance = zeros(numCells,1);
devExp = zeros(numCells,1);
DEdiff = cell(numCells,1);
exclusionER = cell(numCells,1);

whiskerVariableExclusionER = zeros(numCells,wkvNumVar);
whiskerVariableDEdiff = zeros(numCells,wkvNumVar);

parfor ci = 1 : numCells
    if ~isempty(averageCoeff{ci})
        fprintf('Processing JK%03d S%02d: %d/%d\n', mouse, session, ci, numCells)
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

        tempPartialDEsub = zeros(1,length(indPartial));
        tempExclusionER = zeros(1,length(indPartial));

        tempWKVDEdiff = zeros(1,wkvNumVar);
        tempWKVexclusionER = zeros(1,wkvNumVar);

        for pi = 1 : length(indPartial)
            %% exclusion method
            partialInds = setdiff(1:length(coeff), indPartial{pi}+1); % including intercept
            partialCoeffs = coeff(partialInds);
            partialModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,partialInds(2:end)-1)]*partialCoeffs);
            partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
            partialDevExp = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            tempPartialDEsub(pi) = devExplained - partialDevExp;            
            tempExclusionER(pi) = (saturatedLogLikelihood - partialLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);

            if pi == 1
                for j = 1 : wkvNumVar
                    partialInds = setdiff(1:length(coeff), indPartial{pi}((j-1)*wkvPosShift+1:j*wkvPosShift) + 1); % including intercept
                    partialCoeffs = coeff(partialInds);
                    partialModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,partialInds(2:end)-1)]*partialCoeffs);
                    partialLogLikelihood = sum(log(poisspdf(spkTest',partialModel)));
                    partialDevExp = (partialLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                    tempWKVDEdiff(j) = devExplained - partialDevExp;            
                    tempWKVexclusionER(j) = (saturatedLogLikelihood - partialLogLikelihood)/(saturatedLogLikelihood - fullLogLikelihood);
                end                         
            end
        end
        whiskerVariableDEdiff(ci,:) = tempWKVDEdiff;
        whiskerVariableExclusionER(ci,:) = tempWKVexclusionER;
        deviance(ci) = tempDeviance;
        devExp(ci) = devExplained;
        DEdiff{ci} = tempPartialDEsub;
        exclusionER{ci} = tempExclusionER;
    end
end

glmWhisker.cID = cIDAll;
glmWhisker.deviance = deviance;
glmWhisker.devExp = devExp;
glmWhisker.DEdiff = DEdiff;
glmWhisker.exclusionER = exclusionER;
glmWhisker.whiskerVariableExclusionER = whiskerVariableExclusionER;
glmWhisker.whiskerVariableDEdiff = whiskerVariableDEdiff;


