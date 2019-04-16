% for figure 3a
% calcium
% spikes
% full model (without whisker touch variables
% devExp (mean)
% correlation
% predictors
%     touch
%     whisking
%         onset
%         amplitude
%         midpoint
%     licking
%     reward
%     pole up sound

%% basic settings
mouse = 25;
session = 4;
repeat = 10;
%% dependent settings
ufn = sprintf('UberJK%03dS%02d',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_m45_R', mouse, session);
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
allDE = zeros(numCells,1);
calcium = cell(numCells,1);
spikes = cell(numCells,1);
fullModel = cell(numCells,1);
predictors = cell(numCells,1);
corrVal = cell(numCells,1); 
parfor ci = 1 : numCells    
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

    % assigning values within parfor loop
    allDE(ci) = devExplained;
    calcium{ci} = dF;
    spikes{ci} = spkTest;
    fullModel{ci} = model;
    predictors{ci} = testInput;
    corrVal{ci} = corr(spkTest', model);
end

