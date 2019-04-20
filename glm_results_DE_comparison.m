function results = glm_results_DE_comparison(mouse, session, baseDir)

% For d190416_glmDEcomparison_results_summary
% Inputs:
%     mouse, session, baseDir
%
% Files:
%     glmWhiskerTouchVariablesONLY_JKxxxSxx_Rxx
%     glmResponseType_JKxxxSxx_m45_Rxx
%     (use cIDAll from glmWhiskerTouchVariablesONLY files)
%     JKxxxSxxangle_tuning
% 
% Outputs:
%     cIDAll
%     tuned
%     tuned angle
%     complex tuning: multimodal, tuneDirection >= 2, leaveOneOut
%     tuning sharpness
%     ramp tuning
%     touchDE (length(cIDAll), repeat)
%     whiskerDE (length(cIDAll), repeat)
%     touchVSwhiskerRatio: from the average
%     touchDependence: calculated by t-test between whiskerDE and touchDE. if touchDE is statistically larger than whiskerDE.


% 2019/04/18 update
% use whole fitting, instead of 10 repeats. Each repeat can be far off just based on random effects.

% basic setting
repeat = 10;
cd(sprintf('%s%03d',baseDir, mouse))
wtvFnBase = sprintf('glmWhiskerTouchVariablesONLY_JK%03dS%02d_R', mouse, session);
touchFnBase = sprintf('glmResponseType_JK%03dS%02d_m45_R', mouse, session);
tuningFn = sprintf('JK%03dS%02dangle_tuning',mouse,session);
load(tuningFn)
ufn = sprintf('UberJK%03dS%02d',mouse,session);
load(ufn, 'u')
u = u;
% basic output (from touch and tuning information)
results.cIDAll = spk.touchID;
results.tuned = spk.tuned;
results.tunedAngle = spk.tunedAngle;
complex = zeros(length(spk.touchID),1);
complex(union(union(find(spk.multimodal), find(spk.leaveOneOut)), find(spk.tuneDirection >= 2))) = 1;
results.complex = complex;
results.sharpness = spk.sharpness;
results.ramp = spk.ramp;

cIDAll = results.cIDAll;
numCell = length(cIDAll);
touchDE = zeros(numCell, 1);
whiskerDE = zeros(numCell, 1);
touchWhiskerRatio = zeros(numCell, 1);

touchCoeff = cell(numCell, repeat);
wtvCoeff = cell(numCell, repeat);

for ri = 1 : repeat
    wtvFn = sprintf('%s%02d',wtvFnBase, ri);
    wtvdat = load(wtvFn, 'fitCoeffs');
    wtvCoeff(:,ri) = wtvdat.fitCoeffs;
    touchFn = sprintf('%s%02d', touchFnBase, ri);
    tdat = load(touchFn, 'fitCoeffs', 'cIDAll');
    inds = find(ismember(tdat.cIDAll, cIDAll));
    touchCoeff(:,ri) = tdat.fitCoeffs(inds);
end
load(wtvFn, 'allPredictors', 'posShift')
posShift = posShift;
wtvPredictors = allPredictors;
load(touchFn, 'allPredictors')
touchPredictors = allPredictors;

parfor ci = 1 : numCell
    fprintf('Processing %d/%d from JK%03d S%02d\n', ci, numCell, mouse, session);
    cID = cIDAll(ci);

    % spikes
    tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
    spkTemp = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

    planeInd = floor(cID/1000);
    % wtv first
    avgWCoeff = mean(cell2mat(wtvCoeff(ci,:)),2);
    testInput = wtvPredictors{planeInd};

    finiteIndTest = intersect(find(isfinite(spkTemp)), find(isfinite(sum(testInput,2))));
    spkTest = spkTemp(finiteIndTest);

    model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*avgWCoeff);
    mu = mean(spkTest); % null poisson parameter
    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));                    
    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
    
    tempWhiskerDE = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    
    % then touch
    avgTCoeff = mean(cell2mat(touchCoeff(ci,:)),2);
    testInput = touchPredictors{planeInd};
    
    finiteIndTest = intersect(find(isfinite(spkTemp)), find(isfinite(sum(testInput,2))));
    spkTest = spkTemp(finiteIndTest);

    model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*avgTCoeff);
    mu = mean(spkTest); % null poisson parameter
    nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
    saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));                    
    fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
    
    tempTouchDE = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
    

    whiskerDE(ci) = tempWhiskerDE;
    touchDE(ci) = tempTouchDE;
    touchWhiskerRatio(ci) = tempTouchDE / tempWhiskerDE;
end
results.touchDE = touchDE;
results.whiskerDE = whiskerDE;
results.touchWhiskerRatio = touchWhiskerRatio;

