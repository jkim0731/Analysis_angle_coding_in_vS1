function traces = get_traces_per_cell(u, cID, allPredictors, coeff, posShift)

tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
planeInd = floor(cID/1000);

testInput = allPredictors{planeInd};
spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));
dF = cell2mat(cellfun(@(x) [nan(1,posShift), x.dF(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

if length(testInput) ~= length(spkTest)
    error('input matrix and spike length mismatch')
end

% finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
% traces.spikes = spkTest(finiteIndTest);
% traces.calcium = dF(finiteIndTest);
% traces.model = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*coeff);
% traces.predictors = testInput(finiteIndTest,:);


traces.spikes = spkTest;
traces.calcium = dF;
traces.model = exp([ones(size(testInput,1),1),testInput]*coeff);
traces.predictors = testInput;
