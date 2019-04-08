%% building
% 1) a complete model (from m44) and raw spikes data in each cell
% 2) coefficients for angle response (for ANOVA test later on)
% 3) cell ID's that failed at least one loop of glm run
% Save these in a single file

baseDir = 'C:\JK\';

% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
mice = 25;
sessions = {[4]};

angles = 45:15:135;

for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d', baseDir, mouse))
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        fnBase = sprintf('glmResponseType_JK%03dS%02d_m44_R', mouse, session);
        load([fnBase, '01'], 'cIDAll', 'indPartial', 'allPredictors', 'testTn', 'trainingTn', 'posShift')
        tuneCoeffInd = setdiff(indPartial{1}, find(mod(indPartial{1},length(angles)+1)==0));
        angleCoeffInd = cell(length(angles),1);
        for ai = 1 : length(angles)
            angleCoeffInd{ai} = find(mod(indPartial{1},length(angles)+1) == ai);
        end

        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load(ufn);
        
        
        numCell = length(cIDAll);
        coefficients = cell(numCell,10);
        allDone = zeros(numCell,10);
        for ri = 1 : 10
            fn = sprintf('%s%02d',fnBase, ri);
            load(fn, 'fitCoeffs', 'done')
            coefficients(:,ri) = fitCoeffs;
            allDone(:,ri) = done;
        end
        averageCoeff = zeros(numCell, indPartial{end}(end)+1);
        model = cell(numCell,1);
        spikes = cell(numCell,1);
        angleCoeffs = cell(numCell, 1);
        failed = zeros(numCell,1);
        parfor ci = 1 : numCell
            fprintf('Processing %d/%d of JK%03d S%02d\n', ci, numCell, mouse, session)
            tempCoeff = mean(cell2mat(coefficients(ci,:)),2);
            tempCoeff(abs(tempCoeff)<0.01) = deal(0);
            averageCoeff(ci,:) = tempCoeff';
            
            cID = cIDAll(ci);
            tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
            
            if sum(ismember(u.trialNums(tindCell), union(testTn{ci}, trainingTn{ci}))) < length(tindCell)
                error('trial numbers mismatch')
            end
            
            cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
            planeInd = floor(cID/1000);

            testInput = allPredictors{planeInd};
            spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

            if length(testInput) ~= length(spkTest)
                error('input matrix and spike length mismatch')
            end

            finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInput,2))));
            spkTest = spkTest(finiteIndTest);
            spikes{ci} = (spkTest - min(spkTest)) / (max(spkTest) - min(spkTest));

            coeff = averageCoeff(ci,:);
            tempModel = exp([ones(length(finiteIndTest),1),testInput(finiteIndTest,:)]*coeff');
            model{ci} = (tempModel - min(tempModel)) / (max(tempModel) - min(tempModel));
            
            tempAngleCoeff = zeros(10,length(angles));
            for ai = 1 : length(angles)
                for ri = 1 : 10
                    tempAngleCoeff(ri,ai) = sum(coefficients{ci,ri}(angleCoeffInd{ai}+1));
                end
            end
            angleCoeffs{ci} = tempAngleCoeff;
            if length(find(allDone(:) == ci)) < 10
                failed(ci) = ci;
            end
        end
        failed = unique(failed);
    end
end
            

%%
ci = 39;
a = zeros(1,length(angles));
for ai = 1 : length(angles)
    a(ai) = sum(averageCoeff(ci,angleCoeffInd{ai}+1));
end
a
