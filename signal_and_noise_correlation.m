%% signal correlation & noise correlation within touch responsive cells
%% from the "full" model, using both touch and WTV

clear
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56, 70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3], [6],[4],[4],[4]};  
naiveInd = 1:12;
matchingNaiveInd = [1:4,7,9];
L4ind = 13:16;

cd(baseDir)
glm = load('glm_results_withWTV');
tune = load('angle_tuning_summary');

sigCorrMap = cell(length(naiveInd),2);
noiseCorrMap = cell(length(naiveInd),2);
% for i = naiveInd
for i = 1    
    mouse = mice(i);
    session = sessions{i}(1);
    load(sprintf('%s%03d\\UberJK%03dS%02d.mat',baseDir,mouse,mouse,session))
    load(sprintf('%s%03d\\glmWithWhiskerTouchVariables_JK%03dS%02d_R01.mat',baseDir,mouse,mouse,session), 'allPredictors', 'posShift')
    
    cIDlist = glm.naive(i).cellID;
    cGroup = cell(2,1);
    cGroup{1} = cIDlist(find(cIDlist<5000));
    cGroup{2} = cIDlist(find(cIDlist>5000));
    for cgi = 1 : 2
        tempSignal = cell(1,length(cGroup{cgi}));
        tempNoise = cell(1,length(cGroup{cgi}));
        for ci = 1 : length(cGroup{cgi})
            cID = cGroup{cgi}(ci);
            tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
            cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
            planeInd = floor(cID/1000);
            
            glmInd = find(glm.naive(i).cellID == cID);
            coeff = glm.naive(i).coeffs{glmInd};
            input = allPredictors{planeInd};
            spike = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));
            model = exp([ones(size(input,1),1),input]*coeff);
            
            tempSignal{ci} = model;
            tempNoise{ci} = min_max_normalization(spike') - min_max_normalization(model);
        end
        temp = tune.naive(i);
        if cgi == 1
            sharpness = temp.sharpness(intersect(find(temp.tuned), find(temp.touchID<5000) ));
            tunedAngle = temp.tunedAngle(intersect(find(temp.tuned), find(temp.touchID<5000) ));
        else
            sharpness = temp.sharpness(intersect(find(temp.tuned), find(temp.touchID>5000) ));
            tunedAngle = temp.tunedAngle(intersect(find(temp.tuned), find(temp.touchID>5000) ));
        end
        
        [~, indsortSharpness] = sort(sharpness, 'descend');
        tempTunedAngle = tunedAngle(indsortSharpness);
        [~, indsortAngle] = sort(tempTunedAngle);
        indsort = indsortSharpness(indsortAngle);

        sigMat = cell2mat(tempSignal);
        sigMat = sigMat(:,indsort);
        finiteInd = find(isfinite(sum(sigMat,2)));
        sigCorrMap{i,cgi} = corr(sigMat(finiteInd,:)).*(1-eye(size(sigMat,2)));
        noiseMat = cell2mat(tempNoise);
        noiseMat = noiseMat(:,indsort);
        finiteInd = find(isfinite(sum(noiseMat,2)));
        noiseCorrMap{i,cgi} = corr(noiseMat(finiteInd,:)).*(1-eye(size(noiseMat,2)));
    end
end
%
figure
subplot(2,2,1), imagesc(sigCorrMap{1,1})
subplot(2,2,2), imagesc(sigCorrMap{1,2})
subplot(2,2,3), imagesc(noiseCorrMap{1,1})
subplot(2,2,4), imagesc(noiseCorrMap{1,2})

%%
sigCorrMap = cell(length(naiveInd),2);
noiseCorrMap = cell(length(naiveInd),2);
% for i = expertInd
for i = 1    
    mouse = mice(i);
    session = sessions{i}(2);
    load(sprintf('%s%03d\\UberJK%03dS%02d.mat',baseDir,mouse,mouse,session))
    load(sprintf('%s%03d\\glmWithWhiskerTouchVariables_JK%03dS%02d_R01.mat',baseDir,mouse,mouse,session), 'allPredictors', 'posShift')
    
    cIDlist = glm.expert(i).cellID;
    cGroup = cell(2,1);
    cGroup{1} = cIDlist(find(cIDlist<5000));
    cGroup{2} = cIDlist(find(cIDlist>5000));
    for cgi = 1 : 2
        tempSignal = cell(1,length(cGroup{cgi}));
        tempNoise = cell(1,length(cGroup{cgi}));
        for ci = 1 : length(cGroup{cgi})
            cID = cGroup{cgi}(ci);
            tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
            cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
            planeInd = floor(cID/1000);
            
            glmInd = find(glm.expert(i).cellID == cID);
            coeff = glm.expert(i).coeffs{glmInd};
            input = allPredictors{planeInd};
            spike = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));
            model = exp([ones(size(input,1),1),input]*coeff);
            
            tempSignal{ci} = model;
            tempNoise{ci} = min_max_normalization(spike') - min_max_normalization(model);
        end
        temp = tune.expert(i);
        if cgi == 1
            sharpness = temp.sharpness(intersect( find(temp.tuned), find(temp.touchID<5000) ));
            tunedAngle = temp.tunedAngle(intersect( find(temp.tuned), find(temp.touchID<5000) ));
            
        else
            sharpness = temp.sharpness(intersect( find(temp.tuned), find(temp.touchID>5000) ));
            tunedAngle = temp.tunedAngle(intersect( find(temp.tuned), find(temp.touchID>5000) ));
        end
        
        [~, indsortSharpness] = sort(sharpness, 'descend');
        tempTunedAngle = tunedAngle(indsortSharpness);
        [~, indsortAngle] = sort(tempTunedAngle);
        indsort = indsortSharpness(indsortAngle);

        sigMat = cell2mat(tempSignal);
        sigMat = sigMat(:,indsort);
        finiteInd = find(isfinite(sum(sigMat,2)));
        sigCorrMap{i,cgi} = corr(sigMat(finiteInd,:)).*(1-eye(size(sigMat,2)));
        noiseMat = cell2mat(tempNoise);
        noiseMat = noiseMat(:,indsort);
        finiteInd = find(isfinite(sum(noiseMat,2)));
        noiseCorrMap{i,cgi} = corr(noiseMat(finiteInd,:)).*(1-eye(size(noiseMat,2)));
    end
end
%
figure
subplot(2,2,1), imagesc(sigCorrMap{1,1})
subplot(2,2,2), imagesc(sigCorrMap{1,2})
subplot(2,2,3), imagesc(noiseCorrMap{1,1})
subplot(2,2,4), imagesc(noiseCorrMap{1,2})