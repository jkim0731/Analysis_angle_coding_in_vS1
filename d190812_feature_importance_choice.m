%% Feature importance in choice
%% (1) feature distribution (2) individual AUC of ROC curve (3) prediction performance (4) abs coefficients (5) variance importance (drop out method )

%% individual touch - before the first lick
%% individual touch - before answer lick
%% mean touch - before the first lick
%% mean touch - before answer lick

%% naive 7 angles, expert 7 angles
%% (naive 2 angles, expert 2 angles)
%% (radial distance)

%% individual touch - before the first lick, naive 7 angles

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing];
% fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing, 'grouped'];
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';

load([dir, fn])


%% Draw each features


choice = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(choice), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ci = 1 : length(choice)
        tempInd = find(groupMdl{mi}.io.Y == choice(ci));
        for i = 1 : size(features,3)
            features{mi,ci,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end
figure,
index = [1,2,5,6,9,10,3,4,7,8,11,12];
colorList = {'c','m'};
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(choice));
    for ci = 1 : length(choice)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ci) = histcounts(tempCellFeature{mi,ci}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ci));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList{ci}})
    end
    
    
    title(groupMdl{1}.fitCoeffsFields{i})
end
legend({'Right', 'Left'})

%% Calculate ROC curve AUC
AUCval = zeros(length(groupMdl),size(features,3));

for mi = 1 : length(groupMdl)
    for fi = 1 : size(features,3)
        pred = groupMdl{mi}.io.X(:,fi);
        resp = groupMdl{mi}.io.Y;
        mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        [~,~,~,AUCtemp] = perfcurve(resp,scores,1);
        if AUCtemp < 0.5
            AUCval(mi,fi) = 1-AUCtemp;
        else
            AUCval(mi,fi) = AUCtemp;
        end
    end
end

performance = cellfun(@(x) mean(x.gof.modelAccuracy), groupMdl);

figure, hold on
bar([1:6, 8:size(features,3)+1], mean(AUCval),'k')
errorbar([1:6, 8:size(features,3)+1], mean(AUCval), std(AUCval)/sqrt(length(groupMdl)), 'k.')
ylim([0.5 1])
ylabel('AUC')

yyaxis right
bar(size(features,3) + 3, mean(performance), 'facecolor', ones(1,3)*0.4);
errorbar(size(features,3) + 3, mean(performance), std(performance)/sqrt(length(groupMdl)), '.', 'color', ones(1,3)*0.4);
ylim([0.5 1])
ylabel('Performance')
set(gca,'YColor', ones(1,3)*0.4)

xticks([1:6, 8:size(features,3)+1, size(features,3)+3])
xticklabels({groupMdl{1}.fitCoeffsFields{:}, 'Full model'})
xtickangle(45)

%% Show prediction (mean from 10 iterations)

performance = zeros(length(groupMdl), 10);
% confMat = zeros([size(groupMdl{1}.gof.confusionMatrix), length(groupMdl)]);
for mi = 1 : length(groupMdl)    
    performance(mi,:) = groupMdl{mi}.gof.modelAccuracy;
%     confMat(:,:,mi) = groupMdl{mi}.gof.confusionMatrix ./ sum(groupMdl{mi}.gof.confusionMatrix);
end
meanPerformance = round(mean(mean(performance))*100,2)
semPerformance = round(std(mean(performance,2))/sqrt(length(groupMdl))*100,2)


%% Coefficients 

coeffs = zeros([length(groupMdl), size(groupMdl{1}.fitCoeffs,1)]);
for i = 1 : length(groupMdl)    
    coeffs(i,:) = abs(mean(groupMdl{i}.fitCoeffs,2));
end
[~,ITorderIndTemp] = sort(mean(coeffs(:,2:7)), 'descend');
[~,DTorderIndTemp] = sort(mean(coeffs(:,8:end)), 'descend');
ITorderInd = ITorderIndTemp + 1;
DTorderInd = DTorderIndTemp + 7;
orderInd = [1,ITorderInd,DTorderInd];

figure, hold on
bar([1, 3:8, 10:size(coeffs,2)+2], mean(coeffs(:,orderInd)), 'k')
errorbar([1, 3:8, 10:size(coeffs,2)+2], mean(coeffs(:,orderInd)), std(coeffs(:,orderInd))/sqrt(size(coeffs,1)), 'k.')

xticks([1, 3:8, 10:size(coeffs,2)+2])
xticklabels({'Bias', groupMdl{1}.fitCoeffsFields{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylabel('Mean absolute coefficients')


%% Variable importance - exclusion method, deviance explained

nGroup = length(groupMdl);
nIter = size(groupMdl{1}.fitCoeffs,2);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
% for i = 1
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = mean(groupMdl{gi}.fitCoeffs,2);
    
    predMat = dataX * coeffs;    
    predY = exp(predMat)./(1+exp(predMat));
    LL(gi,1) = sum(log(binopdf(dataY,1,predY)));
    nuLL = sum(log(binopdf(dataY,1,0.5)));
    satLL = 0;
    DE(gi,1) = (LL(gi,1)-nuLL)/(satLL-nuLL);
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat)./(1+exp(predMat));
        
        LL(gi,fi+1) = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,fi+1) = (LL(gi,fi+1)-nuLL)/(satLL-nuLL);
    end
end

figure, 
tempMat = (DE(:,1) - DE(:,2:end)) ./ DE(:,1);
tempMat = tempMat(:,orderInd(2:length(orderInd))-1);
bar([1:6, 8:size(DE,2)], mean(tempMat), 'k'), hold on
errorbar([1:6, 8:size(DE,2)], mean(tempMat), std(tempMat)/sqrt(nGroup), 'k.')
xticklabels({groupMdl{1}.fitCoeffsFields{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylabel('Variable importance')
box('off')


figure, bar([1, 3:8, 10:size(DE,2)+2], mean(DE(:,orderInd)), 'k'), hold on
errorbar([1, 3:8, 10:size(DE,2)+2], mean(DE(:,orderInd)), std(DE(:,orderInd))/sqrt(nGroup), 'k.')
ylabel('Deviance explained')
tempLabelCell = cellfun(@(x) ['-',x], groupMdl{1}.fitCoeffsFields, 'uniformoutput', false);
xticks([1, 3:8, 10:size(DE,2)+2])
xticklabels({'Full model', tempLabelCell{orderInd(2:length(orderInd))-1}})
xtickangle(45)


%% variable necessity test
%% re-training w/o each variable

coeffNames = {'wo_dTheta','wo_dPhi','wo_dKh','wo_dKv','wo_slideDistance','wo_duration','wo_theta','wo_phi','wo_Kh','wo_Kv','wo_radialDistance','wo_count'};
DE = zeros(length(groupMdl), length(coeffNames));

for ci = 1 : length(coeffNames)
    fnCoeff = [fn, '_', coeffNames{ci}];
    singleCoeffData = load([dir, fnCoeff]);
    for gi = 1 : length(singleCoeffData.groupMdl)
        dataX = [ones(size(singleCoeffData.groupMdl{gi}.io.X,1),1),singleCoeffData.groupMdl{gi}.io.X];
        dataY = singleCoeffData.groupMdl{gi}.io.Y;
        coeffs = mean(singleCoeffData.groupMdl{gi}.fitCoeffs,2);

        predMat = dataX * coeffs;    
        predY = exp(predMat)./(1+exp(predMat));
        LL = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,ci) = (LL-nuLL)/(satLL-nuLL);
    end
end

fullDE = zeros(length(groupMdl),1);
for gi = 1 : length(groupMdl)
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    coeffs = mean(groupMdl{gi}.fitCoeffs,2);

    predMat = dataX * coeffs;    
    predY = exp(predMat)./(1+exp(predMat));
    LL = sum(log(binopdf(dataY,1,predY)));
    nuLL = sum(log(binopdf(dataY,1,0.5)));
    satLL = 0;
    fullDE(gi) = (LL-nuLL)/(satLL-nuLL);
end

necessityLevel = (fullDE - DE) ./ fullDE ;

figure, hold on
tempMat = necessityLevel(:,orderInd(2:length(orderInd))-1);
bar([1:6, 8:size(tempMat,2)+1], mean(tempMat), 'k')
errorbar([1:6, 8:size(tempMat,2)+1], mean(tempMat), std(tempMat)/sqrt(size(tempMat,1)), 'k+')
xticks([1:6, 8:size(tempMat,2)+1])
nameList = {groupMdl{1}.fitCoeffsFields{orderInd(2:length(orderInd))-1}};
xticklabels(nameList)
xtickangle(45)
ylabel('Necessity level')


%% Comparing performance between full model and partial model. Until the top 2. No re-training
%% From previous results, use (1) slide distance + dKv, (2) slide distance, and (3) dKv.
%% Compare these along with the full model
%% The order of variables are (1) 4 & 5, (2) 5, and (3) 4.

nGroup = length(groupMdl);
indGroups = {[2],[4],[2,4],[1:12]};
LL = zeros(nGroup, length(indGroups));
DE = zeros(nGroup, length(indGroups));
performance = zeros(nGroup, length(indGroups));
angles = unique(groupMdl{1}.io.Y);
for gi = 1 : nGroup
% for i = 1
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    coeffs = mean(cell2mat(cellfun(@(x) mean(x.fitCoeffs,2), groupMdl,'uniformoutput', false)'),2);
    
    for ii = 1 : length(indGroups)
        tempDataX = dataX(:,[1, indGroups{ii}+1]);
        tempCoeffs = coeffs([1, indGroups{ii}+1]);
        predMat = tempDataX * tempCoeffs;
        
        predY = exp(predMat)./(1+exp(predMat));
        LL = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,ii) = (LL-nuLL)/(satLL-nuLL);
        performance(gi,ii) = length(find((predY>0.5) - dataY == 0)) / length(dataY) * 100;
    end
end
pDE = zeros(1,length(indGroups)-1);
pPerformance = zeros(1,length(indGroups)-1);

for ii = 1 : length(indGroups)-1
    [~,pDE(ii)] = ttest(DE(:,ii), DE(:,end));
    [~,pPerformance(ii)] = ttest(performance(:,ii), performance(:,end));
end

figure, hold on
for i = 1 : length(groupMdl)
    plot(1:length(indGroups), DE(i,:), 'ko-','MarkerFaceColor','black')
end
ylabel('Deviance explained')
xticks(1:length(indGroups))
xticklabels({'\Delta\phi', '\Delta\kappa_V', '\Delta\phi + \Delta\kappa_V', 'Full model'})
xtickangle(45)

figure, hold on
for i = 1 : length(groupMdl)
    plot(1:length(indGroups), performance(i,:), 'ko-','MarkerFaceColor','black')
end
xticks(1:length(indGroups))
xticklabels({'\Delta\phi', '\Delta\kappa_V', '\Delta\phi + \Delta\kappa_V', 'Full model'})
xtickangle(45)
ylabel('Performance (%)')




%% Comparing performance between full model and partial model. Until the top 2. From the re-training
%% From previous results, use (1) slide distance + dKv, (2) slide distance, and (3) dKv.
%% Compare these along with the full model
%% The order of variables are (1) 4 & 5, (2) 5, and (3) 4.
fnList = {'_dPhi', '_dKv','_top2'};
performance = zeros(length(groupMdl), length(fnList)+1);
for fi = 1 : length(fnList)
    result = load([fn, fnList{fi}]);    
    for mi = 1 : length(groupMdl)        
        performance(mi,fi) = mean(result.groupMdl{mi}.gof.modelAccuracy);        
    end
end
performance(:,end) = cellfun(@(x) mean(x.gof.modelAccuracy), groupMdl);

pVal = zeros(1,length(fnList));
for fi = 1 : length(fnList)
    [~, pVal(fi)] = ttest(performance(:,fi), performance(:,end));
end

figure, hold on
for mi = 1 : length(groupMdl)
    plot(performance(mi,:), 'ko-')
end
xlim([0.5 4.5])
xticks([1:4])
xticklabels({'\Delta\phi', '\Delta\kappa_V', '\Delta\phi + \Delta\kappa_V', 'Full model'})

xtickangle(45)
ylabel('Performance (%)')























%% for 90 degrees only


%% Draw each features


choice = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(choice), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    outcomeMat = groupMdl{mi}.outcomes.matrix;
    ind = intersect(find(outcomeMat(5,:)>-1), find(outcomeMat(7,:)));
    ind90 = find(outcomeMat(6,ind) == 90);
    for ci = 1 : length(choice)
        tempInd = find(groupMdl{mi}.io.Y(ind90) == choice(ci));
        for i = 1 : size(features,3)
            features{mi,ci,i} = groupMdl{mi}.io.X(ind90(tempInd),i);
        end
    end
end
figure,
colorList = {'c','m'};
for i = 1 : size(features,3)
    subplot(3,4,i), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(choice));
    for ci = 1 : length(choice)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ci) = histcounts(tempCellFeature{mi,ci}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ci));
        
        nonanLength = size(tempMat,1) - length(find(isnan(sum(tempMat,2))));
        shadedErrorBar(histRange(1:end-1), nanmean(tempMat), nanstd(tempMat)/sqrt(nonanLength), 'lineprop', {'color', colorList{ci}})
    end
    
    
    title(groupMdl{1}.fitCoeffsFields{i})
end
legend({'Right', 'Left'})

%% Calculate ROC curve AUC
AUCval = nan(length(groupMdl),size(features,3));

for mi = 1 : length(groupMdl)    
    outcomeMat = groupMdl{mi}.outcomes.matrix;
    ind = intersect(find(outcomeMat(5,:)>-1), find(outcomeMat(7,:)));
    ind90 = find(outcomeMat(6,ind) == 90);
    for fi = 1 : size(features,3)
        pred = groupMdl{mi}.io.X(ind90,fi);
        resp = groupMdl{mi}.io.Y(ind90);
        if length(unique(resp)) == 2
            mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
            scores = mdl.Fitted.Probability;

            [~,~,~,AUCtemp] = perfcurve(resp,scores,1);
            if AUCtemp < 0.5
                AUCval(mi,fi) = 1-AUCtemp;
            else
                AUCval(mi,fi) = AUCtemp;
            end
        end
    end
end

figure, hold on
nonanLength = size(AUCval,1) - length(find(isnan(sum(AUCval,2))));
bar(1:size(features,3), nanmean(AUCval),'k')
errorbar(1:size(features,3), nanmean(AUCval), nanstd(AUCval)/sqrt(nonanLength), 'k.')

xticks(1:size(features,3))
xticklabels(groupMdl{1}.fitCoeffsFields)
xtickangle(45)
ylim([0.5 1])
ylabel('AUC')



%% for 90 degrees only
%% 

Xhow = ''; 
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'answer'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing];
% fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing, 'grouped'];
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';

load([dir, fn])




%% Show prediction (mean from 10 iterations)

Xhow = 'Mean'; 
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'answer'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing, '_90degrees'];
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
load([dir, fn])

performance = zeros(length(groupMdl), 10);
% confMat = zeros([size(groupMdl{1}.gof.confusionMatrix), length(groupMdl)]);
for mi = 1 : length(groupMdl)    
    performance(mi,:) = groupMdl{mi}.gof.modelAccuracy;
%     confMat(:,:,mi) = groupMdl{mi}.gof.confusionMatrix ./ sum(groupMdl{mi}.gof.confusionMatrix);
end
meanPerformance = round(mean(mean(performance))*100,2)
semPerformance = round(std(mean(performance,2))/sqrt(length(groupMdl))*100,2)


%% Coefficients 

coeffs = zeros([length(groupMdl), size(groupMdl{1}.fitCoeffs,1)]);
for i = 1 : length(groupMdl)    
    coeffs(i,:) = abs(mean(groupMdl{i}.fitCoeffs,2));
end
    
figure, hold on
bar(1:size(coeffs,2), mean(coeffs), 'k')
errorbar(1:size(coeffs,2), mean(coeffs), std(coeffs)/sqrt(size(coeffs,1)), 'k.')

xticks(1:size(coeffs,2))
xticklabels({'Bias', groupMdl{1}.fitCoeffsFields{:}})
xtickangle(45)
ylabel('Mean absolute coefficients')


%% Variable importance - exclusion method, deviance explained

nGroup = length(groupMdl);
nIter = size(groupMdl{1}.fitCoeffs,2);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
% for i = 1
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = mean(groupMdl{gi}.fitCoeffs,2);
    
    predMat = dataX * coeffs;    
    predY = exp(predMat)./(1+exp(predMat));
    LL(gi,1) = sum(log(binopdf(dataY,1,predY)));
    nuLL = sum(log(binopdf(dataY,1,0.5)));
    satLL = 0;
    DE(gi,1) = (LL(gi,1)-nuLL)/(satLL-nuLL);
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat)./(1+exp(predMat));
        
        LL(gi,fi+1) = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,fi+1) = (LL(gi,fi+1)-nuLL)/(satLL-nuLL);
    end
end

figure, bar(0:nFeatures, mean(DE), 'k'), hold on
errorbar(0:nFeatures, mean(DE), std(DE)/sqrt(nGroup), 'k.')
ylabel('Deviance explained')
tempLabelCell = cellfun(@(x) ['-',x], groupMdl{1}.fitCoeffsFields, 'uniformoutput', false);
xticks(0:nFeatures)
xticklabels({'Full model', tempLabelCell{:}})
xtickangle(45)

figure, 
tempMat = (DE(:,1) - DE(:,2:end)) ./ DE(:,1);
bar(1:nFeatures, mean(tempMat), 'k'), hold on
errorbar(1:nFeatures, mean(tempMat), std(tempMat)/sqrt(nGroup), 'k.')
xticklabels(groupMdl{1}.fitCoeffsFields)
xtickangle(45)
ylabel('Variable importance')
box('off')
