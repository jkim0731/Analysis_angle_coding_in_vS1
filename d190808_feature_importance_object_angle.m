%% Feature importance in object angle prediction
%% (1) correlation (2) prediction performance (3) abs coefficients (4) drop out method (5) use one method
%% individual touch - before the first lick
%% individual touch - before answer lick
%% mean touch - before the first lick
%% mean touch - before answer lick

%% naive 7 angles, expert 7 angles
%% (naive 2 angles, expert 2 angles)
%% (radial distance)







%% for 7 angles (discrete)

%% individual touch - before the first lick, naive 7 angles

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'answer'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing];
% fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing, 'grouped'];
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';

load([dir, fn])


%% Draw each features
% for mi = 1 : length(groupMdl)
%     figure, 
%     for i = 1 : length(groupMdl{mi}.fitCoeffsFields)
%         subplot(3,4,i), scatter(groupMdl{mi}.io.Y, groupMdl{mi}.io.X(:,i), 'k.')
%         title(groupMdl{mi}.fitCoeffsFields{i})
%     end
% end
angles = unique(groupMdl{1}.io.Y);
features = zeros(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features(mi,ai,i) = mean(groupMdl{mi}.io.X(tempInd,i));
        end
    end
end
figure,
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i))
    tempMat = squeeze(features(:,:,i));
    shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)))
    title(groupMdl{1}.fitCoeffsFields{i})
    if length(angles) > 2
        xticks(1:2:length(angles))
        xticklabels(angles(1:2:end))
    else
        xticks(1:length(angles))
        xticklabels([45,135])
    end
    if i == 3
        xlabel('Object angle (\circ)')
    end
end
%% Calculate correlation
% corrVal = zeros(length(groupMdl), length(groupMdl{1}.fitCoeffsFields));
% for mi = 1 : length(groupMdl)
%     for i = 1 : length(groupMdl{mi}.fitCoeffsFields)
%         corrVal(mi,i) = corr(groupMdl{mi}.io.Y, groupMdl{mi}.io.X(:,i));
%     end
% end
% figure, hold on
% bar(1:size(corrVal,2), mean(corrVal), 'k')
% errorbar(1:size(corrVal,2), mean(corrVal), std(corrVal)/sqrt(size(corrVal,1)), 'k.')
% xticks(1:size(corrVal,2))
% xticklabels(groupMdl{1}.fitCoeffsFields)
% xtickangle(45)
% ylabel('Correlation')

%% Instead of calculating correlation, it makes much more sense to calculate classification performance from each variable only
% I only have mean touch results for this
if strcmpi(Xhow, 'individual')
    coeffNames = {'dTheta','dPhi','dKh','dKv','slideDistance','duration','theta','phi','Kh','Kv','radialDistance'};
else
    coeffNames = {'dTheta','dPhi','dKh','dKv','slideDistance','duration','theta','phi','Kh','Kv','radialDistance','count'};
end
coeffAccuracy = zeros(length(groupMdl), length(coeffNames)+1);
for ci = 1 : length(coeffNames)
    fnCoeff = [fn, '_', coeffNames{ci}];
    singleCoeffData = load([dir, fnCoeff]);
    coeffAccuracy(:,ci) = cellfun(@(x) mean(x.gof.modelAccuracy), singleCoeffData.groupMdl');
end
coeffAccuracy(:,end) = cellfun(@(x) mean(x.gof.modelAccuracy), groupMdl');
xpoints = [1:6, 8:length(coeffNames)+1, length(coeffNames)+3];
figure, hold on
bar(xpoints, mean(coeffAccuracy), 'k')
errorbar(xpoints, mean(coeffAccuracy), std(coeffAccuracy)/sqrt(size(coeffAccuracy,1)), 'k+')
xticks(xpoints)
xticklabels({coeffNames{:}, 'All'})
xtickangle(45)
ylabel('Classification performance')
xlimrange = [0 length(coeffNames)+4];
xlim(xlimrange)
plot(xlimrange , ones(1,length(xlimrange))*1/7, '--', 'color', ones(1,3)*0.4)
%% Show prediction (mean from 10 iterations)
%% For 7 angles
performance = zeros(length(groupMdl), 10);
% confMat = zeros([size(groupMdl{1}.gof.confusionMatrix), length(groupMdl)]);
for mi = 1 : length(groupMdl)    
    performance(mi,:) = groupMdl{mi}.gof.modelAccuracy;
%     confMat(:,:,mi) = groupMdl{mi}.gof.confusionMatrix ./ sum(groupMdl{mi}.gof.confusionMatrix);
end
meanPerformance = mean(mean(performance))
semPerformance = std(mean(performance,2))/sqrt(length(groupMdl))

% confusion matrix
figure, imagesc(mean(confMat,3)), colorbar, 
title([sprintf('%.2f',meanPerformance*100), ' \pm ', sprintf('%.2f',semPerformance*100), ' %'])
xlabel('Prediction'), ylabel('Data')
xticklabels(unique(groupMdl{1}.io.Y))
yticklabels(unique(groupMdl{1}.io.Y))

%% Coefficients (for 7 angles)

coeffs = zeros([size(groupMdl{1}.fitCoeffs{1}), length(groupMdl)]);
for i = 1 : length(groupMdl)
    tempCoeff = zeros(size(groupMdl{i}.fitCoeffs{1}));
    for j = 1 : 10
        tempCoeff = tempCoeff + groupMdl{i}.fitCoeffs{j};
    end
    coeffs(:,:,i) = tempCoeff / 10;
end

absCoeffs = zeros([size(groupMdl{1}.fitCoeffs{1},1), length(groupMdl)]);
for i = 1 : length(groupMdl)
    absCoeffs(:,i) = mean(abs(groupMdl{i}.fitCoeffs{1}),2);
end

[~,ITorderIndTemp] = sort(mean(absCoeffs(2:7,:),2), 'descend');
[~,DTorderIndTemp] = sort(mean(absCoeffs(8:end,:),2), 'descend');
ITorderInd = ITorderIndTemp + 1;
DTorderInd = DTorderIndTemp + 7;
orderInd = [1;ITorderInd;DTorderInd];

figure, hold on
tempColor = jet(7);
if size(groupMdl{1}.fitCoeffs{1},2) == 2
    colorList = tempColor([1,7],:);
else
    colorList = tempColor;
end
% for i = 1 : size(groupMdl{1}.fitCoeffs{1},2)
%     tempMat = squeeze(coeffs(:,i,:));
%     plot(1:size(tempMat,1), mean(tempMat,2), 'o', 'color', colorList(i,:))
% end

% for i = 1 : size(groupMdl{1}.fitCoeffs{1},2)
%     tempMat = squeeze(coeffs(:,i,:));    
%     shadedErrorBar(1:size(tempMat,1), mean(tempMat,2), std(tempMat,[],2)/sqrt(size(tempMat,2)), 'lineProp',{'color', colorList(i,:)})
% end
for i = 1 : size(groupMdl{1}.fitCoeffs{1},2)
    tempMat = squeeze(coeffs(orderInd,i,:));
    errorbar([1, 3:8, 10:size(tempMat,1)+2], mean(tempMat,2), std(tempMat,[],2)/sqrt(size(tempMat,2)), 'o', 'color', colorList(i,:))
end

% for i = 1 : size(groupMdl{1}.fitCoeffs{1},2)
%     tempMat = squeeze(coeffs(:,i,:));
%     plot(1:size(tempMat,1), mean(tempMat,2), 'color', colorList(i,:))
% end

xticks([1, 3:8, 10:size(tempMat,1)+2])
xticklabels({'Bias', groupMdl{1}.fitCoeffsFields{orderInd(2:end)-1}})
xtickangle(45)
ylabel('Coefficients')
angles = unique(groupMdl{1}.io.Y);
legendList = cell(1,length(angles));
for i = 1 : length(legendList)
    legendList{i} = [num2str(angles(i)), '\circ'];
end
legend(legendList)

figure, hold on
tempMat = absCoeffs(orderInd,:)';
bar([1, 3:8, 10:size(tempMat,2)+2], mean(tempMat), 'k')
errorbar([1, 3:8, 10:size(absCoeffs,1)+2], mean(tempMat), std(tempMat)/sqrt(size(tempMat,1)), 'k.')
xticks([1, 3:8, 10:size(absCoeffs,1)+2])
xticklabels({'Bias', groupMdl{1}.fitCoeffsFields{orderInd(2:end)-1}})
xtickangle(45)
ylabel('Mean absolute coefficients')


%% Variable importance - exclusion method

nGroup = length(groupMdl);
nIter = length(groupMdl{1}.fitCoeffs);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
% for i = 1
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = zeros(size(groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
%     for ii = 1
        coeffs = coeffs + groupMdl{gi}.fitCoeffs{ii}/nIter;
    end
    
    predMat = dataX * coeffs;
    predY = exp(predMat) ./ sum(exp(predMat),2);
    if find(abs(sum(predY,2)-1)>10^(-10))
        error('prob does not sum to 1')
    end
    
    patternY = zeros(length(dataY),length(listOfY));
    for pi = 1 : size(patternY,1)
        patternY(pi,find(angles==dataY(pi))) = 1;
    end
    
    nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
    satLL = 0;    
    LL(gi,1) = sum(log(mnpdf(patternY,predY)));    
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat) ./ sum(exp(predMat),2);
        if find(abs(sum(predY,2)-1)>10^(-10))
            error('prob does not sum to 1')
        end
        LL(gi,fi+1) = sum(log(mnpdf(patternY,predY)));
    end
    DE(gi,:) = (LL(gi,:) - nuLL)/(satLL-nuLL);
end

figure, 
tempMat = (DE(:,1)-DE(:,2:end)) ./ DE(:,1);
tempMat = tempMat(:,orderInd(2:length(orderInd))-1);
bar([1:6, 8:nFeatures+1], mean(tempMat), 'k'), hold on
errorbar([1:6, 8:nFeatures+1], mean(tempMat), std(tempMat)/sqrt(nGroup), 'k.')
xticks([1:6, 8:nFeatures+1])
xticklabels({groupMdl{1}.fitCoeffsFields{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylabel('Variable importance')
box('off')

tempMat = DE(:,orderInd);
figure, 
bar([1, 3:8, 10:size(tempMat,2)+2], mean(tempMat), 'k'), hold on
errorbar([1, 3:8, 10:size(tempMat,2)+2], mean(tempMat), std(tempMat)/sqrt(nGroup), 'k.')
ylabel('Deviance explained')
tempLabelCell = cellfun(@(x) ['-',x], groupMdl{1}.fitCoeffsFields, 'uniformoutput', false);
xticks([1, 3:8, 10:size(tempMat,2)+2])
xticklabels({'Full model', tempLabelCell{orderInd(2:length(orderInd))-1}})
xtickangle(45)

%% variable necessity test
%% re-training w/o each variable

coeffNames = {'wo_dTheta','wo_dPhi','wo_dKh','wo_dKv','wo_slideDistance','wo_duration','wo_theta','wo_phi','wo_Kh','wo_Kv','wo_radialDistance','wo_count'};
DE = zeros(length(groupMdl), length(coeffNames));
angles = unique(groupMdl{1}.io.Y);
nIter = length(groupMdl{1}.fitCoeffs);
for ci = 1 : length(coeffNames)
    fnCoeff = [fn, '_', coeffNames{ci}];
    singleCoeffData = load([dir, fnCoeff]);
    for gi = 1 : length(singleCoeffData.groupMdl)
        dataX = [ones(size(singleCoeffData.groupMdl{gi}.io.X,1),1),singleCoeffData.groupMdl{gi}.io.X];
        dataY = singleCoeffData.groupMdl{gi}.io.Y;
        listOfY = unique(dataY);
        coeffs = zeros(size(singleCoeffData.groupMdl{gi}.fitCoeffs{1}));
        for ii = 1 : nIter
    %     for ii = 1
            coeffs = coeffs + singleCoeffData.groupMdl{gi}.fitCoeffs{ii}/nIter;
        end

        predMat = dataX * coeffs;
        predY = exp(predMat) ./ sum(exp(predMat),2);
        if find(abs(sum(predY,2)-1)>10^(-10))
            error('prob does not sum to 1')
        end

        patternY = zeros(length(dataY),length(listOfY));
        for pi = 1 : size(patternY,1)
            patternY(pi,find(angles==dataY(pi))) = 1;
        end

        nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
        satLL = 0;    
        LL = sum(log(mnpdf(patternY,predY)));
        DE(gi,ci) = (LL - nuLL)/(satLL-nuLL);
    end
end

fullDE = zeros(length(groupMdl),1);
for gi = 1 : length(groupMdl)
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = zeros(size(groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
%     for ii = 1
        coeffs = coeffs + groupMdl{gi}.fitCoeffs{ii}/nIter;
    end
    
    predMat = dataX * coeffs;
    predY = exp(predMat) ./ sum(exp(predMat),2);
    if find(abs(sum(predY,2)-1)>10^(-10))
        error('prob does not sum to 1')
    end
    
    patternY = zeros(length(dataY),length(listOfY));
    for pi = 1 : size(patternY,1)
        patternY(pi,find(angles==dataY(pi))) = 1;
    end
    
    nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
    satLL = 0;    
    LL = sum(log(mnpdf(patternY,predY)));    
    fullDE(gi,1) = (LL - nuLL)/(satLL-nuLL);
end

necessityLevel = (fullDE - DE) ./ fullDE ;
tempMat = necessityLevel(:,orderInd(2:length(orderInd))-1);
figure, hold on
bar([1:6, 8:size(tempMat,2)+1], mean(tempMat), 'k')
errorbar([1:6, 8:size(tempMat,2)+1], mean(tempMat), std(tempMat)/sqrt(size(tempMat,1)), 'k+')
xticks([1:6, 8:size(tempMat,2)+1])
xticklabels({groupMdl{1}.fitCoeffsFields{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylabel('Necessity level')

%% Comparing performance between full model and partial model. Until the top 2. No re-training
%% From previous results, use (1) slide distance + dKv, (2) slide distance, and (3) dKv.
%% Compare these along with the full model
%% The order of variables are (1) 4 & 5, (2) 5, and (3) 4.

nGroup = length(groupMdl);
nIter = length(groupMdl{1}.fitCoeffs);
indGroups = {[4],[5],[4,5],[3:5],[2,4,5],[2:5],[1:6],[1:12]};
LL = zeros(nGroup, length(indGroups));
DE = zeros(nGroup, length(indGroups));
performance = zeros(nGroup, length(indGroups));
angles = unique(groupMdl{1}.io.Y);
for gi = 1 : nGroup
% for i = 1
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = zeros(size(groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
%     for ii = 1
        coeffs = coeffs + groupMdl{gi}.fitCoeffs{ii}/nIter;
    end
    
    predMat = dataX * coeffs;
    predY = exp(predMat) ./ sum(exp(predMat),2);
    if find(abs(sum(predY,2)-1)>10^(-10))
        error('prob does not sum to 1')
    end
    
    patternY = zeros(length(dataY),length(listOfY));
    for pi = 1 : size(patternY,1)
        patternY(pi,find(angles==dataY(pi))) = 1;
    end
    
    [~, dataYind] = max(patternY,[],2);
    
    nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
    satLL = 0;    
    LL(gi,1) = sum(log(mnpdf(patternY,predY)));    
    for ii = 1 : length(indGroups)
        tempDataX = dataX(:,[1, indGroups{ii}+1]);
        tempCoeffs = coeffs([1, indGroups{ii}+1],:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat) ./ sum(exp(predMat),2);
        if find(abs(sum(predY,2)-1)>10^(-10))
            error('prob does not sum to 1')
        end
        LL(gi,ii) = sum(log(mnpdf(patternY,predY)));
        [~, predYind] = max(predY, [], 2);
        performance(gi,ii) = length(find(predYind - dataYind == 0)) / length(dataYind) * 100;
    end
    DE(gi,:) = (LL(gi,:) - nuLL)/(satLL-nuLL);
end

figure, bar(1:length(indGroups), mean(DE), 'k'), hold on
errorbar(1:length(indGroups), mean(DE), std(DE)/sqrt(length(indGroups)), 'k.')
ylabel('Deviance explained')
xticks(1:length(indGroups))
% xticklabels({'\Delta\kappa_V', 'Slide distance', '\Delta\kappa_V + slide distance', 'Full model'})
xticklabels({'\Delta\kappa_V', 'Slide distance', 'Top2', 'Top3_1', 'Top3_2', 'Top4', 'All kinetics', 'Full model'})
xtickangle(45)

figure, 
bar(1:length(indGroups), mean(performance), 'k'), hold on
errorbar(1:length(indGroups), mean(performance), std(performance)/sqrt(length(indGroups)), 'k.')
% xticklabels({'\Delta\kappa_V', 'Slide distance', '\Delta\kappa_V + slide distance', 'Full model'})
xticklabels({'\Delta\kappa_V', 'Slide distance', 'Top2', 'Top3_1', 'Top3_2', 'Top4', 'All kinetics', 'Full model'})
xtickangle(45)
ylabel('Performance (%)')




%% Comparing performance between full model and partial model. Until the top 2. From the re-training
%% From previous results, use (1) slide distance + dKv, (2) slide distance, and (3) dKv.
%% Compare these along with the full model
%% The order of variables are (1) 4 & 5, (2) 5, and (3) 4.
fnList = {'_top2', '_top3_1', '_top3_2', '_top4', '_kinetics'};
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
xlim([0.5 6.5])
xticks([1:6])
xticklabels({'Top2', 'Top3_1', 'Top3_2', 'Top4', 'All kinetics', 'Full model'})

xtickangle(45)
ylabel('Performance (%)')












%% for 2 angles (Two or RadialDistance)

%% individual touch - before the first lick, naive 7 angles

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'answer'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing];
% fn = ['mdl', task, learned, Xhow, Yout, '_12features_', timing, 'grouped'];
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';

load([dir, fn])

%% Draw each features

angles = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features{mi,ai,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end
figure,

tempColorList = jet(7);
colorList = tempColorList([1,7],:);
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(angles));
    for ai = 1 : length(angles)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ai) = histcounts(tempCellFeature{mi,ai}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ai));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList(ai,:)})
    end
    
    
    title(groupMdl{1}.fitCoeffsFields{i})
end
legend({'45\circ', '135\circ'})




%% Compare ROC curve AUC and single parameter binomial classification

featNames = {'dTheta','dPhi','dKh','dKv','slideDistance','duration','theta','phi','Kh','Kv','radialDistance', 'count'};
accuracy = zeros(6, length(featNames));
for i = 1 : length(featNames)
    data = load([dir, fn, '_', featNames{i}]);
    accuracy(:,i) = cellfun(@(x) mean(x.gof.modelAccuracy), data.groupMdl)';
end

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

figure
hold on
bar([1:length(featNames)]-0.2, mean(accuracy), 0.3, 'k')
bar([1:length(featNames)]+0.2, mean(AUCval), 0.3, 'w')
errorbar([1:length(featNames)]-0.2, mean(accuracy), std(accuracy), 'k.')
errorbar([1:length(featNames)]+0.2, mean(AUCval), std(AUCval), 'k.')
legend({'Binomial', 'AUC'})

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

grayColor = ones(1,3) * 0.4;
figure; hold on
bar([1:6, 8:size(features,3)+1], mean(AUCval),'k')
errorbar([1:6, 8:size(features,3)+1], mean(AUCval), std(AUCval)/sqrt(length(groupMdl)-1), 'k.')
ylim([0.5 1])
ylabel('AUC')

yyaxis right
ylim([0.5 1])
bar(size(features,3) + 3, mean(performance), 'facecolor', grayColor)
errorbar(size(features,3) + 3, mean(performance), std(performance)/sqrt(length(groupMdl)-1), '.', 'color', grayColor)
ylabel('Performance')

xticks([1:6, 8:size(features,3)+1, size(features,3)+3])
xticklabels({groupMdl{1}.fitCoeffsFields{:}, 'Full model'})
xtickangle(45)
set(gca,'YColor', grayColor)






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
