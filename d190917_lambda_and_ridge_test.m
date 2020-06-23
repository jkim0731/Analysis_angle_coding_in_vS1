% saving lambda values from already-run glm models.
% Needs a little bit modification of Jon's codes 
% Test in just expert two angles
% 
% 1. compare lambda values in predicting angle and choice
% 2. Change lambda values and look at the prediction performance and ratio between coefficients of dKv and dPhi
% 3. Run ridge regression and compare performance, lambda values, and the coefficient ratios.

%% basic settings
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
decisionPoint = 'lick'; % 'answer' or 'lick', added 2019/07/29 JK, for using touches before the first lick (within the answer period) or the answer lick
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

%% re-run using the same dataset (from saved files)
for targetType = 1 : 2 % touch or choice
    if targetType == 1
        if contains(sessionGroupName, 'Discrete')
            yOut = 'Discrete';
        else
            yOut = 'Ttype';
        end
        yOutForName = 'Touch';
    else
        yOut = 'Choice';
        yOutForName = 'Choice';
    end
    
    loadfn = ['mdl', sessionGroupName, Xhow, yOutForName, '_new_', decisionPoint];
    savefn = ['CV', sessionGroupName, Xhow, yOutForName, '_new_', decisionPoint];
    load([baseDir, loadfn])
    tempMdl = groupMdl;
    
    for gi = 1 : length(tempMdl)
        mdl = tempMdl{gi};
        DmatX = mdl.io.X;
        DmatY = mdl.io.Y; 
        if numel(unique(DmatY))==2 % BINOMIAL GLM MODEL 
            mdl = binomialModel_test(mdl,DmatX,DmatY,glmnetOpt);
            mdl.logDist = 'binomial';
            groupMdl{gi} = mdl;

        else % MULTINOMIAL GLM MODEL 
            mdl = multinomialModel_test(mdl,DmatX,DmatY,glmnetOpt);
            mdl.logDist = 'multinomial';
            groupMdl{gi} = mdl;
        end
    end
    save(savefn, 'info', 'groupMdl')
end


%% Example traces of how the lambda are selected.

lambda = zeros(6,2); % (:,1) for angle, (:,2) for choice
for targetType = 1 : 2 % touch or choice
    if targetType == 1
        if contains(sessionGroupName, 'Discrete')
            yOut = 'Discrete';
        else
            yOut = 'Ttype';
        end
        yOutForName = 'Touch';
    else
        yOut = 'Choice';
        yOutForName = 'Choice';
    end
    
    loadfn = ['CV', sessionGroupName, Xhow, yOutForName, '_new_', decisionPoint];
    load([baseDir, loadfn])
    for i = 1 : 6
        lambda(i,targetType) = mean([groupMdl{i}.cv.lambda_1se]);
    end
end

%%
cvglmnetPlot(groupMdl{1}.cv(1))
%%
fit = groupMdl{1}.cv(1).glmnet_fit;
figure
plot(log(groupMdl{1}.cv(1).lambda), groupMdl{1}.cv(1).cvm, 'ro-'), hold on
errorbar(log(groupMdl{1}.cv(1).lambda), groupMdl{1}.cv(1).cvm, groupMdl{1}.cv(1).cvsd/sqrt(5), 'color', [0.7 0.7 0.7], 'linestyle', 'none')
ylimVal = ylim();
plot(ones(1,2)*log(groupMdl{1}.cv(1).lambda_min), ylimVal, 'b--', 'linewidth', 1)
plot(ones(1,2)*log(groupMdl{1}.cv(1).lambda_1se), ylimVal, 'b--', 'linewidth', 2)
%%
ind__1se = find(groupMdl{1}.cv(1).lambda == groupMdl{1}.cv(1).lambda_1se);
groupMdl{1}.cv(1).cvm(ind__1se) 
indMin = find(groupMdl{1}.cv(1).lambda == groupMdl{1}.cv(1).lambda_min);
groupMdl{1}.cv(1).cvm(indMin) + groupMdl{1}.cv(1).cvsd(indMin)
groupMdl{1}.cv(1).cvm(ind__1se-1) 



%% Example trace of how coefficients change (dPhi and dKappaV)

load('CVTwoExpertMeanTouch_new_lick')
figure
for mi = 1 : 6
    cv = groupMdl{mi}.cv(1);
    fit = cv.glmnet_fit;
    subplot(2,3,mi),
    plot(log(fit.lambda), fit.beta(2,:)), hold on
    plot(log(fit.lambda), fit.beta(4,:))
    ylimVal = ylim();
    plot(ones(1,2)*log(cv.lambda_1se), ylimVal, 'b--')
    if mi == 3
        legend({'\Delta\phi', '\Delta\kappa_V'}, 'location', 'southeast')
    end
    if mi == 2
        title('Angle')
    end
    if mi == 4
        xlabel('log(\lambda)')
        ylabel('Coefficients')
    end
end
load('CVTwoExpertMeanChoice_new_lick')
figure
for mi = 1 : 6
    cv = groupMdl{mi}.cv(1);
    fit = cv.glmnet_fit;
    subplot(2,3,mi),
    plot(log(fit.lambda), fit.beta(2,:)), hold on
    plot(log(fit.lambda), fit.beta(4,:))
    ylimVal = ylim();
    plot(ones(1,2)*log(cv.lambda_1se), ylimVal, 'b--')
    if mi == 3
        legend({'\Delta\phi', '\Delta\kappa_V'}, 'location', 'southeast')
    end
    if mi == 2
        title('Choice')
    end
    if mi == 4
        xlabel('log(\lambda)')
        ylabel('Coefficients')
    end
end



%% Compare the values of lambda
% mean from 10 iterations
lambda = zeros(6,2); % (:,1) for angle, (:,2) for choice
for targetType = 1 : 2 % touch or choice
    if targetType == 1
        if contains(sessionGroupName, 'Discrete')
            yOut = 'Discrete';
        else
            yOut = 'Ttype';
        end
        yOutForName = 'Touch';
    else
        yOut = 'Choice';
        yOutForName = 'Choice';
    end
    
    loadfn = ['CV', sessionGroupName, Xhow, yOutForName, '_new_', decisionPoint];
    load([baseDir, loadfn])
    for i = 1 : 6
        lambda(i,targetType) = mean([groupMdl{i}.cv.lambda_1se]);
    end
end

%%
figure,
hold on
bar(1:2, mean(log(lambda)), 'k')
errorbar(1:2, mean(log(lambda)), std(log(lambda))/sqrt(6), 'k', 'linestyle', 'none')
xticks([1 2]), xticklabels({'Angle', 'Choice'}), ylabel('log(\lambda)')

%% Results: higher lambda in choice prediction.
%% Does it affect coeff ratio between dKv (4) and dPhi (2)? 

%% first, from each mouse
data = struct;
for targetType = 1 : 2 % touch or choice
    if targetType == 1
        if contains(sessionGroupName, 'Discrete')
            yOut = 'Discrete';
        else
            yOut = 'Ttype';
        end
        yOutForName = 'Touch';
    else
        yOut = 'Choice';
        yOutForName = 'Choice';
    end
    
    loadfn = ['CV', sessionGroupName, Xhow, yOutForName, '_new_', decisionPoint];
    temp = load([baseDir, loadfn]);
    data(targetType).data = temp;
end

%%
angleXlim = zeros(6,2); % min max
choiceXlim = zeros(6,2); % min max
lambda = zeros(6,2); % angle choice
coeffRatioInAngle = zeros(6,2); % 1: from original lambda, 2: from the other side of lambda
for targetType = 1 : 2
    for i = 1 : 6
        cv = data(targetType).data.groupMdl{i}.cv;
        lambdaMean = mean([cv.lambda_1se]);
        lambdaStd = std([cv.lambda_1se]);
        angleXlim(i,:) = [lambdaMean - 2 * lambdaStd, lambdaMean + 2 * lambdaStd];
        lambda(i,targetType) = lambdaMean;
    end
end
for targetType = 1 : 2
    figure, hold on
    for i = 1 : 6
        subplot(2,3,i), hold on
        cv = data(targetType).data.groupMdl{i}.cv;
        for j = 1 : 10
            ratio = abs(cv(j).glmnet_fit.beta(2,:)) ./ abs(cv(j).glmnet_fit.beta(4,:));
            h = plot(cv(j).lambda, ratio);
            plot(cv(j).lambda_1se, ratio(find(cv(j).lambda == cv(j).lambda_1se)), '.', 'markersize', 20, 'color', get(h,'color'))
            [~,closestInd] = min(abs(cv(j).lambda - lambda(i,3-targetType)));
            plot(cv(j).lambda(closestInd), ratio(closestInd), 'o', 'markersize', 5, 'color', get(h,'color'))
        end
        xlim([min(angleXlim(i,1), choiceXlim(i,1)), max(angleXlim(i,2), choiceXlim(i,2))])
        if i == 2
            if targetType == 1
                title('Angle')
            else
                title('Choice')
            end
        end
        if i == 4
            ylabel('coeff(\Delta\phi) / coeff(\Delta\kappa_V) (abs)')
        end
    end
end

%% Then, compare the averaged abs coefficients

angleXlim = zeros(6,2); % min max
choiceXlim = zeros(6,2); % min max
lambda = zeros(6,2); % angle choice
for targetType = 1 : 2
    for i = 1 : 6
        cv = data(targetType).data.groupMdl{i}.cv;
        lambdaMean = mean([cv.lambda_1se]);
        lambdaStd = std([cv.lambda_1se]);
        angleXlim(i,:) = [lambdaMean - 2 * lambdaStd, lambdaMean + 2 * lambdaStd];
        lambda(i,targetType) = lambdaMean;
    end
end
absCoeff = zeros(6,4,2); % mice, [dkv from og, dphi from os, dkv from the other side, dphi from tos], [angle, choice]
% absCoeffAnglefromAnglelambda = zeros(6,2); % 1: dkv, 2: dphi 
% absCoeffAnglefromChoicelambda = zeros(6,2); % 1: dkv, 2: dphi 
% absCoeffChoicefromChoicelambda = zeros(6,2); % 1: dkv, 2: dphi
% absCoeffChoicefromAnglelambda = zeros(6,2); % 1: dkv, 2: dphi
for targetType = 1 : 2
    for i = 1 : 6
        subplot(2,3,i), hold on
        cv = data(targetType).data.groupMdl{i}.cv;
        temp = zeros(10,4);        
        for j = 1 : 10
            lambdaInd = find(cv(j).lambda == cv(j).lambda_1se);
            temp(j,1) = abs(cv(j).glmnet_fit.beta(4,lambdaInd)); % dKv
            temp(j,2) = abs(cv(j).glmnet_fit.beta(2,lambdaInd)); % dPhi
            
            [~,closestInd] = min(abs(cv(j).lambda - lambda(i,3-targetType)));
            temp(j,3) = abs(cv(j).glmnet_fit.beta(4,closestInd)); % dKv
            temp(j,4) = abs(cv(j).glmnet_fit.beta(2,closestInd)); % dPhi
        end
        absCoeff(i,:,targetType) = mean(temp);
    end
end

%%
figure,
i = 1;
    subplot(1,2,i)
    bar([1,2,5,6], mean(squeeze(absCoeff(:,:,i))), 'k'), hold on
    errorbar([1,2,5,6], mean(squeeze(absCoeff(:,:,i))), std(squeeze(absCoeff(:,:,i)))/sqrt(6), 'k', 'linestyle', 'none')
    xticks([1,2,5,6])
    xticklabels({'\Delta\kappa_V lambda(Angle)', '\Delta\phi lambda(Angle)', '\Delta\kappa_V lambda(Choice)', '\Delta\phi lambda(Choice)'})
    xtickangle(45)
    title('Angle')
i = 2;
    subplot(1,2,i)
    bar([1,2,5,6], mean(squeeze(absCoeff(:,:,i))), 'k'), hold on
    errorbar([1,2,5,6], mean(squeeze(absCoeff(:,:,i))), std(squeeze(absCoeff(:,:,i)))/sqrt(6), 'k', 'linestyle', 'none')
    xticks([1,2,5,6])
    xticklabels({'\Delta\kappa_V lambda(Choice)', '\Delta\phi lambda(Choice)', '\Delta\kappa_V lambda(Angle)', '\Delta\phi lambda(Angle)'})
    xtickangle(45)
    title('Choice')
p = zeros(1,4);
[~,p(1)] = ttest(absCoeff(:,1,1), absCoeff(:,2,1));
[~,p(2)] = ttest(absCoeff(:,3,1), absCoeff(:,4,1));
[~,p(3)] = ttest(absCoeff(:,1,2), absCoeff(:,2,2));
[~,p(4)] = ttest(absCoeff(:,3,2), absCoeff(:,4,2));

%% Results: Lambda affects coefficients and their difference a little bit, but not enough to change the differential effect on coefficients of ??V and ?? between the angle model and the choice model.




%% Is correct rate the sole reason in difference between ??V and ?? between the angle model and the choice model?
% Generate random wrong trials, as same as the actual performance, and see if their difference are similar.
% Increase wrong trials and see how it shapes the relationship.

%% method 1. 
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

simNum = 100;

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

tempMdl = touchData.groupMdl;

betaSim = zeros(numMice,simNum,2); % 1 for dKv, 2 for dPhi
modelAccuracy = zeros(numMice, simNum);
performance = zeros(numMice,1);
lambdaSim = zeros(numMice,simNum);
for mi = 1 : numMice
    mdl = touchData.groupMdl{mi};
    DmatX = mdl.io.X;
    DmatY = mdl.io.Y; 
    
    choiceY = choiceData.groupMdl{mi}.io.Y;    
    numDiff = length(find(abs(choiceY-DmatY)));
    performance(mi) = 1-mean(abs(choiceY-DmatY));

    for si = 1 : simNum
        fprintf('sim %03d/%03d in mouse #%d/%d\n',si, simNum, mi, numMice)
        tempY = DmatY;
        ind2flip = randperm(length(DmatY), numDiff);
        tempY(ind2flip) = 1-tempY(ind2flip);

        randMdl = binomialModel_test(mdl,DmatX,tempY,glmnetOpt);
        betaSim(mi,si,1) = mean(randMdl.fitCoeffs(5,:)); % it's 5 because beta_0 is already included
        betaSim(mi,si,2) = mean(randMdl.fitCoeffs(3,:));
        modelAccuracy(mi,si) = mean(randMdl.gof.modelAccuracy);
        lambdaSim(mi,si) = mean([randMdl.cv.lambda_1se]);
    end
end

%%
saveDir = 'C:\Users\shires\Dropbox\Works\Data analysis\Object angle coding and learning\';
saveFn = 'd190917_data_1';
save([saveDir, saveFn], 'betaSim', 'modelAccuracy', 'lambdaSim')


%%
figure, hold on
plot(cellfun(@(x) mean(x.gof.modelAccuracy),choiceData.groupMdl), 'ro')
plot(cellfun(@(x) mean(x.gof.modelAccuracy),touchData.groupMdl), 'bo')

bar(mean(modelAccuracy,2), 'k')
errorbar(mean(modelAccuracy,2), std(modelAccuracy,[],2)/sqrt(6), 'k', 'linestyle', 'none')
plot(cellfun(@(x) mean(x.gof.modelAccuracy),choiceData.groupMdl), 'ro')
plot(cellfun(@(x) mean(x.gof.modelAccuracy),touchData.groupMdl), 'bo')

legend({'Choice model', 'Angle model'})
xlabel('Mouse #'), ylabel('Model accuracy')
ylim([0.5 1])

%% 
figure, hold on
mi = 1;
plot([mi-0.2, mi+0.2], mean(abs(choiceData.groupMdl{mi}.fitCoeffs([5,3],:)),2), 'ro')
plot([mi-0.2, mi+0.2], mean(abs(touchData.groupMdl{mi}.fitCoeffs([5,3],:)),2), 'bo')

for mi = 1 : numMice
    
    mat = squeeze(abs(betaSim(mi,:,:)));
    bar([mi-0.2, mi+0.2], mean(mat), 0.8, 'k')
    errorbar([mi-0.2, mi+0.2], mean(mat), std(mat), 'k', 'linestyle', 'none')
    
    plot([mi-0.2, mi+0.2], mean(abs(choiceData.groupMdl{mi}.fitCoeffs([5,3],:)),2), 'ro-')
    plot([mi-0.2, mi+0.2], mean(abs(touchData.groupMdl{mi}.fitCoeffs([5,3],:)),2), 'bo-')
    
    
end
    
legend({'Choice model', 'Angle model'})
xlabel('Mouse #'), ylabel('Coefficients')
    

%%

figure, 
subplot(1,3,1), hold on
mat = zeros(mi,2);
for mi = 1 : 6
    mat(mi,:) = mean(abs(touchData.groupMdl{mi}.fitCoeffs([5,3],:)),2)';
end
bar(mean(mat), 'k')
errorbar(mean(mat), std(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
title('Angle model')
[~,p] = ttest(mat(:,1), mat(:,2))

subplot(1,3,2), hold on
mat = zeros(mi,2);
for mi = 1 : 6
    mat(mi,:) = mean(abs(choiceData.groupMdl{mi}.fitCoeffs([5,3],:)),2)';
end
bar(mean(mat), 'k')
errorbar(mean(mat), std(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
title('Choice model')
[~,p] = ttest(mat(:,1), mat(:,2))

subplot(1,3,3), hold on
mat = zeros(mi,2);
for mi = 1 : 6
    mat(mi,:) = mean(squeeze(abs(betaSim(mi,:,:))));
end
bar(mean(mat), 'k')
errorbar(mean(mat), std(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks([1 2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
title('Simulated choice')
[~,p] = ttest(mat(:,1), mat(:,2))




%% Results: The difference is NOT shown by randomly distributin wrong trials. 



%% Simulatin #2
% But, do I really need this?


baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

crList = zeros(11,1);
for i = 1 : 11
    crList(i) = (100-5*(i-1))/100;
end
simNum = 100;

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95;
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

tempMdl = touchData.groupMdl;

betaSim = zeros(numMice, length(crList),2); % 1 for dKv, 2 for dPhi
modelAccuracy = zeros(numMice, length(crList));
lambdaSim = zeros(numMice, length(crList)); 

for mi = 1 : numMice
    mdl = touchData.groupMdl{mi};
    DmatX = mdl.io.X;
    DmatY = mdl.io.Y; 
    
    for ci = 1 : length(crList)
        fprintf('sim repeat %03d/%03d in mouse #%d/%d\n',ci, length(crList), mi, numMice)
        numDiff = round(length(DmatY)*(1-crList(ci)));
        
        beta = zeros(simNum,2);
        ma = zeros(simNum,1);
        lambda = zeros(simNum,1);
        for si = 1 : simNum
            
            tempY = DmatY;
            ind2flip = randperm(length(DmatY), numDiff);
            tempY(ind2flip) = 1-tempY(ind2flip);

            randMdl = binomialModel_test(mdl,DmatX,tempY,glmnetOpt);
            beta(si,1) = mean(randMdl.fitCoeffs(5,:)); % it's 5 because beta_0 is already included
            beta(si,2) = mean(randMdl.fitCoeffs(3,:));
            ma(si) = mean(randMdl.gof.modelAccuracy);
            lambda(si) = mean([randMdl.cv.lambda_1se]);
        end
        
        betaSim(mi,ci,1) = mean(beta(:,1)); 
        betaSim(mi,ci,2) = mean(beta(:,2));
        modelAccuracy(mi,ci) = mean(ma);
        lambdaSim(mi,ci) = mean(lambda);

    end
end

%%
saveDir = 'C:\Users\shires\Dropbox\Works\Data analysis\Object angle coding and learning\';
saveFn = 'd190917_data_2';
save([saveDir, saveFn], 'betaSim', 'modelAccuracy', 'lambdaSim')

%%
colorList = zeros(11,3);
for i = 1 : 11
    colorList(i,:) = ones(1,3)*(i-1)*0.08;
end

p = zeros(13,1);
figure, 
subplot(121), hold on
mat = cell2mat(cellfun(@(x) abs(mean(x.fitCoeffs([5,3],:),2)), touchData.groupMdl', 'uniformoutput', false))';
errorbar([1,2], mean(mat), std(mat)/sqrt(6), 'bo-', 'linewidth', 1)
[~,p(1)] = ttest(mat(:,1), mat(:,2));
mat = cell2mat(cellfun(@(x) abs(mean(x.fitCoeffs([5,3],:),2)), choiceData.groupMdl', 'uniformoutput', false))';
errorbar([1,2], mean(mat), std(mat)/sqrt(6), 'ro-', 'linewidth', 1)
[~,p(2)] = ttest(mat(:,1), mat(:,2));
for i = 1 : 11
    mat = abs(squeeze(betaSim(:,i,:)));
    errorbar([1,2], mean(mat), std(mat)/sqrt(6), '-o', 'color', colorList(i,:), 'linewidth', 1)
    [~,p(i+2)] = ttest(mat(:,1), mat(:,2));
end

xlim([0 3])
xticks([1,2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
ylabel('Absolute coefficient')

subplot(122), hold on
mat = cell2mat(cellfun(@(x) abs(mean(x.fitCoeffs([5,3],:),2)), touchData.groupMdl', 'uniformoutput', false))';
ratio = mat(:,2)./mat(:,1);
errorbar(1, mean(ratio), std(ratio)/sqrt(6), 'bo', 'linestyle', 'none')

mat = cell2mat(cellfun(@(x) abs(mean(x.fitCoeffs([5,3],:),2)), choiceData.groupMdl', 'uniformoutput', false))';
ratio = mat(:,2)./mat(:,1);
errorbar(3.2, mean(ratio), std(ratio)/sqrt(6), 'ro', 'linestyle', 'none')

for i = 1 : 11
    mat = abs(squeeze(betaSim(:,i,:)));
    ratio = mat(:,2)./mat(:,1);
    errorbar(i, mean(ratio), std(ratio)/sqrt(6), 'o', 'color', colorList(i,:), 'linestyle', 'none')
end

xticks([1:11])
xtlCell = cell(1,11);
for i = 1 : 11
    xtlCell{i} = sprintf('%d',100-(i-1)*5);
end
xticklabels(xtlCell)
xlabel('Correct rate (%)')
ylabel('Absolute coefficient ratio (\Delta\phi / \Delta\kappa_V)')

legendList = cell(1,13);
legendList{1} = sprintf('Angle p = %.2f', p(1));
legendList{2} = sprintf('Choice p = %.2f', p(2));
for i = 3 : 13
    legendList{i} = sprintf('%d %% p = %.2f',100-(i-3)*5, p(i));
end

legend(legendList)

yyaxis right
errorbar(1:11, mean(lambdaSim), std(lambdaSim)/sqrt(6), 'g.')

%% Results: again showing specially low ratio in real choice data

%% What is so special about these wrong trials?
% look at the feature distribution
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'



histNormMeth = 'cdf'; % or 'probability'

crList = zeros(11,1);
for i = 1 : 11
    crList(i) = (100-5*(i-1))/100;
end
simNum = 100;

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

histRange = -3:0.2:3.1;
distDkv = zeros(6,length(histRange)-1,2,3); % 3rd dim, 1: 135 or left, 2: 45 or right.   4th dim, 1: angle, 2: correct, 3: wrong
distDp = zeros(6,length(histRange)-1,2,3); % 3rd dim, 1: 135 or left, 2: 45 or right.   4th dim, 1: angle, 2: correct, 3: wrong
for mi = 1 : 6
    angleY = touchData.groupMdl{mi}.io.Y;
    choiceY = choiceData.groupMdl{mi}.io.Y;
    
    correctTrialInds = find((angleY-choiceY)==0);
    wrongTrialInds = find(abs(angleY-choiceY)==1);
    
    for i = 1 : 2
        Inds = find(angleY == i-1);
        cInds = intersect(correctTrialInds, Inds);
        wInds = intersect(wrongTrialInds, Inds);
        
        distDkv(mi,:,i,1) = histcounts(touchData.groupMdl{mi}.io.X(Inds,4), histRange, 'normalization', histNormMeth) * length(Inds) / length(angleY);
        distDp(mi,:,i,1) = histcounts(touchData.groupMdl{mi}.io.X(Inds,2), histRange, 'normalization', histNormMeth) * length(Inds) / length(angleY);

        distDkv(mi,:,i,2) = histcounts(touchData.groupMdl{mi}.io.X(cInds,4), histRange, 'normalization', histNormMeth) * length(cInds) / length(correctTrialInds);
        distDp(mi,:,i,2) = histcounts(touchData.groupMdl{mi}.io.X(cInds,2), histRange, 'normalization', histNormMeth) * length(cInds) / length(correctTrialInds);

        distDkv(mi,:,i,3) = histcounts(touchData.groupMdl{mi}.io.X(wInds,4), histRange, 'normalization', histNormMeth) * length(wInds) / length(wrongTrialInds);
        distDp(mi,:,i,3) = histcounts(touchData.groupMdl{mi}.io.X(wInds,2), histRange, 'normalization', histNormMeth) * length(wInds) / length(wrongTrialInds);
    end
end
%%
set(0,'defaultaxesfontname', 'Arial')
set(0,'defaultaxesfontsize', 12)
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.5]) 
subplot(121), hold on % for correct trials
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), std(squeeze(distDkv(:,:,2,1)))/sqrt(6), 'lineprop', 'b')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), std(squeeze(distDkv(:,:,1,1)))/sqrt(6), 'lineprop', 'r')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,2))), std(squeeze(distDkv(:,:,2,2)))/sqrt(6), 'lineprop', 'c')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,2))), std(squeeze(distDkv(:,:,1,2)))/sqrt(6), 'lineprop', 'm')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), 'b')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), 'r')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,2))), 'c')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,2))), 'm')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), std(squeeze(distDkv(:,:,2,1)))/sqrt(6), 'b')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), std(squeeze(distDkv(:,:,1,1)))/sqrt(6), 'r')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,2))), std(squeeze(distDkv(:,:,2,2)))/sqrt(6), 'c')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,2))), std(squeeze(distDkv(:,:,1,2)))/sqrt(6), 'm')
if strcmp(histNormMeth,'cdf')
    ylim([0 0.8])
    ylabel('Cumulative proportion')
else
    ylim([0 0.15])
    ylabel('Proportion')
end
xlabel('\Delta\kappa_V (standardized)')
legend({'45\circ', '135\circ', '45\circ correct', '135\circ correct'}, 'box', 'off')

subplot(122), hold on % for wrong trials
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), std(squeeze(distDkv(:,:,2,1)))/sqrt(6), 'lineprop', 'b')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), std(squeeze(distDkv(:,:,1,1)))/sqrt(6), 'lineprop', 'r')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,3))), std(squeeze(distDkv(:,:,2,3)))/sqrt(6), 'lineprop', 'm')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,3))), std(squeeze(distDkv(:,:,1,3)))/sqrt(6), 'lineprop', 'c')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), 'b')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), 'r')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,3))), 'm')
plot(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,3))), 'c')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,1))), std(squeeze(distDkv(:,:,2,1)))/sqrt(6), 'b')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,1))), std(squeeze(distDkv(:,:,1,1)))/sqrt(6), 'r')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,2,3))), std(squeeze(distDkv(:,:,2,3)))/sqrt(6), 'm')
boundedline(histRange(1:end-1), mean(squeeze(distDkv(:,:,1,3))), std(squeeze(distDkv(:,:,1,3)))/sqrt(6), 'c')
if strcmp(histNormMeth,'cdf')
    ylim([0 0.8])
    ylabel('Cumulative proportion')
else
    ylim([0 0.15])
    ylabel('Proportion')
end
xlabel('\Delta\kappa_V (standardized)')

legend({'45\circ', '135\circ', '45\circ wrong', '135\circ wrong'}, 'box', 'off')

figure('units', 'normalized', 'outerposition', [0.1 0.1 0.7 0.5]) 
subplot(121), hold on % for correct trials
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), std(squeeze(distDp(:,:,2,1)))/sqrt(6), 'lineprop', 'b')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), std(squeeze(distDp(:,:,1,1)))/sqrt(6), 'lineprop', 'r')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,2,2))), std(squeeze(distDp(:,:,2,2)))/sqrt(6), 'lineprop', 'c')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,1,2))), std(squeeze(distDp(:,:,1,2)))/sqrt(6), 'lineprop', 'm')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), 'b')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), 'r')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,2,2))), 'c')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,1,2))), 'm')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), std(squeeze(distDp(:,:,2,1)))/sqrt(6), 'b')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), std(squeeze(distDp(:,:,1,1)))/sqrt(6), 'r')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,2,2))), std(squeeze(distDp(:,:,2,2)))/sqrt(6), 'c')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,1,2))), std(squeeze(distDp(:,:,1,2)))/sqrt(6), 'm')
if strcmp(histNormMeth,'cdf')
    ylim([0 0.8])
    ylabel('Cumulative proportion')
else
    ylim([0 0.15])
    ylabel('Proportion')
end
xlabel('\Delta\phi (standardized)')

legend({'45\circ', '135\circ', '45\circ correct', '135\circ correct'}, 'box', 'off')

subplot(122), hold on % for wrong trials
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), std(squeeze(distDp(:,:,2,1)))/sqrt(6), 'lineprop', 'b')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), std(squeeze(distDp(:,:,1,1)))/sqrt(6), 'lineprop', 'r')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,2,3))), std(squeeze(distDp(:,:,2,3)))/sqrt(6), 'lineprop', 'm')
% shadedErrorBar(histRange(1:end-1), mean(squeeze(distDp(:,:,1,3))), std(squeeze(distDp(:,:,1,3)))/sqrt(6), 'lineprop', 'c')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), 'b')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), 'r')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,2,3))), 'm')
plot(histRange(1:end-1), mean(squeeze(distDp(:,:,1,3))), 'c')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,2,1))), std(squeeze(distDp(:,:,2,1)))/sqrt(6), 'b')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,1,1))), std(squeeze(distDp(:,:,1,1)))/sqrt(6), 'r')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,2,3))), std(squeeze(distDp(:,:,2,3)))/sqrt(6), 'm')
boundedline(histRange(1:end-1), mean(squeeze(distDp(:,:,1,3))), std(squeeze(distDp(:,:,1,3)))/sqrt(6), 'c')
if strcmp(histNormMeth,'cdf')
    ylim([0 0.8])
    ylabel('Cumulative proportion')
else
    ylim([0 0.15])
    ylabel('Proportion')
end
xlabel('\Delta\phi (standardized)')
legend({'45\circ', '135\circ', '45\circ wrong', '135\circ wrong'}, 'box', 'off')




%%

saveDir = 'C:\Users\shires\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig2-S5 Wrong trials dkv and dphi\';
fn = 'dkv_cumprop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
fix_eps_fonts([saveDir, fn])


%%

saveDir = 'C:\Users\shires\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig2-S5 Wrong trials dkv and dphi\';
fn = 'dphi_cumprop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Results: wrong trials in 135 degrees have dKv distribution closer to 45 degrees, which should have made it easier to be fooled.
% Wrong trials in 45 degrees have dPhi distribution farther from 135 degrees, which should have made it harder to be fooled.
% This result corroborates with their variable importance and coefficients,
% explaining how very similar distribution between choice and angle led to
% different results in variable importance and coefficients.

%% Last check - what about in case of ridge regression?

baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

simNum = 100;

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0; % for ridge
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

tempMdl = touchData.groupMdl;

betaDkv = zeros(numMice,2); % 1 for angle model, 2 for choice model
betaDp = zeros(numMice,2); % 1 for angle model, 2 for choice model
modelAccuracy = zeros(numMice, 2); % 1 for angle model, 2 for choice model
lambda = zeros(numMice,2); % 1 for angle model, 2 for choice model

for mi = 1 : numMice
    mdl = touchData.groupMdl{mi};
    DmatX = mdl.io.X;
    DmatY = mdl.io.Y; 
    
    mdl = binomialModel_test(mdl, DmatX, DmatY, glmnetOpt);
    
    betaDkv(mi,1) = mean(mdl.fitCoeffs(5,:));
    betaDp(mi,1) = mean(mdl.fitCoeffs(3,:));
    modelAccuracy(mi,1) = mean(mdl.gof.modelAccuracy);
    lambda(mi,1) = mean([mdl.cv.lambda_1se]);
    
    mdl = choiceData.groupMdl{mi};
    DmatX = mdl.io.X;
    DmatY = mdl.io.Y; 
    
    mdl = binomialModel_test(mdl, DmatX, DmatY, glmnetOpt);
    
    betaDkv(mi,2) = mean(mdl.fitCoeffs(5,:));
    betaDp(mi,2) = mean(mdl.fitCoeffs(3,:));
    modelAccuracy(mi,2) = mean(mdl.gof.modelAccuracy);
    lambda(mi,2) = mean([mdl.cv.lambda_1se]);
end

%%
figure,
subplot(221), hold on
bar(mean(modelAccuracy), 'k')
errorbar(mean(modelAccuracy), std(modelAccuracy)/sqrt(6), 'k', 'linestyle', 'none')
xlim([0 3])
xticks([1 2])
xticklabels({'Angle', 'Choice'})
xtickangle(45)
title('Model accuracy')

subplot(222), hold on
bar(mean(lambda), 'k')
errorbar(mean(lambda), std(lambda)/sqrt(6), 'k', 'linestyle', 'none')
xlim([0 3])
xticks([1 2])
xticklabels({'Angle', 'Choice'})
xtickangle(45)
title('\lambda')

subplot(223), hold on
bar(mean([abs(betaDkv(:,1)), abs(betaDp(:,1))]), 'k')
errorbar(mean([abs(betaDkv(:,1)), abs(betaDp(:,1))]), std([abs(betaDkv(:,1)), abs(betaDp(:,1))])/sqrt(6), 'k', 'linestyle', 'none')
xlim([0 3])
xticks([1 2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
xtickangle(45)
ylim([0 1.5])
ylabel('Absolute coefficient')
title('Angle model')

subplot(224), hold on
bar(mean([abs(betaDkv(:,2)), abs(betaDp(:,2))]), 'k')
errorbar(mean([abs(betaDkv(:,2)), abs(betaDp(:,2))]), std([abs(betaDkv(:,2)), abs(betaDp(:,2))])/sqrt(6), 'k', 'linestyle', 'none')
xlim([0 3])
xticks([1 2])
xticklabels({'\Delta\kappa_V', '\Delta\phi'})
xtickangle(45)
ylim([0 1.5])
ylabel('Absolute coefficient')
title('Choice model')
%%
[~,p] = ttest(betaDkv(:,1), betaDp(:,1))
[~,p] = ttest(betaDkv(:,2), betaDp(:,2))
title('Angle model')

%% Results: The difference is increased in choice model, but it is not significant


%% what about in ratio?
figure,
mat = betaDp./betaDkv;
bar(mean(mat), 'k'), hold on
errorbar(mean(mat), std(mat)/sqrt(6), 'k', 'linestyle', 'none')
xlim([0 3])
xticks([1 2])
xticklabels({'Angle', 'Choice'})
xtickangle(45)
ylabel('Ratio abs(\Delta\phi)/abs(\Delta\kappa_V)')


%% What if I train the model only in wrong trials?

baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0; % for ridge
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

tempMdl = touchData.groupMdl;

betaDkv = zeros(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
betaDp = zeros(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
modelAccuracy = zeros(numMice, 3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
lambda = zeros(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model

for mi = 1 : numMice
    mdl = touchData.groupMdl{mi};
    DmatXAngle = mdl.io.X;
    DmatYAngle = mdl.io.Y; 
    
    mdl = binomialModel_test(mdl, DmatXAngle, DmatYAngle, glmnetOpt); % I have to re-train for lambda
    betaDkv(mi,1) = mean(mdl.fitCoeffs(5,:));
    betaDp(mi,1) = mean(mdl.fitCoeffs(3,:));
    modelAccuracy(mi,1) = mean(mdl.gof.modelAccuracy);
    lambda(mi,1) = mean([mdl.cv.lambda_1se]);
    
    mdl = choiceData.groupMdl{mi};
    DmatXChoice = mdl.io.X;
    DmatYChoice = mdl.io.Y; 
    
    wrongInds = find(abs(DmatYChoice - DmatYAngle) == 1);
    correctInds = find(abs(DmatYChoice - DmatYAngle) == 0);
    
    tempBetaDkv = nan(10,1);
    tempBetaDp = nan(10,1);
    tempModelAccuracy = nan(10,1);
    tempLambda = nan(10,1);
    for ri = 1 : 10
        correctInds2train = correctInds(randperm(length(correctInds), length(wrongInds)));
        try
            mdl = binomialModel_test(mdl, DmatXChoice(correctInds2train,:), DmatYChoice(correctInds2train), glmnetOpt);
            tempBetaDkv(ri) = mean(mdl.fitCoeffs(5,:));
            tempBetaDp(ri) = mean(mdl.fitCoeffs(3,:));
            tempModelAccuracy(ri) = mean(mdl.gof.modelAccuracy);
            tempLambda(ri) = mean([mdl.cv.lambda_1se]);
        catch
        end
    end
    
    betaDkv(mi,2) = nanmean(tempBetaDkv);
    betaDp(mi,2) = nanmean(tempBetaDp);
    modelAccuracy(mi,2) = nanmean(tempModelAccuracy);
    lambda(mi,2) = nanmean(tempLambda);

    
    mdl = binomialModel_test(mdl, DmatXChoice(wrongInds,:), DmatYChoice(wrongInds), glmnetOpt);
    
    betaDkv(mi,3) = mean(mdl.fitCoeffs(5,:));
    betaDp(mi,3) = mean(mdl.fitCoeffs(3,:));
    modelAccuracy(mi,3) = mean(mdl.gof.modelAccuracy);
    lambda(mi,3) = mean([mdl.cv.lambda_1se]);
end

%%
xval = [1 4 7];
figure, 
subplot(211), hold on
bar(xval-0.4, mean(abs(betaDkv)), 0.2, 'k')
errorbar(xval-0.4, mean(abs(betaDkv)), std(abs(betaDkv))/sqrt(6), 'k', 'linestyle', 'none')

bar(xval+0.4, mean(abs(betaDp)), 0.2, 'k')
errorbar(xval+0.4, mean(abs(betaDp)), std(abs(betaDp))/sqrt(6), 'k', 'linestyle', 'none')

xticks(xval)
xticklabels({'', '', ''})
ylabel('Absolute coefficient')
title('Ridge')


[~,p] = ttest(abs(betaDkv(:,1)), abs(betaDp(:,1)))
[~,p] = ttest(abs(betaDkv(:,2)), abs(betaDp(:,2)))
[~,p] = ttest(abs(betaDkv(:,3)), abs(betaDp(:,3)))

subplot(212), hold on
mat = abs(betaDp./betaDkv);
bar(xval, mean(mat), 0.2, 'k')
errorbar(xval, mean(mat), std(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks(xval)
xticklabels({'Angle model', 'Correct model', 'Wrong model'})
xtickangle(45)
ylabel('Ratio abs(\Delta\phi)/abs(\Delta\kappa_V)')

[~,p] = ttest(mat(:,1), mat(:,2))
[~,p] = ttest(mat(:,1), mat(:,3))
[~,p] = ttest(mat(:,2), mat(:,3))

%%
xval = [1 4 7];
figure, 
subplot(211), hold on
bar(xval-0.4, nanmean(abs(betaDkv)), 0.2, 'k')
errorbar(xval-0.4, nanmean(abs(betaDkv)), nanstd(abs(betaDkv))/sqrt(6), 'k', 'linestyle', 'none')

bar(xval+0.4, nanmean(abs(betaDp)), 0.2, 'k')
errorbar(xval+0.4, nanmean(abs(betaDp)), nanstd(abs(betaDp))/sqrt(6), 'k', 'linestyle', 'none')

xticks(xval)
xticklabels({'', '', ''})
ylabel('Absolute coefficient')
title('Ridge')


[~,p] = ttest(abs(betaDkv(:,1)), abs(betaDp(:,1)))
[~,p] = ttest(abs(betaDkv(:,2)), abs(betaDp(:,2)))
[~,p] = ttest(abs(betaDkv(:,3)), abs(betaDp(:,3)))

subplot(212), hold on
mat = abs(betaDp) - abs(betaDkv);
bar(xval, nanmean(mat), 0.2, 'k')
errorbar(xval, nanmean(mat), nanstd(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks(xval)
xticklabels({'Angle model', 'Correct model', 'Wrong model'})
xtickangle(45)
ylabel('Abs(\Delta\phi) - abs(\Delta\kappa_V)')

[~,p] = ttest(mat(:,1), mat(:,2))
[~,p] = ttest(mat(:,1), mat(:,3))
[~,p] = ttest(mat(:,2), mat(:,3))

%% lasso
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([baseDir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([baseDir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);

% GLM model parameters
glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; %set to 0 b/c already standardized
glmnetOpt.alpha = 0.95; % for ridge
glmnetOpt.xfoldCV = 5;
glmnetOpt.numIterations = 10;

tempMdl = touchData.groupMdl;

betaDkv = nan(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
betaDp = nan(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
modelAccuracy = nan(numMice, 3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model
lambda = nan(numMice,3); % 1 for angle model, 2 for correct choice model with the same # of wrong trials, 3 for wrong choice model

for mi = 1 : numMice
    mdl = touchData.groupMdl{mi};
    DmatXAngle = mdl.io.X;
    DmatYAngle = mdl.io.Y; 
    
    mdl = binomialModel_test(mdl, DmatXAngle, DmatYAngle, glmnetOpt); % I have to re-train for lambda
    betaDkv(mi,1) = mean(mdl.fitCoeffs(5,:));
    betaDp(mi,1) = mean(mdl.fitCoeffs(3,:));
    modelAccuracy(mi,1) = mean(mdl.gof.modelAccuracy);
    lambda(mi,1) = mean([mdl.cv.lambda_1se]);
    
    mdl = choiceData.groupMdl{mi};
    DmatXChoice = mdl.io.X;
    DmatYChoice = mdl.io.Y; 
    
    wrongInds = find(abs(DmatYChoice - DmatYAngle) == 1);
    correctInds = find(abs(DmatYChoice - DmatYAngle) == 0);
    
    tempBetaDkv = nan(10,1);
    tempBetaDp = nan(10,1);
    tempModelAccuracy = nan(10,1);
    tempLambda = nan(10,1);
    for ri = 1 : 10
        correctInds2train = correctInds(randperm(length(correctInds), length(wrongInds)));
        try
            mdl = binomialModel_test(mdl, DmatXChoice(correctInds2train,:), DmatYChoice(correctInds2train), glmnetOpt);
            tempBetaDkv(ri) = mean(mdl.fitCoeffs(5,:));
            tempBetaDp(ri) = mean(mdl.fitCoeffs(3,:));
            tempModelAccuracy(ri) = mean(mdl.gof.modelAccuracy);
            tempLambda(ri) = mean([mdl.cv.lambda_1se]);
        catch
        end
    end
    
    betaDkv(mi,2) = nanmean(tempBetaDkv);
    betaDp(mi,2) = nanmean(tempBetaDp);
    modelAccuracy(mi,2) = nanmean(tempModelAccuracy);
    lambda(mi,2) = nanmean(tempLambda);

    try
        mdl = binomialModel_test(mdl, DmatXChoice(wrongInds,:), DmatYChoice(wrongInds), glmnetOpt);

        betaDkv(mi,3) = mean(mdl.fitCoeffs(5,:));
        betaDp(mi,3) = mean(mdl.fitCoeffs(3,:));
        modelAccuracy(mi,3) = mean(mdl.gof.modelAccuracy);
        lambda(mi,3) = mean([mdl.cv.lambda_1se]);
    catch
    end
end

%%
xval = [1 4 7];
figure, 
subplot(211), hold on
bar(xval-0.4, nanmean(abs(betaDkv)), 0.2, 'k')
errorbar(xval-0.4, nanmean(abs(betaDkv)), nanstd(abs(betaDkv))/sqrt(6), 'k', 'linestyle', 'none')

bar(xval+0.4, nanmean(abs(betaDp)), 0.2, 'k')
errorbar(xval+0.4, nanmean(abs(betaDp)), nanstd(abs(betaDp))/sqrt(6), 'k', 'linestyle', 'none')

xticks(xval)
xticklabels({'', '', ''})
ylabel('Absolute coefficient')
title('Lasso')


[~,p] = ttest(abs(betaDkv(:,1)), abs(betaDp(:,1)))
[~,p] = ttest(abs(betaDkv(:,2)), abs(betaDp(:,2)))
[~,p] = ttest(abs(betaDkv(:,3)), abs(betaDp(:,3)))

subplot(212), hold on
mat = abs(betaDp./betaDkv);
bar(xval, nanmean(mat), 0.2, 'k')
errorbar(xval, nanmean(mat), nanstd(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks(xval)
xticklabels({'Angle model', 'Correct model', 'Wrong model'})
xtickangle(45)
ylabel('Ratio abs(\Delta\phi)/abs(\Delta\kappa_V)')

[~,p] = ttest(mat(:,1), mat(:,2))
[~,p] = ttest(mat(:,1), mat(:,3))
[~,p] = ttest(mat(:,2), mat(:,3))


%%
xval = [1 4 7];
figure, 
subplot(211), hold on
bar(xval-0.4, nanmean(abs(betaDkv)), 0.2, 'k')
errorbar(xval-0.4, nanmean(abs(betaDkv)), nanstd(abs(betaDkv))/sqrt(6), 'k', 'linestyle', 'none')

bar(xval+0.4, nanmean(abs(betaDp)), 0.2, 'k')
errorbar(xval+0.4, nanmean(abs(betaDp)), nanstd(abs(betaDp))/sqrt(6), 'k', 'linestyle', 'none')

xticks(xval)
xticklabels({'', '', ''})
ylabel('Absolute coefficient')
title('Lasso')


[~,p] = ttest(abs(betaDkv(:,1)), abs(betaDp(:,1)))
[~,p] = ttest(abs(betaDkv(:,2)), abs(betaDp(:,2)))
[~,p] = ttest(abs(betaDkv(:,3)), abs(betaDp(:,3)))

subplot(212), hold on
mat = abs(betaDp) - abs(betaDkv);
bar(xval, nanmean(mat), 0.2, 'k')
errorbar(xval, nanmean(mat), nanstd(mat)/sqrt(6), 'k', 'linestyle', 'none')
xticks(xval)
xticklabels({'Angle model', 'Correct model', 'Wrong model'})
xtickangle(45)
ylabel('Abs(\Delta\phi) - abs(\Delta\kappa_V)')

[~,p] = ttest(mat(:,1), mat(:,2))
[~,p] = ttest(mat(:,1), mat(:,3))
[~,p] = ttest(mat(:,2), mat(:,3))