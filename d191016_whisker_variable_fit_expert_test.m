%% What's wrong with expert mice in reconstruction from whisker model??
%% They have very high probability of being tuned even without all whisker features.
%% They are highly correlated with angle tuning responses from inferred spikes.
%% But whisker models of expert mice have high deviance explained ONLY in whisker variables (others are not even fit).
%% How come??

%% Target JK025 S19 for now.

cd('Y:\Whiskernas\JK\suite2p\025\')
data1 = load('angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK025S19');
data2 = load('angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK025S19_2ndRound');

%%
corr(cellfun(@mean, data1.spkValAllCell{1,1}), cellfun(@mean, data2.spkValAllCell{1,1}))
cellfun(@mean, data1.spkValAllCell{1,3})
%%
corr(cellfun(@mean, data1.spkValAllCell{1,1}), cellfun(@mean, data2.spkValAllCell{1,3}))

%%

cd('Y:\Whiskernas\JK\suite2p\025\')
load('angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK025S19_categoryTest.mat')

%%
corrValsNaive = zeros(size(spkValAllCell));
for ci = 1 : size(spkValAllCell,1)
    for i = 1 : size(spkValAllCell,2)
        tempCorr = corr(cellfun(@mean, data1.spkValAllCell{ci,1}), cellfun(@mean, spkValAllCell{ci,i}));
        if isnan(tempCorr) || tempCorr < 0
            corrValsNaive(ci,i) = 0;
        else
            corrValsNaive(ci,i) = tempCorr;
        end
    end
end

%%
figure, hold on
bar(mean(corrValsNaive), 'k')
errorbar(mean(corrValsNaive), sem(corrValsNaive), 'k', 'lines', 'no')
%%
xticks(1:size(spkValAllCell,1))
xticklabels(featureNames)
xtickangle(45)
%%
ylabel('Correlation')


%% from preAnswer
% clear
featureNames = {'noWhisker','noSound','noReward','noWhisking','noLicking','noWhisker&Sound','noWhisker&Reward','noWhisker&Whisking','noWhisker&Licking'};
baseDir = 'Y:\Whiskernas\JK\suite2p\';
data1 = load([baseDir, 'modelAngleTuning_NC']);
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
learnerInd = find(cellfun(@(x) length(x) == 2, sessions));
for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        if si == 1
            data.naive(mi) = load(sprintf('%s%03d\\angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK%03dS%02d_categoryTest', baseDir, mouse, mouse, session));
        else
            ei = find(learnerInd == mi);
            data.expert(ei) = load(sprintf('%s%03d\\angle_tuning_model_touchCell_NC_preAnswer_perTouch_JK%03dS%02d_categoryTest', baseDir, mouse, mouse, session));
        end
    end
end

%%
corrValNaive = cell(length(mice),1);
for mi = 1 : length(mice)
    indTuned = find(data1.naive(mi).tunedAllCell(:,1));
    indTemp = find(data1.naive(mi).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrValsNaive = zeros(length(indWhisker), size(data.naive(mi).spkValAllCell,2));
    for ci = 1 : length(indWhisker)
        for i = 1 : size(data.naive(mi).spkValAllCell,2)
            if data.naive(mi).tunedAllCell(indWhisker(ci),i)
                tempCorr = corr(cellfun(@mean, data1.naive(mi).spkValAllCell{indWhisker(ci),1}), cellfun(@mean, data.naive(mi).spkValAllCell{indWhisker(ci),i}));
                if isnan(tempCorr) || tempCorr < 0
                    corrValsNaive(ci,i) = 0;
                else
                    corrValsNaive(ci,i) = tempCorr;
                end
            else
                corrValsNaive(ci,i) = 0;
            end
        end
    end
    corrValNaive{mi} = corrValsNaive;
end

corrValExpert = cell(length(learnerInd),1);
for mi = 1 : length(learnerInd)
    indTuned = find(data1.expert(mi).tunedAllCell(:,1));
    indTemp = find(data1.expert(mi).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrValsNaive = zeros(length(indWhisker), size(data.expert(mi).spkValAllCell,2));
    for ci = 1 : length(indWhisker)
        for i = 1 : size(data.expert(mi).spkValAllCell,2)
            if data.expert(mi).tunedAllCell(indWhisker(ci),i)
                tempCorr = corr(cellfun(@mean, data1.expert(mi).spkValAllCell{indWhisker(ci),1}), cellfun(@mean, data.expert(mi).spkValAllCell{indWhisker(ci),i}));
                if isnan(tempCorr) || tempCorr < 0
                    corrValsNaive(ci,i) = 0;
                else
                    corrValsNaive(ci,i) = tempCorr;
                end
            else
                corrValsNaive(ci,i) = 0;
            end
        end
    end
    corrValExpert{mi} = corrValsNaive;
end

%%
posAdj = 0.2;
barWidth = 0.4;
xpos = 1 : length(featureNames);
naiveMat = cell2mat(cellfun(@mean, corrValNaive, 'un', 0));
expertMat = cell2mat(cellfun(@mean, corrValExpert, 'un', 0));
figure, hold on
bar(xpos - posAdj, mean(naiveMat), barWidth, 'k')
bar(xpos + posAdj, mean(expertMat), barWidth, 'c')
errorbar(xpos - posAdj, mean(naiveMat), sem(naiveMat), 'k', 'lines', 'no')
errorbar(xpos + posAdj, mean(expertMat), zeros(1,length(xpos)), sem(expertMat), 'k', 'lines', 'no')
legend({'Naive', 'Expert'})
xticks(xpos)
xticklabels(featureNames)
xtickangle(45)
ylabel('Correlation')


%%
posAdj = 0.2;
barWidth = 0.4;
xpos = 1 : length(featureNames);
naiveMat = cell2mat(cellfun(@mean, corrValNaive(learnerInd), 'un', 0));
expertMat = cell2mat(cellfun(@mean, corrValExpert, 'un', 0));
figure, hold on
bar(xpos - posAdj, mean(naiveMat), barWidth, 'k')
bar(xpos + posAdj, mean(expertMat), barWidth, 'c')
errorbar(xpos - posAdj, mean(naiveMat), sem(naiveMat), 'k', 'lines', 'no')
errorbar(xpos + posAdj, mean(expertMat), zeros(1,length(xpos)), sem(expertMat), 'k', 'lines', 'no')
legend({'Naive', 'Expert'})
xticks(xpos)
xticklabels(featureNames)
xtickangle(45)
ylabel('Correlation')



%% from pre-lick
featureNames = {'Full model', 'noWhisker','noSound','noReward','noWhisking','noLicking','noWhisker&Sound','noWhisker&Reward','noWhisker&Whisking','noWhisker&Licking'};
baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
learnerInd = find(cellfun(@(x) length(x) == 2, sessions));
for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        if si == 1
            data.naive(mi) = load(sprintf('%s%03d\\angle_tuning_model_touchCell_NC_preLick_perTouch_JK%03dS%02d', baseDir, mouse, mouse, session));
%             data1.naive(mi) = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preLick_perTouch_spkOnly_NC', baseDir, mouse, mouse, session));
        else
            ei = find(learnerInd == mi);
            data.expert(ei) = load(sprintf('%s%03d\\angle_tuning_model_touchCell_NC_preLick_perTouch_JK%03dS%02d', baseDir, mouse, mouse, session));
%             data1.expert(ei) = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preLick_perTouch_spkOnly_NC', baseDir, mouse, mouse, session));
        end
    end
end

%%
testInds = [3,16,37:44];
corrValNaive = cell(length(mice),1);
for mi = 1 : length(mice)
    indTuned = find(data.naive(mi).tunedAllCell(:,1));
    indTemp = find(data.naive(mi).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrValsNaive = zeros(length(indWhisker), length(testInds));
    for ci = 1 : length(indWhisker)
        for i = 1 : length(testInds)
            ti = testInds(i);
            if data.naive(mi).tunedAllCell(indWhisker(ci),ti)
                tempCorr = corr(cellfun(@mean, data.naive(mi).spkValAllCell{indWhisker(ci),1}), cellfun(@mean, data.naive(mi).spkValAllCell{indWhisker(ci),ti}));
                if isnan(tempCorr) || tempCorr < 0
                    corrValsNaive(ci,i) = 0;
                else
                    corrValsNaive(ci,i) = tempCorr;
                end
            else
                corrValsNaive(ci,i) = 0;
            end
        end
    end
    corrValNaive{mi} = corrValsNaive;
end

corrValExpert = cell(length(learnerInd),1);
for mi = 1 : length(learnerInd)
    indTuned = find(data.expert(mi).tunedAllCell(:,1));
    indTemp = find(data.expert(mi).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    corrValsNaive = zeros(length(indWhisker), length(testInds));
    for ci = 1 : length(indWhisker)
        for i = 1 : length(testInds)
            ti = testInds(i);
            if data.expert(mi).tunedAllCell(indWhisker(ci),ti)
                tempCorr = corr(cellfun(@mean, data.expert(mi).spkValAllCell{indWhisker(ci),1}), cellfun(@mean, data.expert(mi).spkValAllCell{indWhisker(ci),ti}));
                if isnan(tempCorr) || tempCorr < 0
                    corrValsNaive(ci,i) = 0;
                else
                    corrValsNaive(ci,i) = tempCorr;
                end
            else
                corrValsNaive(ci,i) = 0;
            end
        end
    end
    corrValExpert{mi} = corrValsNaive;
end

%%
posAdj = 0.2;
barWidth = 0.4;
xpos = 1 : length(featureNames);
naiveMat = cell2mat(cellfun(@mean, corrValNaive, 'un', 0));
expertMat = cell2mat(cellfun(@mean, corrValExpert, 'un', 0));
figure, hold on
bar(xpos - posAdj, mean(naiveMat), barWidth, 'k')
bar(xpos + posAdj, mean(expertMat), barWidth, 'c')
errorbar(xpos - posAdj, mean(naiveMat), sem(naiveMat), 'k', 'lines', 'no')
errorbar(xpos + posAdj, mean(expertMat), zeros(1,length(xpos)), sem(expertMat), 'k', 'lines', 'no')
legend({'Naive', 'Expert'})
xticks(xpos)
xticklabels(featureNames)
xtickangle(45)
ylabel('Correlation')


%%
posAdj = 0.2;
barWidth = 0.4;
xpos = 1 : length(featureNames);
naiveMat = cell2mat(cellfun(@mean, corrValNaive(learnerInd), 'un', 0));
expertMat = cell2mat(cellfun(@mean, corrValExpert, 'un', 0));
figure, hold on
bar(xpos - posAdj, mean(naiveMat), barWidth, 'k')
bar(xpos + posAdj, mean(expertMat), barWidth, 'c')
errorbar(xpos - posAdj, mean(naiveMat), sem(naiveMat), 'k', 'lines', 'no')
errorbar(xpos + posAdj, mean(expertMat), zeros(1,length(xpos)), sem(expertMat), 'k', 'lines', 'no')
legend({'Naive', 'Expert'})
xticks(xpos)
xticklabels(featureNames)
xtickangle(45)
ylabel('Correlation')



%% An example from pre-lick negShift test
naive = load('angle_tuning_model_touchCell_NC_preLick_negShift_perTouch_JK025S04');
expert = load('angle_tuning_model_touchCell_NC_preLick_negShift_perTouch_JK025S19');

featureNames = {'Full model', 'noWhisker','noSound','noReward','noWhisking','noLicking','noWhisker&Sound','noWhisker&Reward','noWhisker&Whisking','noWhisker&Licking'};
indTuned = find(naive.tunedAllCell(:,1));
indTemp = find(naive.tunedAllCell(indTuned,3));
indWhisker = indTuned(indTemp);
testInds = [3,16,37:44];
corrValsNaive = zeros(length(indWhisker), length(testInds));
for ci = 1 : length(indWhisker)
    for ti = 1 : length(testInds)
        if naive.tunedAllCell(indWhisker(ci), testInds(ti))
            tempVal = corr(cellfun(@mean, naive.spkValAllCell{indWhisker(ci), 1}), cellfun(@mean, naive.spkValAllCell{indWhisker(ci), testInds(ti)}));
%             if isfinite(tempVal) && tempVal > 0
                corrValsNaive(ci,ti) = tempVal;
%             end
        end
    end
end

indTuned = find(expert.tunedAllCell(:,1));
indTemp = find(expert.tunedAllCell(indTuned,3));
indWhisker = indTuned(indTemp);
testInds = [3,16,37:44];
corrValsExpert = zeros(length(indWhisker), length(testInds));
for ci = 1 : length(indWhisker)
    for ti = 1 : length(testInds)
        if expert.tunedAllCell(indWhisker(ci), testInds(ti))
            tempVal = corr(cellfun(@mean, expert.spkValAllCell{indWhisker(ci), 1}), cellfun(@mean, expert.spkValAllCell{indWhisker(ci), testInds(ti)}));
%             if isfinite(tempVal) && tempVal > 0
                corrValsExpert(ci,ti) = tempVal;
%             end
        end
    end
end

posAdj = 0.2;
barWidth = 0.4;
xpos = 1 : length(testInds);
figure, hold on
bar(xpos - posAdj, mean(corrValsNaive), barWidth, 'k')
bar(xpos + posAdj, mean(corrValsExpert), barWidth, 'c')
errorbar(xpos - posAdj, mean(corrValsNaive), zeros(1,length(xpos)), sem(corrValsNaive), 'k', 'lines', 'no')
errorbar(xpos + posAdj, mean(corrValsExpert), zeros(1,length(xpos)), sem(corrValsNaive), 'k', 'lines', 'no')
