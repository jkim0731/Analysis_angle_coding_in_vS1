% purpose: to look at changes in wkv during subsequent touches, vs object
% angle

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK037'},{'JK038'},{'JK039'},{'JK041'},{'JK052'},{'JK053'},{'JK054'},{'JK056'}];

% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'Y:\Whiskernas\JK\suite2p\';
nFeatures = 11; % 11 features (wkv excluding touch count)
featureNames = {'\Delta\theta', '\Delta\phi', '\Delta\kappa_H', '\Delta\kappa_V', 'Slide distance', 'touch duration', ...
    '\theta', '\phi', '\kappa_H', '\kappa_V', 'Arc length'};
nMice = length(mns);
wkvChanges = cell(nMice,1);
maxTouchNum = zeros(nMice,1);
for mi = 1 : nMice
%     wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]);
    load([baseDir, mns{mi}(3:5), filesep, 'Uber', mns{mi}, sns{mi}, '_NC'])
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
    touchTrialInd = find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
    
    featureMat = cell(nFeatures, 1);
    featureMat{1} = cellfun(@(x) x.protractionTouchDThetaByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{3} = cellfun(@(x) x.protractionTouchDKappaHByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{4} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{5} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{6} = cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    
    featureMat{7} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{8} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{9} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{10} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{11} = cellfun(@(x) x.arcLength(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    
    angles = unique(cellfun(@(x) x.angle, u.trials(touchTrialInd)));
%     angles = unique(cellfun(@(x) x.poleAngle, wfa.trials(touchTrialInd)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.angle == angles(ai), u.trials(touchTrialInd)));
%         angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialInd)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            wkvChanges{mi}{ai,fi} = cell2mat(featureCell);
        end
    end
end

%
colors = jet(length(angles));
maxTn = max(maxTouchNum);

figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
for fi = 1 : nFeatures
    subplot(3,4,fi), hold on
    
%     for ai = 1 : length(angles)
%         tempMat = nan(length(mns), maxTn);
%         for mi = 1 : nMice
%             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
%         end        
% %         errorbar(mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), '-', 'color', colors(ai,:))
%         shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
% %         plot(mean(tempMat), '-', 'color', colors(ai,:))
%     end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==8
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles, ones(length(angles),1)), 'uniformoutput', false), 'location', 'northeastoutside')

%% Results:
%% (1) theta and phi at touch changes during subsequent touches -> They might discriminate the angle by gathering information from multiple touches
%%      -> Test this by looking at # of touches before answer lick, first lick, putting dphi dtheta(across touches)/dphi(across touches) as an input in glm.
%% (2) at 135 degrees, theta stays the same when phi goes up... how could this be happening? (It happens only in expert mice)
%%      -> Figure out this first. Somethings fishy.




%% Relative change in every touch.

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK037'},{'JK038'},{'JK039'},{'JK041'},{'JK052'},{'JK053'},{'JK054'},{'JK056'}];

% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
% sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'Y:\Whiskernas\JK\suite2p\';
nFeatures = 11; % 11 features (wkv excluding touch count)
featureNames = {'\Delta\theta', '\Delta\phi', '\Delta\kappa_H', '\Delta\kappa_V', 'Slide distance', 'touch duration', ...
    '\theta', '\phi', '\kappa_H', '\kappa_V', 'Arc length'};
nMice = length(mns);
wkvChanges = cell(nMice,1);
maxTouchNum = zeros(nMice,1);
for mi = 1 : nMice
%     wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]);
    load([baseDir, mns{mi}(3:5), filesep, 'Uber', mns{mi}, sns{mi}, '_NC'])
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
    touchTrialInd = find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
    
    featureMat = cell(nFeatures, 1);
    featureMat{1} = cellfun(@(x) x.protractionTouchDThetaByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{3} = cellfun(@(x) x.protractionTouchDKappaHByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{4} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{5} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{6} = cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials(touchTrialInd), 'uniformoutput', false);
    
    featureMat{7} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{8} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{9} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{10} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{11} = cellfun(@(x) x.arcLength(cellfun(@(y) y(1)-1, x.protractionTouchChunksByWhisking))', u.trials(touchTrialInd), 'uniformoutput', false);
    
    angles = unique(cellfun(@(x) x.angle, u.trials(touchTrialInd)));
%     angles = unique(cellfun(@(x) x.poleAngle, wfa.trials(touchTrialInd)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.angle == angles(ai), u.trials(touchTrialInd)));
%         angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialInd)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvMat = tempMat - tempMat(:,1);
            wkvChanges{mi}{ai,fi} = wkvMat;
        end
    end
end

%
colors = jet(length(angles));
maxTn = max(maxTouchNum);

figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
for fi = 1 : nFeatures
    subplot(3,4,fi), hold on
    
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end        
%         errorbar(mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), '-', 'color', colors(ai,:))
        shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
%         plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==8
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles, ones(length(angles),1)), 'uniformoutput', false), 'location', 'northeastoutside')




%% Look at other cases - expert and radial distance

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK037'},{'JK038'},{'JK039'},{'JK041'},{'JK052'},{'JK053'},{'JK054'},{'JK056'}];

sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
% sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
nFeatures = 11; % 11 features (wkv excluding touch count)
featureNames = {'\Delta\theta', '\Delta\phi', '\Delta\kappa_H', '\Delta\kappa_V', 'Slide distance', 'touch duration', ...
    '\theta', '\phi', '\kappa_H', '\kappa_V', 'Arc length'};
nMice = length(mns);
wkvChanges = cell(nMice,1);
maxTouchNum = zeros(nMice,1);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]);
    wfa.trials = wfa.trials';
%     load([baseDir, mns{mi}(3:5), filesep, 'Uber', mns{mi}, sns{mi}, '_NC'])
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialInd = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    
    featureMat = cell(nFeatures, 1);
    featureMat{1} = cellfun(@(x) cellfun(@(y) x.theta(y(find(nanmax(abs(x.theta(y)-x.theta(y(1))))==abs(x.theta(y)-x.theta(y(1))),1)))-x.theta(y(1)), ...
        x.protractionTFchunksByWhisking(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)))))), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) cellfun(@(y) x.phi(y(find(nanmax(abs(x.phi(y)-x.phi(y(1))))==abs(x.phi(y)-x.phi(y(1))),1)))-x.phi(y(1)), ...
        x.protractionTFchunksByWhisking(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)))))), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{3} = cellfun(@(x) cellfun(@(y) x.kappaH(y(find(nanmax(abs(x.kappaH(y)-x.kappaH(y(1))))==abs(x.kappaH(y)-x.kappaH(y(1))),1)))-x.kappaH(y(1)), ...
        x.protractionTFchunksByWhisking(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)))))), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{4} = cellfun(@(x) cellfun(@(y) x.kappaV(y(find(nanmax(abs(x.kappaV(y)-x.kappaV(y(1))))==abs(x.kappaV(y)-x.kappaV(y(1))),1)))-x.kappaV(y(1)), ...
        x.protractionTFchunksByWhisking(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunksByWhisking)))))), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{5} = cellfun(@(x) cellfun(@(y) nanmax(y), x.protractionSlideByWhisking), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{6} = cellfun(@(x) x.protractionTouchDurationByWhisking, wfa.trials(touchTrialInd), 'uniformoutput', false);
    
    featureMat{7} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{8} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{9} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{10} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{11} = cellfun(@(x) x.arcLength(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    
%     angles = unique(cellfun(@(x) x.angle, wfa.trials(touchTrialInd)));
    angles = unique(cellfun(@(x) x.poleAngle, wfa.trials(touchTrialInd)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
%         angleInd = find(cellfun(@(x) x.angle == angles(ai), u.trials(touchTrialInd)));
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialInd)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChanges{mi}{ai,fi} = tempMat;
%             wkvMat = tempMat - tempMat(:,1);
%             wkvChanges{mi}{ai,fi} = wkvMat;
        end
    end
end

%
colorsList = jet(7);
colors = colorsList([1,7],:);
% colors = colorsList([1:7],:);
maxTn = max(maxTouchNum);

figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
for fi = 1 : nFeatures
    subplot(3,4,fi), hold on
    
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end        
%         errorbar(mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), '-', 'color', colors(ai,:))
        shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
%         plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==8
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles, ones(length(angles),1)), 'uniformoutput', false), 'location', 'northeastoutside')


%%
%% Result: the trend in theta/phi looks the same in other cases
%%





%% How does the touch height on the pole change?
% calculated from the matching wst time point, whiskerPoleIntersection{:,2}(2)

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
% sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
nMice = length(mns);
maxTouchNum = zeros(nMice,1);
touchHeight = cell(nMice,1);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';
    wsa = Whisker.WhiskerSignalTrialArray_2pad([baseDir, mns{mi}, sns{mi}]); % to get the touch height information
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
    touchTrialIndWS = find(ismember(wsa.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, wsa.trialNums(touchTrialIndWS))); % assume sorted
%     touchOnsetInds = cellfun(@(x) cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialInd));
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsWS = cellfun(@(x,y) cellfun(@(z) find(abs(x.time{2}-z) == nanmin(abs(x.time{2} - z)),1), y, 'uniformoutput', false), wsa.trials(touchTrialIndWS), touchOnsetTimes', 'uniformoutput', false);
    touchHeightValid = cellfun(@(x,y) cellfun(@(z) length(x.whiskerPoleIntersection{z,2}), y), wsa.trials(touchTrialIndWS), touchOnsetIndsWS, 'uniformoutput', false);
    validInd = find(cellfun(@(x) isempty(find(x==0)), touchHeightValid));
    touchHeightCell = cellfun(@(x,y) cellfun(@(z) double(x.whiskerPoleIntersection{z,2}(2))/x.pxPerMm, y), wsa.trials(touchTrialIndWS(validInd)), touchOnsetIndsWS(validInd), 'uniformoutput', false);
%%    
    angles = unique(cellfun(@(x) x.angle, wsa.trials(touchTrialIndWS)));
    touchHeight{mi}  = cell(length(angles), 1); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.angle == angles(ai), wsa.trials(touchTrialIndWS(validInd))));
        touchHeight{mi}{ai} = nan(length(angleInd),maxTouchNum(mi));
        
        featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], touchHeightCell(angleInd)', 'uniformoutput', false);
        tempMat = cell2mat(featureCell);
        touchHeight{mi}{ai} = tempMat;
        
    end
end

%
colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in the length of ''angles''')
end
maxTn = max(maxTouchNum);

figure
hold on
for ai = 1 : length(angles)
    tempMat = nan(length(mns), maxTn);
    for mi = 1 : nMice
        tempMat(mi,1:maxTouchNum(mi)) = nanmean(touchHeight{mi}{ai} - touchHeight{mi}{ai}(:,1));
    end        
    shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
end
for ai = 1 : length(angles)
    tempMat = nan(length(mns), maxTn);
    for mi = 1 : nMice
        tempMat(mi,1:maxTouchNum(mi)) = nanmean(touchHeight{mi}{ai} - touchHeight{mi}{ai}(:,1));
    end        
    plot(1:size(tempMat,2), mean(tempMat), 'color',colors(ai,:))
end
xlabel('Touch order')
ylabel('Relative height (mm)')
title('Change in touch height on the pole')

%%
%% Results: Height works as expected from theta change. (goes up in 45 degrees, stays the same in 135 degrees)
%%

%% How does the follicle height change?

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
% sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
nMice = length(mns);
maxTouchNum = zeros(nMice,1);
baseHeight = cell(nMice,1);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';
    w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch height information
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
    touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);
    baseHeightCell = cellfun(@(x,y) cellfun(@(z) x.base(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3), touchOnsetIndsW3, 'uniformoutput', false);
%%    
    angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
    baseHeight{mi}  = cell(length(angles), 1); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.servoAngle == angles(ai), w3a.trials(touchTrialIndW3)));
        baseHeight{mi}{ai} = nan(length(angleInd),maxTouchNum(mi));
        
        featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], baseHeightCell(angleInd)', 'uniformoutput', false);
        tempMat = cell2mat(featureCell);
        baseHeight{mi}{ai} = tempMat;
    end
end

%
colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in the length of ''angles''')
end
maxTn = max(maxTouchNum);

figure
hold on
for ai = 1 : length(angles)
    tempMat = nan(length(mns), maxTn);
    for mi = 1 : nMice
        tempMat(mi,1:maxTouchNum(mi)) = nanmean(baseHeight{mi}{ai} - baseHeight{mi}{ai}(:,1));
    end        
    shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
end
for ai = 1 : length(angles)
    tempMat = nan(length(mns), maxTn);
    for mi = 1 : nMice
        tempMat(mi,1:maxTouchNum(mi)) = nanmean(baseHeight{mi}{ai} - baseHeight{mi}{ai}(:,1));
    end        
    plot(1:size(tempMat,2), mean(tempMat), 'color',colors(ai,:))
end
xlabel('Touch order')
ylabel('Relative height (mm)')
title('Change in whisker base height')


%%
%% Results: Whisker base height remains the same at 135°, while it goes up at 45°.
%% It suggests their differential control of whisker pad in different pole angle. (to maintain the same height of the base while elevating the base angle requires the whole whisker pad to go down)
%%



%% How does it look like in each touch, not divided by whisking?

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK037'},{'JK038'},{'JK039'},{'JK041'},{'JK052'},{'JK053'},{'JK054'},{'JK056'}];

% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
nFeatures = 11; % 11 features (wkv excluding touch count)
featureNames = {'\Delta\theta', '\Delta\phi', '\Delta\kappa_H', '\Delta\kappa_V', 'Slide distance', 'touch duration', ...
    '\theta', '\phi', '\kappa_H', '\kappa_V', 'Arc length'};
nMice = length(mns);
wkvChanges = cell(nMice,1);
maxTouchNum = zeros(nMice,1);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]);
    wfa.trials = wfa.trials';
%     load([baseDir, mns{mi}(3:5), filesep, 'Uber', mns{mi}, sns{mi}, '_NC'])
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunks), wfa.trials));
    touchTrialInd = find(cellfun(@(x) length(x.protractionTFchunks), wfa.trials));
    
    featureMat = cell(nFeatures, 1);
    featureMat{1} = cellfun(@(x) cellfun(@(y) x.theta(y(find(nanmax(abs(x.theta(y)-x.theta(y(1))))==abs(x.theta(y)-x.theta(y(1))),1)))-x.theta(y(1)), ...
        x.protractionTFchunks( intersect(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunks)))), find(isfinite(x.kappaH(cellfun(@(z) z(1), x.protractionTFchunks))))) )), ...
        wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) cellfun(@(y) x.phi(y(find(nanmax(abs(x.phi(y)-x.phi(y(1))))==abs(x.phi(y)-x.phi(y(1))),1)))-x.phi(y(1)), ...
        x.protractionTFchunks( intersect(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunks)))), find(isfinite(x.kappaH(cellfun(@(z) z(1), x.protractionTFchunks))))) )), ...
        wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{3} = cellfun(@(x) cellfun(@(y) x.kappaH(y(find(nanmax(abs(x.kappaH(y)-x.kappaH(y(1))))==abs(x.kappaH(y)-x.kappaH(y(1))),1)))-x.kappaH(y(1)), ...
        x.protractionTFchunks( intersect(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunks)))), find(isfinite(x.kappaH(cellfun(@(z) z(1), x.protractionTFchunks))))) )), ...
        wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{4} = cellfun(@(x) cellfun(@(y) x.kappaV(y(find(nanmax(abs(x.kappaV(y)-x.kappaV(y(1))))==abs(x.kappaV(y)-x.kappaV(y(1))),1)))-x.kappaV(y(1)), ...
        x.protractionTFchunks( intersect(find(isfinite(x.theta(cellfun(@(z) z(1), x.protractionTFchunks)))), find(isfinite(x.kappaH(cellfun(@(z) z(1), x.protractionTFchunks))))) )), ...
        wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{5} = cellfun(@(x) cellfun(@(y) nanmax(y), x.protractionSlide), wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{6} = cellfun(@(x) x.protractionTouchDuration, wfa.trials(touchTrialInd), 'uniformoutput', false);
    
    featureMat{7} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunks))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{8} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunks))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{9} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunks))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{10} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunks))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    featureMat{11} = cellfun(@(x) x.arcLength(cellfun(@(y) y(1)-1, x.protractionTFchunks))', wfa.trials(touchTrialInd), 'uniformoutput', false);
    
%     angles = unique(cellfun(@(x) x.angle, wfa.trials(touchTrialInd)));
    angles = unique(cellfun(@(x) x.poleAngle, wfa.trials(touchTrialInd)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
%         angleInd = find(cellfun(@(x) x.angle == angles(ai), u.trials(touchTrialInd)));
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialInd)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChanges{mi}{ai,fi} = tempMat;
%             wkvMat = tempMat - tempMat(:,1);
%             wkvChanges{mi}{ai,fi} = wkvMat;
        end
    end
end

%
colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in length of angles')
end
maxTn = max(maxTouchNum);

figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
for fi = 1 : nFeatures
    subplot(3,4,fi), hold on
    
%     for ai = 1 : length(angles)
%         tempMat = nan(length(mns), maxTn);
%         for mi = 1 : nMice
%             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
%         end        
%         shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
%     end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==8
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles, ones(length(angles),1)), 'uniformoutput', false), 'location', 'northeastoutside')

%%
%% Results: It's not because of splitting touches into different whiskings
%%

%% Now, show the summarizing figure of theta, phi, kappaV and base height.


mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
nMice = length(mns);
maxTouchNum = zeros(nMice,1);
baseHeight = cell(nMice,1);
nFeatures = 7; % theta, phi, base height, kappaV, touch height
featureNames = {'\theta', '\phi', 'Whisker base height','\kappa_V', 'Touch pole height', 'Touch pole position', '\kappa_H'};
maxTouchOrder = 10;
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';

    w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
    touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);
    
    featureMat{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMat{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMat{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    
    featureMat{3} = cellfun(@(x,y) cellfun(@(z) x.base(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    featureMat{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    featureMat{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%%    
    angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChanges{mi}{ai,fi} = tempMat;
        end
    end
end

%%
colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in the length of ''angles''')
end
maxTn = max(maxTouchNum);

% figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
figure,
for fi = 1 : nFeatures
    subplot(3,3,fi), hold on
    
%     for ai = 1 : length(angles)
%         tempMat = nan(length(mns), maxTn);
%         for mi = 1 : nMice
%             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi} - wkvChanges{mi}{ai,fi}(:,1));
%         end
%         tempMat = tempMat(:,1:maxTouchOrder);
%         shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
%     end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi} - wkvChanges{mi}{ai,fi}(:,1));
        end
        tempMat = tempMat(:,1:maxTouchOrder);
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==7
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles', ones(1,length(angles))), 'uniformoutput', false), 'location', 'northeastoutside')


%% Find the best representing mouse 
%% From matching phi, theta, and whisker base height in 45 and 135 degrees of both naive and expert 7 angles

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];

sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
nMice = length(mns);
maxTouchNum = zeros(nMice,1);
baseHeight = cell(nMice,1);
nFeatures = 3; % theta, phi, base height
% featureNames = {'\theta', '\phi', 'Whisker base height','\kappa_V', 'Touch pole height', 'Touch pole position', '\kappa_H'};
maxTouchOrder = 10;
wkvChanges = cell(1,nMice);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';

    w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
    touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);
    
    featureMat{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
%     featureMat{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
%     featureMat{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    
    featureMat{3} = cellfun(@(x,y) cellfun(@(z) x.base(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%     featureMat{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%     featureMat{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%%    
    angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChanges{mi}{ai,fi} = tempMat(:,1:maxTouchOrder);
        end
    end
end
%%


templateNaiveEach = cellfun(@(x) cellfun(@nanmean, [x(1,:),x(7,:)], 'uniformoutput', false), wkvChanges, 'uniformoutput', false);
templateNaive = zeros(maxTouchOrder, nFeatures*2);
for i = 1 : nFeatures*2
    templateNaive(:,i) = mean(cell2mat(cellfun(@(x) x{i}, templateNaiveEach', 'uniformoutput', false)))';
end

%%
mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];

sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles

wkvChangesExpert = cell(1,nMice);
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';

    w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
    touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);
    
    featureMat{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMat{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
%     featureMat{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
%     featureMat{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    
    featureMat{3} = cellfun(@(x,y) cellfun(@(z) x.base(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%     featureMat{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%     featureMat{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
%%    
    angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
    wkvChangesExpert{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
        for fi = 1 : nFeatures
            wkvChangesExpert{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChangesExpert{mi}{ai,fi} = tempMat(:,1:maxTouchOrder);
        end
    end
end
%%

templateExpertEach = cellfun(@(x) cellfun(@nanmean, [x(1,:),x(7,:)], 'uniformoutput', false), wkvChangesExpert, 'uniformoutput', false);
templateExpert = zeros(maxTouchOrder, nFeatures*2);
for i = 1 : nFeatures*2
    templateExpert(:,i) = mean(cell2mat(cellfun(@(x) x{i}, templateExpertEach', 'uniformoutput', false)))';
end


%% same naive and expert mouse

templateAll = [templateNaive, templateExpert];
templateAllEach = cellfun(@(x,y) [x,y], templateNaiveEach, templateExpertEach, 'uniformoutput', false);

corrVal = zeros(nMice,nFeatures*2*2);
for i = 1 : nFeatures*2*2
    
    temp = corr(templateAll(:,i), cell2mat(cellfun(@(x) x{i}', templateAllEach, 'uniformoutput', false)));
    corrVal(:,i) = (temp-min(temp))./(max(temp)-min(temp));
end

[~, corrSortedInd] = sort(sum(corrVal,2),'descend');

%% result: 3,2,6,1,4,5

%% only from naive

corrValNaive = zeros(nMice,nFeatures*2);
for i = 1 : nFeatures*2
    
    temp = corr(templateNaive(:,i), cell2mat(cellfun(@(x) x{i}', templateNaiveEach, 'uniformoutput', false)));
    corrValNaive(:,i) = (temp-min(temp))./(max(temp)-min(temp));
end
[~, corrNaiveSortedInd] = sort(sum(corrValNaive,2),'descend');

%% result: 3,1,6,2,5,4


%% only from expert

corrValExpert = zeros(nMice,nFeatures*2);
for i = 1 : nFeatures*2
    
    temp = corr(templateExpert(:,i), cell2mat(cellfun(@(x) x{i}', templateExpertEach, 'uniformoutput', false)));
    corrValExpert(:,i) = (temp-min(temp))./(max(temp)-min(temp));
end
[~, corrExpertSortedInd] = sort(sum(corrValExpert,2),'descend');

%% result: 2,4,3,6,1,5


%% show summary graphs

mi = 3;
figure
subplot(231), hold on
plot(templateNaiveEach{mi}{1}, 'b-')
plot(templateNaiveEach{mi}{4}, 'r-')
title('\theta')

subplot(232), hold on
plot(templateNaiveEach{mi}{2}, 'b-')
plot(templateNaiveEach{mi}{5}, 'r-')
title('\phi')

subplot(233), hold on
plot(templateNaiveEach{mi}{3}, 'b-')
plot(templateNaiveEach{mi}{6}, 'r-')
title('whisker base height')


subplot(234), hold on
plot(templateExpertEach{mi}{1}, 'b-')
plot(templateExpertEach{mi}{4}, 'r-')

subplot(235), hold on
plot(templateExpertEach{mi}{2}, 'b-')
plot(templateExpertEach{mi}{5}, 'r-')

subplot(236), hold on
plot(templateExpertEach{mi}{3}, 'b-')
plot(templateExpertEach{mi}{6}, 'r-')


%% Show video and consecutive touch kinematics of whisker in 3D

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance
baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';

mi = 3;

nFeatures = 7; % theta, phi, base height, kappaV, touch height
featureNames = {'\theta', '\phi', 'Whisker base height','\kappa_V', 'Touch pole height', 'Touch pole position', '\kappa_H'};
maxTouchOrder = 10;

wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
wfa.trials = wfa.trials';

w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
maxTouchNum = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);

featureMat{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);

featureMat{3} = cellfun(@(x,y) cellfun(@(z) x.base(z-1,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
featureMat{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
featureMat{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    
angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
angleTn = cell(length(angles),1);

%
wkvChanges = cell(length(angles), nFeatures); 

for ai = 1 : length(angles)
    angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
    angleTn{ai} = wfa.trialNums(touchTrialIndWF(angleInd));
    
    for fi = 1 : nFeatures
        wkvChanges{ai,fi} = nan(length(angleInd),maxTouchNum);
        featureCell = cellfun(@(x) [x, nan(1,maxTouchNum - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
        tempMat = cell2mat(featureCell);
        wkvChanges{ai,fi} = tempMat(:,1:maxTouchOrder);
    end
end
    

%% finding the best correlated trial
% ai = 1; % for 45 degrees
ai = 7; % for 135 degrees
sorti = 3;

tempCell = wkvChanges(ai,:);
templates = cellfun(@nanmean, tempCell, 'uniformoutput', false);
% %%
% figure,
% for i = 1 : 7
%     subplot(3,3,i), plot(templates{i})
% end
% %%
corrValNorm = zeros(size(tempCell{1},1), nFeatures);
for fi = 1 : nFeatures
    tempMat = tempCell{fi}';
    tempMat(isnan(tempMat)) = 0;
    tempCorr = corr(templates{fi}', tempMat);
    corrValNorm(:,fi) = (tempCorr - min(tempCorr)) ./ (max(tempCorr) - min(tempCorr));
end

[~, indSorted] = sort(sum(corrValNorm,2), 'descend');
tn = angleTn{ai}(indSorted(sorti));
w3ind = find(w3a.trialNums == tn);
w3 = w3a.trials{w3ind};
wfind = find(wfa.trialNums == tn);
wf = wfa.trials{wfind};
%
touchTimes = wf.time(cell2mat(cellfun(@(x) x', wf.protractionTFchunks, 'uniformoutput', false)));
touchOnsetTimes = wf.time(cellfun(@(x) x(1), wf.protractionTFchunksByWhisking));
touchIndsW3 = find(ismember(w3.time, touchTimes));
touchOnsetIndsW3 = find(ismember(w3.time, touchOnsetTimes));


% %
% % v = VideoWriter(sprintf('expert%03dd_%02d.avi',angles(ai),sorti),'Uncompressed AVI');
% v = VideoWriter(sprintf('naive%03dd_%02d.avi',angles(ai),sorti),'Uncompressed AVI');
% v.FrameRate = 31;
% open(v)
% 
% maxLim = max(cell2mat(cellfun(@max, w3.trackerData(w3.poleUpFrames), 'uniformoutput', false)));
% minLim = min(cell2mat(cellfun(@min, w3.trackerData(w3.poleUpFrames), 'uniformoutput', false)));
% 
% fh = figure;
% for i = 1 : length(w3.poleUpFrames)
%     clf
%     currFrame = w3.poleUpFrames(i);
%     subplot(221), hold on
%     show_one_3D(w3, currFrame, ismember(currFrame, touchIndsW3))
%     xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
%     
%     subplot(222), hold on
%     show_one_3D(w3, currFrame, ismember(currFrame, touchIndsW3), [0 90])
%     title('Top-down view')
%     xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
%     
%     text(maxLim(1)+10, maxLim(2)-10, [num2str(round( w3.time(currFrame) - w3.time(w3.poleUpFrames(1)) , 2)), ' s'])
% 
%     subplot(223), hold on
%     show_one_3D(w3, currFrame, ismember(currFrame, touchIndsW3), [90 0])
%     title('Front view')
%     xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
% 
%     subplot(224), hold on
%     show_one_3D(w3, currFrame, ismember(currFrame, touchIndsW3), [0 0])
%     title('Lateral view')
%     xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
% 
%     drawnow
%     writeVideo(v,getframe(fh))
% end
% close(v)

%
maxLim = max(cell2mat(cellfun(@max, w3.trackerData(touchOnsetIndsW3), 'uniformoutput', false)));
minLim = min(cell2mat(cellfun(@min, w3.trackerData(touchOnsetIndsW3), 'uniformoutput', false)));
minColorVal = 0.2;
% figure('units','normalized','outerposition', [0.1 0.1 0.66 0.67])
figure,
numTouch = length(touchOnsetIndsW3);
whiskerColor = [linspace(1-minColorVal,0,numTouch); linspace(1-minColorVal,0,numTouch); linspace(1-minColorVal,0,numTouch)]';
baseColor = [linspace(1-minColorVal,minColorVal,numTouch); zeros(1,numTouch); zeros(1,numTouch)]';
touchColor = [zeros(1,numTouch); zeros(1,numTouch); linspace(1-minColorVal,minColorVal,numTouch)]';
for i = 1 : length(touchOnsetIndsW3)
    frameNum = touchOnsetIndsW3(i);    
    x = w3.trackerData{frameNum}(w3.baseInd(frameNum):end,1);
    y = w3.trackerData{frameNum}(w3.baseInd(frameNum):end,2);
    z = w3.trackerData{frameNum}(w3.baseInd(frameNum):end,3);
    
    subplot(221), hold on
    
    plot3(x, y, z, '-', 'color', whiskerColor(i,:))
    plot3(w3.base(frameNum,1), w3.base(frameNum,2), w3.base(frameNum,3), '.', 'markersize', 20, 'color', baseColor(i,:))
    whisker = [x,y,z];
    ind = find(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2) == min(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2)) );
    plot3(x(ind), y(ind), z(ind), '.', 'markersize', 20, 'color', touchColor(i,:))
    axis equal
    xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
    xlabel('A-P'), ylabel('M-L'), zlabel('D-V')
    view(3)
    
    subplot(222), hold on
    
    plot3(x, y, z, '-', 'color', whiskerColor(i,:))
    plot3(w3.base(frameNum,1), w3.base(frameNum,2), w3.base(frameNum,3), '.', 'markersize', 20, 'color', baseColor(i,:))
    whisker = [x,y,z];
    ind = find(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2) == min(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2)) );
    plot3(x(ind), y(ind), z(ind), '.', 'markersize', 20, 'color', touchColor(i,:))
    axis equal
    xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
    xlabel('A-P'), ylabel('M-L'), zlabel('D-V')
    view([0 90]), title('Top-down view')
    
    subplot(223), hold on
    
    plot3(x, y, z, '-', 'color', whiskerColor(i,:))
    plot3(w3.base(frameNum,1), w3.base(frameNum,2), w3.base(frameNum,3), '.', 'markersize', 20, 'color', baseColor(i,:))
    whisker = [x,y,z];
    ind = find(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2) == min(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2)) );
    plot3(x(ind), y(ind), z(ind), '.', 'markersize', 20, 'color', touchColor(i,:))
    axis equal
    xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
    xlabel('A-P'), ylabel('M-L'), zlabel('D-V')
    view([90 0]), title('Front view')
        
    subplot(224), hold on
    
    plot3(x, y, z, '-', 'color', whiskerColor(i,:))
    plot3(w3.base(frameNum,1), w3.base(frameNum,2), w3.base(frameNum,3), '.', 'markersize', 20, 'color', baseColor(i,:))
    whisker = [x,y,z];
    ind = find(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2) == min(sum((whisker-w3.intersectPoint(frameNum,:)).^2,2)) );
    plot3(x(ind), y(ind), z(ind), '.', 'markersize', 20, 'color', touchColor(i,:))
    axis equal
    xlim([minLim(1) maxLim(1)]), ylim([minLim(2) maxLim(2)]), zlim([minLim(3) maxLim(3)])
    xlabel('A-P'), ylabel('M-L'), zlabel('D-V')
    view([0 0]), title('Lateral view')
        
end

% another figure of example trace

colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in the length of ''angles''')
end


figure,
for fi = 1 : nFeatures
    subplot(3,3,fi), hold on
    
    tempMat = wkvChanges{ai,fi}(indSorted(sorti),:);
    
    plot(tempMat, '-', 'color', colors(ai,:))

    title(featureNames{fi})
    if fi==7
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% 
% 
% % check from w3 and wf
% figure,
% subplot(3,3,1), plot(cellfun(@(x) wf.theta(x(1)-1), wf.protractionTFchunksByWhisking(1:maxTouchOrder)), 'r-')
% subplot(3,3,2), plot(cellfun(@(x) wf.phi(x(1)-1), wf.protractionTFchunksByWhisking(1:maxTouchOrder)), 'r-')
% subplot(3,3,3), plot(w3.base(touchOnsetIndsW3(1:maxTouchOrder)-1,3), 'r-')
% subplot(3,3,4), plot(cellfun(@(x) wf.kappaV(x(1)-1), wf.protractionTFchunksByWhisking(1:maxTouchOrder)), 'r-')
% subplot(3,3,5), plot(w3.intersectPoint(touchOnsetIndsW3(1:maxTouchOrder),3), 'r-')
% subplot(3,3,6), plot(w3.intersectPoint(touchOnsetIndsW3(1:maxTouchOrder),1), 'r-')
% subplot(3,3,7), plot(cellfun(@(x) wf.kappaH(x(1)-1), wf.protractionTFchunksByWhisking(1:maxTouchOrder)), 'r-')

%% Is the change in elevation of whisker base angle (phi) related to licking side?

%% First, look at before first lick features.

mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
% sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance

baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';
uDir = 'Y:\Whiskernas\JK\suite2p\';
nMice = length(mns);
maxTouchNum = zeros(nMice,1);
baseHeight = cell(nMice,1);
nFeatures = 7; % theta, phi, base height, kappaV, touch height
featureMatPre = cell(1,nFeatures);
featureMat = cell(1,nFeatures);
featureNames = {'\theta', '\phi', 'Whisker base height','\kappa_V', 'Touch pole height', 'Touch pole position', '\kappa_H'};
maxTouchOrder = 10;
for mi = 1 : nMice
% for mi = 3
    wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
    wfa.trials = wfa.trials';

    w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
    
    mouse = mns{mi}(3:5);
    load(sprintf('%s%s\\Uber%s%s_NC',uDir,mouse,mns{mi}, sns{mi})) % load u
    %%
    
    maxTouchNum(mi) = max(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
    
    touchTrialTotalIndU = find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
%     touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
    touchTrialNum = cellfun(@(x) x.trialNum, u.trials(touchTrialTotalIndU));
    touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndWF = find(ismember(wfa.trialNums, touchTrialNum)); % assume sorted
    touchTrialIndU = find(ismember(u.trialNums, touchTrialNum)); % assume sorted
    touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);
    
    lickTimes = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(touchTrialIndU), 'uniformoutput', false);
    poleUpOnsetTimes = cellfun(@(x) x.poleUpOnsetTime, u.trials(touchTrialIndU), 'uniformoutput', false);
    firstLickTimes = cellfun(@(x,y) x(find(x>y,1)), lickTimes, poleUpOnsetTimes, 'uniformoutput', false);
    noLickInd = find(cellfun(@isempty, firstLickTimes));
    firstLickTimes(noLickInd) = cellfun(@(x) x.poleDownOnsetTime, u.trials(touchTrialIndU(noLickInd)), 'uniformoutput', false);
%     firstLickIndsW3 = cellfun(@(x,y) find(abs(x.time-y) == min(abs(x.time-y)),1), w3a.trials(touchTrialIndW3), firstLickTimes, 'uniformoutput', false);
%     firstLickIndsWF = cellfun(@(x,y) find(abs(x.time-y) == min(abs(x.time-y)),1), wfa.trials(touchTrialIndWF), firstLickTimes, 'uniformoutput', false);
    %%
    featureMatPre{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMatPre{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMatPre{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    featureMatPre{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
    
    featureMatPre{3} = cellfun(@(x,y) cellfun(@(z) x.base(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    featureMatPre{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    featureMatPre{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    
    wfMask = cellfun(@(x,z) [ones(1,sum(x.time(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))<z)), nan(1, sum(x.time(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))>=z))], ...
        wfa.trials(touchTrialIndWF), firstLickTimes, 'uniformoutput', false);
    w3Mask = cellfun(@(x,y,z) [ones(1,sum(x.time(cell2mat(y))<z)), nan(1, sum(x.time(cell2mat(y))>=z))], ...
        w3a.trials(touchTrialIndW3), touchOnsetIndsW3, firstLickTimes', 'uniformoutput', false);
    
    featureMat{1} = cellfun(@(x,y) x.*y, featureMatPre{1}, wfMask, 'uniformoutput', false);
    featureMat{2} = cellfun(@(x,y) x.*y, featureMatPre{2}, wfMask, 'uniformoutput', false);
    featureMat{4} = cellfun(@(x,y) x.*y, featureMatPre{4}, wfMask, 'uniformoutput', false);
    featureMat{7} = cellfun(@(x,y) x.*y, featureMatPre{7}, wfMask, 'uniformoutput', false);
    
    featureMat{3} = cellfun(@(x,y) x.*y, featureMatPre{3}, w3Mask', 'uniformoutput', false);
    featureMat{5} = cellfun(@(x,y) x.*y, featureMatPre{5}, w3Mask', 'uniformoutput', false);
    featureMat{6} = cellfun(@(x,y) x.*y, featureMatPre{6}, w3Mask', 'uniformoutput', false);
%
    angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
    wkvChanges{mi} = cell(length(angles), nFeatures); 
    for ai = 1 : length(angles)
        angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
        for fi = 1 : nFeatures
            wkvChanges{mi}{ai,fi} = nan(length(angleInd),maxTouchNum(mi));
            featureCell = cellfun(@(x) [x, nan(1,maxTouchNum(mi) - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
            tempMat = cell2mat(featureCell);
            wkvChanges{mi}{ai,fi} = tempMat;
        end 
    end
end

%%
colorsList = jet(7);
if length(angles) == 2
    colors = colorsList([1,7],:);
elseif length(angles) == 7
    colors = colorsList([1:7],:);
else
    error('Error in the length of ''angles''')
end
maxTn = max(maxTouchNum);

% figure('units','normalized','outerposition',[0.1 0.1 0.66 0.67])
figure,
for fi = 1 : nFeatures
    subplot(3,3,fi), hold on
    
%     for ai = 1 : length(angles)
%         tempMat = nan(length(mns), maxTn);
%         for mi = 1 : nMice
%             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi} - wkvChanges{mi}{ai,fi}(:,1));
% %             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
%         end
%         tempMat = tempMat(:,1:maxTouchOrder);
%         shadedErrorBar(1:size(tempMat,2), mean(tempMat), std(tempMat)./sqrt(sum(~isnan(tempMat))), 'lineprops', {'color',colors(ai,:)})
%     end
    for ai = 1 : length(angles)
        tempMat = nan(length(mns), maxTn);
        for mi = 1 : nMice
            tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi} - wkvChanges{mi}{ai,fi}(:,1));
%             tempMat(mi,1:maxTouchNum(mi)) = nanmean(wkvChanges{mi}{ai,fi});
        end
        tempMat = tempMat(:,1:maxTouchOrder);
        plot(mean(tempMat), '-', 'color', colors(ai,:))
    end
    title(featureNames{fi})
    if fi==7
        xlabel('Touch order')
    end
%     set(gca,'fontsize',12)
end
% subplot(3,4,nFeatures+1)
legend(cellfun(@(x) [num2str(x), '\circ'], mat2cell(angles', ones(1,length(angles))), 'uniformoutput', false), 'location', 'northeastoutside')


%% See if it is true in other object angles
%% only consider in expert mice, since naive mice have bias or alternation or both.
%% Look at ambiguous angles, 75:105 degrees
%% Use only dedicated lick trials (licking all the same side, > 3 licks before answer lick)
%% Compare it to 135 degrees. (and 45 degrees too)


mns = [{'JK025'},{'JK027'},{'JK030'},{'JK036'},{'JK039'},{'JK052'}];
% sns = [{'S05'},{'S02'},{'S04'},{'S02'},{'S02'},{'S05'}]; %naive
% sns = [{'S18'},{'S07'},{'S20'},{'S16'},{'S21'},{'S20'}]; %expert
% sns = [{'S04'},{'S03'},{'S03'},{'S01'},{'S01'},{'S03'}]; %naive discreteAngles
sns = [{'S19'},{'S10'},{'S21'},{'S17'},{'S23'},{'S21'}]; %expert discreteAngles
% sns = [{'S22'},{'S14'},{'S22'},{'S18'},{'S24'},{'S26'}]; %radial distance
baseDir = 'Y:\Whiskernas\JK\whisker\tracked\';

mi = 1;

nFeatures = 7; % theta, phi, base height, kappaV, touch height
featureNames = {'\theta', '\phi', 'Whisker base height','\kappa_V', 'Touch pole height', 'Touch pole position', '\kappa_H'};
maxTouchOrder = 10;

wfa = Whisker.WhiskerFinal_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the touch information
wfa.trials = wfa.trials';

w3a = Whisker.Whisker3D_2padArray([baseDir, mns{mi}, sns{mi}]); % to get the whisker base height information
maxTouchNum = max(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
touchTrialTotalIndWF = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
touchTrialNum = cellfun(@(x) x.trialNum, wfa.trials(touchTrialTotalIndWF));
touchTrialIndW3 = find(ismember(w3a.trialNums, touchTrialNum)); % assume sorted
touchTrialIndWF = find(ismember(wfa.trialNums, w3a.trialNums(touchTrialIndW3))); % assume sorted
touchOnsetTimes = cellfun(@(x) cellfun(@(y) x.time(y(1)-1), x.protractionTFchunksByWhisking, 'uniformoutput', false), wfa.trials(touchTrialIndWF), 'uniformoutput', false);
touchOnsetIndsW3 = cellfun(@(x,y) cellfun(@(z) find(abs(x.time-z) == nanmin(abs(x.time - z)),1), y, 'uniformoutput', false), w3a.trials(touchTrialIndW3), touchOnsetTimes', 'uniformoutput', false);

featureMat{1} = cellfun(@(x) x.theta(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{2} = cellfun(@(x) x.phi(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{4} = cellfun(@(x) x.kappaV(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);
featureMat{7} = cellfun(@(x) x.kappaH(cellfun(@(y) y(1)-1, x.protractionTFchunksByWhisking))', wfa.trials(touchTrialIndWF), 'uniformoutput', false);

featureMat{3} = cellfun(@(x,y) cellfun(@(z) x.base(z-1,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
featureMat{5} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,3)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
featureMat{6} = cellfun(@(x,y) cellfun(@(z) x.intersectPoint(z,1)/x.pxPerMm, y), w3a.trials(touchTrialIndW3)', touchOnsetIndsW3', 'uniformoutput', false);
    
angles = unique(cellfun(@(x) x.servoAngle, w3a.trials(touchTrialIndW3)));
angleTn = cell(length(angles),1);

%
wkvChanges = cell(length(angles), nFeatures); 

for ai = 1 : length(angles)
    angleInd = find(cellfun(@(x) x.poleAngle == angles(ai), wfa.trials(touchTrialIndWF)));
    angleTn{ai} = wfa.trialNums(touchTrialIndWF(angleInd));
    
    for fi = 1 : nFeatures
        wkvChanges{ai,fi} = nan(length(angleInd),maxTouchNum);
        featureCell = cellfun(@(x) [x, nan(1,maxTouchNum - length(x))], featureMat{fi}(angleInd), 'uniformoutput', false);
        tempMat = cell2mat(featureCell);
        wkvChanges{ai,fi} = tempMat(:,1:maxTouchOrder);
    end
end
    

