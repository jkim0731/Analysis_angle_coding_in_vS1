% data structure to make:
% Whisker kinematic variables. all & divided into planes. Alongside with angles. Better if
% matched with 'vals' of 'angle_tuning_summary'
% averaged in each trial. Multiple touches will have averaged "at touch"
% variables and "during touch" variables
% Use just protraction chunks, and predecision touches.

% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};

clear
    mouse = 27;
    session = 10;

    mouseName = sprintf('JK%03d',mouse);
    sessionName = sprintf('S%02d',session);
    angles = 45:15:135;
    %% Load necessary files
    uberFolder = 'Y:\Whiskernas\JK\suite2p\';        
    ufn = sprintf('UberJK%03dS%02d',mouse, session);
    load([uberFolder, sprintf('%03d',mouse), filesep, ufn])
    
    %%
%     it: touch theta, phi, kappaH, kappaV, radialD, touch count
%     dt: max(dTheta), max(dPhi), max(dKappaH), max(dKappaV), max(touchDuration), max(slideDistance), max(abs(dPhi)), max(abs(dKappaV))
%     In each trial, these values are averaged.
    
    touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunks), u.trials));
    tempData = nan(length(touchTrialInds),15); % 1:6, it. 7:14, dt. 15, angles
    tempData(:,1) = cellfun(@(x) mean(x.theta(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
    tempData(:,2) = cellfun(@(x) mean(x.phi(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
    tempData(:,3) = cellfun(@(x) mean(x.kappaH(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
    tempData(:,4) = cellfun(@(x) mean(x.kappaV(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
    tempData(:,5) = cellfun(@(x) mean(x.arcLength(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
    tempData(:,6) = cellfun(@(x) length(x.protractionTouchChunks), u.trials(touchTrialInds));
    
    tempData(:,7) = cellfun(@(x) mean(x.protractionTouchDTheta), u.trials(touchTrialInds));
    tempData(:,8) = cellfun(@(x) mean(x.protractionTouchDPhi), u.trials(touchTrialInds));
    tempData(:,9) = cellfun(@(x) mean(x.protractionTouchDKappaH), u.trials(touchTrialInds));
    tempData(:,10) = cellfun(@(x) mean(x.protractionTouchDKappaV), u.trials(touchTrialInds));
    tempData(:,11) = cellfun(@(x) mean(x.protractionTouchDuration), u.trials(touchTrialInds));
    tempData(:,12) = cellfun(@(x) mean(x.protractionTouchSlideDistance), u.trials(touchTrialInds));
    tempData(:,13) = cellfun(@(x) mean(x.protractionAbsTouchDPhi), u.trials(touchTrialInds));
    tempData(:,14) = cellfun(@(x) mean(x.protractionAbsTouchDKappaV), u.trials(touchTrialInds));
    
    tempData(:,15) = cellfun(@(x) x.angle, u.trials(touchTrialInds));
    
    % remove rows with Nans
    rmInds = find(isnan(sum(tempData,2)));
    tempData(rmInds,:) = [];
    
    %% Re-order data in the order of angles
    
    [angleVals, I] = sort(tempData(:,end));
    dataAngleOrdered = zeros(size(tempData));
    for ai = 1 : length(angles)
        tempinds = find(angleVals == angles(ai));
        dataAngleOrdered(tempinds,:) = tempData(I(tempinds),:);
    end
    
    
    %% t-SNE
    tsneData = dataAngleOrdered(:,1:end-1);
    Y = tsne(tsneData,'NumDimensions', 3, 'Algorithm','exact','Standardize',true);
    distance = squareform(pdist(Y));
    
    cmap = jet(7);
%%
    figure; hold on;
    for ai = 1 : length(angles)
        tempinds = find(angleVals == angles(ai));
        scatter3(Y(tempinds,1), Y(tempinds,2), Y(tempinds,3), 25, cmap(ai,:), 'filled')
    end

    view(-135, 35);

%%
maxD = max(max(distance)); minD = min(min(distance));
distance = (distance - minD) ./ (maxD - minD);
figure;
imagesc(distance);
set(gca,'YDir','normal');
colormap('hot');
title([mouseName, ' ', sessionName])
xticks(find(diff(angleVals))+0.5)
xticklabels({''})
yticks(find(diff(angleVals))+0.5)
yticklabels({''})
set(gca, 'tickdir', 'out', 'ticklength', [0.03 0], 'box', 'off')

    %% Notes on JK036 S01 and JK039 S01
    
    bDir = 'Y:\Whiskernas\JK\SoloData\';
    
    load([bDir, mouseName, filesep, 'behavior_', mouseName])
    sessionInd = find(cellfun(@(x) strcmp(x.sessionName, sessionName), b));
    bSession = b{sessionInd};
    figure, 
    regTrialInds = find(cellfun(@(x) (x.motorApPosition > 0), bSession.trials));
    plot(cellfun(@(x) x.motorApPosition, bSession.trials(regTrialInds(2:end))))
    xlabel('Trials'), ylabel('AP motor position')
    title([mouseName, ' ', sessionName])
    
    %%

mice = [25,27,30,36,39,52];
sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18],[1,23,24],[3,21,26]};
corrs = cell(3,1); % 1 for position vs theta, 2 for position vs kappaV, 3 for theta vs kappaV
for i = 1 : 3
    corrs{i} = zeros(length(mice),2);
end
for mi = 1 : length(mice)
% for mi = 1
    mouse = mice(mi);
    for si = 1 : 2        
%     for si = 1
        session = sessions{mi}(si);

        mouseName = sprintf('JK%03d',mouse);
        sessionName = sprintf('S%02d',session);

        uberFolder = 'Y:\Whiskernas\JK\suite2p\';        
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load([uberFolder, sprintf('%03d',mouse), filesep, ufn])
        
        bDir = 'Y:\Whiskernas\JK\SoloData\';
        load([bDir, mouseName, filesep, 'behavior_', mouseName])
        sessionInd = find(cellfun(@(x) strcmp(x.sessionName, sessionName), b));
        bSession = b{sessionInd};

        touchTrialInds = find(cellfun(@(x) length(x.protractionTouchChunks), u.trials));
        touchTrialNums = cellfun(@(x) x.trialNum, u.trials(touchTrialInds));
        bTrialInds = find(cellfun(@(x) ismember(x.trialNum, touchTrialNums), bSession.trials));
        
        motorPositions = cellfun(@(x) x.motorApPosition, bSession.trials(bTrialInds));
        meanThetas = cellfun(@(x) mean(x.theta(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
        meanKappaVs = cellfun(@(x) mean(x.kappaV(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds));
        
        mpNanInd = find(isnan(motorPositions));
        mtNanInd = find(isnan(meanThetas));
        mkNanInd = find(isnan(meanKappaVs));
        nanInds = union(union(mpNanInd, mtNanInd), mkNanInd);
        keepInds = setdiff(1:length(motorPositions), nanInds);
%         %%
%         figure, hold on
%         
%         h = plot(cellfun(@(x) x.motorApPosition, bSession.trials(bTrialInds)));
%         xlabel('Trials'), ylabel('AP motor position')
%         yyaxis right
%         plot(cellfun(@(x) mean(x.theta(cellfun(@(y) y(1), x.protractionTouchChunks))), u.trials(touchTrialInds)))
%         ylabel('Mean touch \theta')
%         c = get(h,'Color');
%         ax = gca;
%         ax.YAxis(1).Color = c;        
%         title([mouseName, ' ', sessionName])
        corrs{1}(mi,si) = corr(motorPositions(keepInds)', meanThetas(keepInds));
        corrs{2}(mi,si) = corr(motorPositions(keepInds)', meanKappaVs(keepInds));
        corrs{3}(mi,si) = corr(meanThetas(keepInds), meanKappaVs(keepInds));

    end
end

%% Results
% Anterior-posterior pole position distribution is different from session to session (their shape),
% and it affects mean theta at touch which is highly correlated with ap pole position and not being used for angle discrimination.
% So, it will be better to see the distance between trials with the parameters that are actually being used for angle discrimination.
% Previously, Jon showed it, but only in expert mice.
% Do the same analysis in naive mice.
% Then compare the classifier performance. Collect paramters important in that classification. Look at high-dimensional distance (or in 3D t-SNE space).
% Try other different kinds of classification methods (LDA, QDA, KNN, SVM, trees, ensembles, etc).

%% Classification using multinomial glm
%% Comparing between naive and expert

%% Results: similar between naive and expert. Only top 4 features were good enough (slide distance, dKappaV, dKappaH, and dPhi)
% but only in trial-averaged kinematics. performance was similar between
% naive and expert even for individual touches, but top 5 was not close in
% either cases.

%% distance between angles from top 4 features.
mice = [25,27,30,36,39,52];
sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[1,23,24],[3,21,26]};
% for mi = 1 : length(mice)
for mi = 2
%     for si = 1 : 2
    for si = 2    
        mouse = mice(mi);
        session = sessions{mi}(si);

        mouseName = sprintf('JK%03d',mouse);
        sessionName = sprintf('S%02d',session);
        angles = 45:15:135;
        % Load necessary files
        uberFolder = 'Y:\Whiskernas\JK\suite2p\';        
        ufn = sprintf('Uber%s%s',mouseName, sessionName);
        load([uberFolder, sprintf('%03d',mouse), filesep, ufn])

        touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunks), u.trials));
        tempData = nan(length(touchTrialInds),5); % 

        tempData(:,1) = cellfun(@(x) mean(x.protractionTouchDPhi), u.trials(touchTrialInds));
        tempData(:,2) = cellfun(@(x) mean(x.protractionTouchDKappaH), u.trials(touchTrialInds));
        tempData(:,3) = cellfun(@(x) mean(x.protractionTouchDKappaV), u.trials(touchTrialInds));
        tempData(:,4) = cellfun(@(x) mean(x.protractionTouchSlideDistance), u.trials(touchTrialInds));

        tempData(:,5) = cellfun(@(x) x.angle, u.trials(touchTrialInds));

        % remove rows with Nans
        rmInds = find(isnan(sum(tempData,2)));
        tempData(rmInds,:) = [];

        % Re-order data in the order of angles

        [angleVals, I] = sort(tempData(:,end));
        dataAngleOrdered = zeros(size(tempData));
        for ai = 1 : length(angles)
            tempinds = find(angleVals == angles(ai));
            dataAngleOrdered(tempinds,:) = tempData(I(tempinds),:);
        end

        % t-SNE
        tsneData = dataAngleOrdered(:,1:end-1);
        Y = tsne(tsneData,'NumDimensions', 3, 'Algorithm','exact','Standardize',true);
        
%
%
%
        distance = squareform(pdist(Y));
%         meanTsneData = mean(tsneData);
%         stdTsneData = std(tsneData);
%         tsneData = (tsneData - meanTsneData)./stdTsneData;
%         distance = squareform(pdist(tsneData));        
%
%
%

        cmap = jet(7);
    
%         figure; hold on;
%         for ai = 1 : length(angles)
%             tempinds = find(angleVals == angles(ai));
%             scatter3(Y(tempinds,1), Y(tempinds,2), Y(tempinds,3), 25, cmap(ai,:), 'filled')
%         end
%         title([mouseName, ' ', sessionName])
%         view(-135, 35);

        maxD = max(max(distance)); minD = min(min(distance));
        distance = (distance - minD) ./ (maxD - minD);
        figure;
        imagesc(distance);
        set(gca,'YDir','normal');
        colormap('hot');
        title([mouseName, ' ', sessionName])
        xticks(find(diff(angleVals))+0.5)
        xticklabels({''})
        yticks(find(diff(angleVals))+0.5)
        yticklabels({''})
        set(gca, 'tickdir', 'out', 'ticklength', [0.03 0], 'box', 'off')
    end
end


%% classification performance
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,10,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]};

performance = zeros(length(mice), 2, 3); % (:,1,1) = naive LDA. (:,2,1) = expert LDA. (:,1,2) = naive SVM. (:,2,3) = expert KNN.

for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);

        mouseName = sprintf('JK%03d',mouse);
        sessionName = sprintf('S%02d',session);
        angles = 45:15:135;
        %% Load necessary files
        uberFolder = 'Y:\Whiskernas\JK\suite2p\';
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load([uberFolder, sprintf('%03d',mouse), filesep, ufn])

        %%
    %     it: touch theta, phi, kappaH, kappaV, radialD, touch count
    %     dt: max(dTheta), max(dPhi), max(dKappaH), max(dKappaV), max(touchDuration), max(slideDistance), max(abs(dPhi)), max(abs(dKappaV))
    %     In each trial, these values are averaged.
        touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        answerLickTimeCell = cellfun(@(x) x.answerLickTime, u.trials(touchTrialInds), 'uniformoutput', false);
        for altci = 1 : length(answerLickTimeCell)
            if isempty(answerLickTimeCell{altci})
                answerLickTimeCell{altci} = u.trials{touchTrialInds(altci)}.poleDownOnsetTime;
            end
        end
        
        preDecisionTTI = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y , u.trials(touchTrialInds), answerLickTimeCell));
        preDecisionTouchNum = cellfun(@(x,y) length(cellfun(@(z) find(x.whiskerTime(z(1))<y), x.protractionTouchChunksByWhisking, 'uniformoutput', false)), ...
            u.trials(touchTrialInds(preDecisionTTI)), answerLickTimeCell(preDecisionTTI), 'uniformoutput', false);

        tempData = nan(length(preDecisionTTI),13); % 1:6, it. 7:12, dt. 13, angles
        tempData(:,1) = cellfun(@(x,y) mean(x.theta(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(1:y)))), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,2) = cellfun(@(x,y) mean(x.phi(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(1:y)))), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,3) = cellfun(@(x,y) mean(x.kappaH(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(1:y)))), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,4) = cellfun(@(x,y) mean(x.kappaV(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(1:y)))), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,5) = cellfun(@(x,y) mean(x.arcLength(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(1:y)))), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,6) = cell2mat(preDecisionTouchNum);

        tempData(:,7) = cellfun(@(x,y) mean(x.protractionTouchDThetaByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,8) = cellfun(@(x,y) mean(x.protractionTouchDPhiByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,9) = cellfun(@(x,y) mean(x.protractionTouchDKappaHByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,10) = cellfun(@(x,y) mean(x.protractionTouchDKappaVByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,11) = cellfun(@(x,y) mean(x.protractionTouchDurationByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,12) = cellfun(@(x,y) mean(x.protractionTouchSlideDistanceByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
%         tempData(:,13) = cellfun(@(x,y) mean(x.protractionAbsTouchDPhiByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
%         tempData(:,14) = cellfun(@(x,y) mean(x.protractionAbsTouchDKappaVByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);

        tempData(:,13) = cellfun(@(x) x.angle, u.trials(touchTrialInds(preDecisionTTI)));

        % remove rows with Nans
        rmInds = find(isnan(sum(tempData,2)));
        tempData(rmInds,:) = [];

        %% Re-order data in the order of angles

        [angleVals, I] = sort(tempData(:,end));
        dataAngleOrdered = zeros(size(tempData));
        for ai = 1 : length(angles)
            tempinds = find(angleVals == angles(ai));
            dataAngleOrdered(tempinds,:) = tempData(I(tempinds),:);
        end
        
        [~, performance(mi,si,1)] = angle_tuning_wkv_LDA(tempData, angles);
        [~, performance(mi,si,2)] = angle_tuning_func_reorg_SVM(tempData, angles); % the classification function is the same with that for neuronal activities
        [~, performance(mi,si,3)] = angle_tuning_func_reorg_KNN(tempData, angles); % the classification function is the same with that for neuronal activities
    
    end
end

% %% normality test
for i = 1 : 3
    for j = 1 : 2
        kstest(performance(:,j,i))
    end
end

% %% signed rank test
for i = 1 : 3
    signrank(performance(:,1,i), performance(:,2,i))
end

% %% t-test
for i = 1 : 3
    [~, p] = ttest(performance(:,1,i), performance(:,2,i))
end

% %% making figures
figure,
for i = 1 : 3
    subplot(1,3,i), hold on
    for j = 1 : 6
        plot([1,2], [performance(j,1,i), performance(j,2,i)], 'k-')
    end
    scatter(ones(6,1), performance(:,1,i))
    scatter(ones(6,1)*2, performance(:,2,i))
    xlim([0.5 2.5])
    switch i
        case 1
            title('WKV classfication using LDA')
        case 2
            title('WKV classfication using SVM')
        case 3
            title('WKV classfication using KNN')
    end
end

%%
method = 3;
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.2 0.4]),
hold on
for j = 1 : 6
    plot([1,2], [performance(j,1,method), performance(j,2,method)], 'k-', 'linewidth', 1)
end
scatter(ones(6,1), performance(:,1,method), 'b', 'filled')
scatter(ones(6,1)*2, performance(:,2,method), 'r', 'filled')
xlim([0.5 2.5])
xticks([1,2])
xticklabels({'Naive', 'Expert'})
set(gca, 'linewidth', 1, 'fontsize', 14)
ylabel('Performance')
yticks(0.4:0.1:0.9)



%% Try looking into whisker kinematics changes deeper
% Just some of the key variables: dKappaV, dKappaH, dPhi, slide distance, 
% touch count. In different angles. Compare between before and after.
% From pre-decision kinematics only

mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]};

angles = 45:15:135;
naiveDist = zeros(length(mice), length(angles), 5); % 5 variables
expertDist = zeros(length(mice), length(angles), 5);
for mi = 1 : length(mice)
% for mi = 1
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
%     for si = 1    
        session = sessions{mi}(si);

        mouseName = sprintf('JK%03d',mouse);
        sessionName = sprintf('S%02d',session);
        angles = 45:15:135;

        %% Load necessary files
        uberFolder = 'Y:\Whiskernas\JK\suite2p\';
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load([uberFolder, sprintf('%03d',mouse), filesep, ufn])

        %%
    
        touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        answerLickTimeCell = cellfun(@(x) x.answerLickTime, u.trials(touchTrialInds), 'uniformoutput', false);
        for altci = 1 : length(answerLickTimeCell)
            if isempty(answerLickTimeCell{altci})
                answerLickTimeCell{altci} = u.trials{touchTrialInds(altci)}.poleDownOnsetTime;
            end
        end
        
        preDecisionTTI = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y , u.trials(touchTrialInds), answerLickTimeCell));
        preDecisionTouchNum = cellfun(@(x,y) length(cellfun(@(z) find(x.whiskerTime(z(1))<y), x.protractionTouchChunksByWhisking, 'uniformoutput', false)), ...
            u.trials(touchTrialInds(preDecisionTTI)), answerLickTimeCell(preDecisionTTI), 'uniformoutput', false);
        
        tempData = nan(length(preDecisionTTI),6); % dKappaV, dKappaH, dPhi, slide distance, touch count, and angles
        
        tempData(:,1) = cellfun(@(x,y) mean(x.protractionTouchDKappaVByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,2) = cellfun(@(x,y) mean(x.protractionTouchDKappaHByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,3) = cellfun(@(x,y) mean(x.protractionTouchDPhiByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);        
        tempData(:,4) = cellfun(@(x,y) mean(x.protractionTouchSlideDistanceByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,5) = cellfun(@(x,y) length(x.protractionTouchChunksByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,6) = cellfun(@(x,y) x.angle, u.trials(touchTrialInds(preDecisionTTI)));
        
        % remove rows with Nans
        rmInds = find(isnan(sum(tempData,2)));
        tempData(rmInds,:) = [];
        
        for ai = 1 : length(angles)
            tempInd = find(tempData(:,6) == angles(ai));
            for vi = 1 : 5
                if si == 1
                    naiveDist(mi,ai,vi) = mean(tempData(tempInd,vi));
                else
                    expertDist(mi,ai,vi) = mean(tempData(tempInd,vi));
                end
            end
        end
    end
end
%%
figure,
for vi = 1 : 5
    subplot(2,3,vi), hold on
    tempNaive = squeeze(naiveDist(:,:,vi));
    shadedErrorBar(angles, mean(tempNaive), std(tempNaive)/sqrt(length(mice)), 'lineprops', 'b')
    tempExpert = squeeze(expertDist(:,:,vi));
    shadedErrorBar(angles, mean(tempExpert), std(tempExpert)/sqrt(length(mice)), 'lineprops', 'r')
    xticks(angles)
    xlim([angles(1), angles(end)])
    switch vi
        case 1
            title('Max \Delta\kappa_V')
            ylabel('(mm^-2)')
        case 2
            title('Max \Delta\kappa_H')
            ylabel('(mm^-2)')
        case 3
            title('Max \Delta\phi')
            ylabel('(\circ)')
            xlabel('Object angles (\circ)')
        case 4
            title('Max slide distance')
            xlabel('Object angles (\circ)')
            ylabel('(mm)')
        case 5
            title('Touch count')
            legend({'Naive', 'Expert'}, 'location', 'northeastoutside', 'box', 'off')
            xlabel('Object angles (\circ)')
            ylabel('Number')
    end
end



%% Try looking into whisker kinematics changes deeper
%% Trying just "during touch" features"
% All touch features: dKappaV, dKappaH, dPhi, dTheta, touch duration, 
% slide distance. In different angles. Compare between before and after.
% From pre-decision kinematics only

mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]};

angles = 45:15:135;
naiveDist = zeros(length(mice), length(angles), 6); % 6 variables
expertDist = zeros(length(mice), length(angles), 6);
for mi = 1 : length(mice)
% for mi = 1
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
%     for si = 1    
        session = sessions{mi}(si);

        mouseName = sprintf('JK%03d',mouse);
        sessionName = sprintf('S%02d',session);
        angles = 45:15:135;

        %% Load necessary files
        uberFolder = 'Y:\Whiskernas\JK\suite2p\';
        ufn = sprintf('UberJK%03dS%02d',mouse, session);
        load([uberFolder, sprintf('%03d',mouse), filesep, ufn])

        %%
    
        touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunksByWhisking), u.trials));
        answerLickTimeCell = cellfun(@(x) x.answerLickTime, u.trials(touchTrialInds), 'uniformoutput', false);
        for altci = 1 : length(answerLickTimeCell)
            if isempty(answerLickTimeCell{altci})
                answerLickTimeCell{altci} = u.trials{touchTrialInds(altci)}.poleDownOnsetTime;
            end
        end
        
        preDecisionTTI = find(cellfun(@(x,y) x.whiskerTime(x.protractionTouchChunksByWhisking{1}(1)) < y , u.trials(touchTrialInds), answerLickTimeCell));
        preDecisionTouchNum = cellfun(@(x,y) length(cellfun(@(z) find(x.whiskerTime(z(1))<y), x.protractionTouchChunksByWhisking, 'uniformoutput', false)), ...
            u.trials(touchTrialInds(preDecisionTTI)), answerLickTimeCell(preDecisionTTI), 'uniformoutput', false);
        
        tempData = nan(length(preDecisionTTI),7); % dKappaV, dKappaH, dPhi, slide distance, touch count, and angles
        
        tempData(:,1) = cellfun(@(x,y) mean(x.protractionTouchDThetaByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum); 
        tempData(:,2) = cellfun(@(x,y) mean(x.protractionTouchDPhiByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum); 
        tempData(:,3) = cellfun(@(x,y) mean(x.protractionTouchDKappaHByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,4) = cellfun(@(x,y) mean(x.protractionTouchDKappaVByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,5) = cellfun(@(x,y) mean(x.protractionTouchDurationByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,6) = cellfun(@(x,y) mean(x.protractionTouchSlideDistanceByWhisking(1:y)), u.trials(touchTrialInds(preDecisionTTI)), preDecisionTouchNum);
        tempData(:,7) = cellfun(@(x,y) x.angle, u.trials(touchTrialInds(preDecisionTTI)));
        
        % remove rows with Nans
        rmInds = find(isnan(sum(tempData,2)));
        tempData(rmInds,:) = [];
        
        for ai = 1 : length(angles)
            tempInd = find(tempData(:,7) == angles(ai));
            for vi = 1 : 6
                if si == 1
                    naiveDist(mi,ai,vi) = mean(tempData(tempInd,vi));
                else
                    expertDist(mi,ai,vi) = mean(tempData(tempInd,vi));
                end
            end
        end
    end
end
%
figure,
for vi = 1 : 6
    subplot(2,3,vi), hold on
    tempNaive = squeeze(naiveDist(:,:,vi));
    shadedErrorBar(angles, mean(tempNaive), std(tempNaive)/sqrt(length(mice)), 'lineprops', 'b')
    tempExpert = squeeze(expertDist(:,:,vi));
    shadedErrorBar(angles, mean(tempExpert), std(tempExpert)/sqrt(length(mice)), 'lineprops', 'r')
    xticks(angles)
    xlim([angles(1), angles(end)])
    switch vi
        case 1
            title('Max \Delta\theta')
            ylabel('(\circ)')
        case 2
            title('Max \Delta\phi')
            ylabel('(\circ)')
            xlabel('Object angles (\circ)')
        case 3
            title('Max \Delta\kappa_H')
            ylabel('(mm^-2)')
        case 4
            title('Max \Delta\kappa_V')
            ylabel('(mm^-2)')        
        case 5
            title('Touch Duration')
            xlabel('Object angles (\circ)')
            ylabel('Seconds')            
        case 6
            title('Max slide distance')
            xlabel('Object angles (\circ)')
            ylabel('(mm)')
            legend({'Naive', 'Expert'}, 'location', 'northeast', 'box', 'off')
    end
end