load('Y:/Whiskernas/JK/suite2p/angle_tuning_summary');
processedResults = struct();

%% params
dataSet{1} = nonlearner(1:6);
dataSet{2} = naive;
dataSet{3} = expert;
% dataSet = expert([1,3:6]);
% dataSet = naive([1,3,4,7,9]);
% miceName = {'JK025','JK027','JK030','JK036','JK039','JK052'};
% animalID = 4; % which animal in the dataset to analyse
angles = linspace(45, 135, 7); % how the dataset angle IDs should be classed
% layers = [0 5000]; % choose layer ID limits
layers = [5000 9000]; % choose layer ID limits
% depth = [350 600]; % L4
% depth = [136 350]; % L3
% depth = [1 136]; % L2
% depth = [1 350]; % L2/3
depth = [136 600]; % L3&4
angleDiscrim = [45; 60; 75; 90; 105; 120; 135]; % which classes should classifier separate?
reps = 10; % repetitions for classifier training

classificationPerformance = cell(3,1);
classificationPerformance{1} = nan(length(dataSet{1}),3);
classificationPerformance{2} = nan(length(dataSet{2}),3);
classificationPerformance{3} = nan(length(dataSet{3}),3);

for i = 1 : 3
    currDataSet = dataSet{i};
    for animalID = 1:length(currDataSet)
    % for animalID = 2
        disp(['Starting analysis for animal ' num2str(animalID)]);

        %% reshape data
        tableData = [];
        angleVals = [];
        for roiID = 1:length(currDataSet(animalID).touchID)
            roiName = currDataSet(animalID).touchID(roiID);
            roiDepth = currDataSet(animalID).depth(roiID);
            if (roiName > layers(1) && roiName < layers(2))
                if (roiDepth >= depth(1) && roiDepth < depth(2))
                    responses = currDataSet(animalID).val{roiID};

                    angleVals = [];
                    responseVals = [];
                    for a = 1:length(responses)
                        angleVals = [angleVals; repmat(angles(a), length(responses{a}), 1)];
                        responseVals = [responseVals; responses{a}];
                    end
                    tableData = [tableData, responseVals];
                end
            end
        end
        classificationSetForParameterSetting = [tableData, angleVals];

        %% t-sne
    %     runTsne;

        %% classification
        if ~isempty(classificationSetForParameterSetting)
            [~, classificationPerformance{i}(animalID,1)] = angle_tuning_func_reorg_LDA(classificationSetForParameterSetting, angles);
            [~, classificationPerformance{i}(animalID,2)] = angle_tuning_func_reorg_SVM(classificationSetForParameterSetting, angles);
            [~, classificationPerformance{i}(animalID,3)] = angle_tuning_func_reorg_KNN(classificationSetForParameterSetting, angles);
        end
    end
end

save('AngleFuncReorg_L3and4', 'classificationPerformance')

%%

% load('AngleFuncReorg_L23')
% load('AngleFuncReorg_L2')
% load('AngleFuncReorg_L3UpperLayer')
load('AngleFuncReorg_L3LowerLayer')
% load('AngleFuncReorg_L4')
% load('AngleFuncReorg_L3and4')
methodNames = {'LDA', 'SVM', 'KNN'};
colors = {'k', 'b', 'r'};
%
method = 3; % 1 LDA, 2 SVM, 3 KNN
data = [];
for i = 1 : 3
    data = [data, classificationPerformance{i}(:,method)];
end
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.2 0.4])
hold on
for i = 1 : 6
    plot([2 3], [data(:,2), data(:,3)], 'k-', 'linewidth', 1)
end
for i = 1 : 3
%     if i == 2
%         scatter(ones(1,5)*i, data([1,3:6],i), colors{i}, 'filled')
%     else
        scatter(ones(1,6)*i, data(:,i), colors{i}, 'filled')
%     end
end
xlim([0.5 3.5])
xticks(1:3)
yticks(0.3:0.1:0.8)
xticklabels({'Nonlearner', 'Naive', 'Expert'})
xtickangle(45)
set(gca, 'linewidth', 1, 'fontsize', 14)
ylabel('Performance')
%%
[~, p] = ttest(data(:,3)-data(:,2))


%%
clear

load('Y:/Whiskernas/JK/suite2p/angle_tuning_summary');
processedResults = struct();


dataSet{1} = nonlearner(1:6);
dataSet{2} = naive;
dataSet{3} = expert;


angles = linspace(45, 135, 7); % how the dataset angle IDs should be classed
layers = [0 5000]; % choose layer ID limits
% layers = [5000 9000]; % choose layer ID limits
% depth = [350 600]; % L4
% depth = [136 350]; % L3
% depth = [1 136]; % L2
depth = [1 350]; % L2/3
% depth = [136 600]; % L3&4
angleDiscrim = [45; 60; 75; 90; 105; 120; 135]; % which classes should classifier separate?
angleDiffs = 0:15:90;
tsneDistances = cell(3,1);

for i = 1 : 3

    tsneDistances{i} = cell(length(angleDiffs),1);
    for j = 1 : length(tsneDistances{i})
        tsneDistances{i}{j} = [];
    end
    
    currDataSet = dataSet{i};
    for animalID = 1:length(currDataSet)    
        disp(['Starting analysis for animal ' num2str(animalID)]);

        % reshape data
        tableData = [];
        angleVals = [];
        for roiID = 1:length(currDataSet(animalID).touchID)
            roiName = currDataSet(animalID).touchID(roiID);
            roiDepth = currDataSet(animalID).depth(roiID);
            if (roiName > layers(1) && roiName < layers(2))
                if (roiDepth >= depth(1) && roiDepth < depth(2))
                    responses = currDataSet(animalID).val{roiID};

                    angleVals = [];
                    responseVals = [];
                    for a = 1:length(responses)
                        angleVals = [angleVals; repmat(angles(a), length(responses{a}), 1)];
                        responseVals = [responseVals; responses{a}];
                    end
                    tableData = [tableData, responseVals];
                end
            end
        end
        
        if size(tableData,2)>10
%             Y = tsne(tableData,'NumDimensions',3);
            Y = tableData;
            % distance t-sne space (change Y to tableData to look in response space)
            distance = squareform(pdist(Y));
            % distance = squareform(pdist(tableData));
            maxD = max(max(distance)); minD = min(min(distance));
            distance = (distance - minD) ./ (maxD - minD);
            angleDistance = squareform(pdist(angleVals));
            tempInds = triu(ones(size(distance,1)),1);
            tsneDist = distance(find(tempInds(:)));
            angleDist = angleDistance(find(tempInds(:)));
            for j = 1 : length(angleDiffs)
                tempInd = find(angleDist == angleDiffs(j));
%                 tsneDistances{i}{j} = [tsneDistances{i}{j}; tsneDist(tempInd)];
                tsneDistances{i}{j} = [tsneDistances{i}{j}; mean(tsneDist(tempInd))];
            end
        end
    end
end


%

save('populationDistanceGroup_L23.mat', 'tsneDistances')
% save('populationTsneDistanceGroup_L23.mat', 'tsneDistances')
%%

% load('populationTsneDistanceGroup_L3LL.mat', 'tsneDistances')
load('populationDistanceGroup_L3LL.mat', 'tsneDistances')
%%
colors = {'k', 'b', 'r'};
figure, hold on
for i = 1 : 3
    shadedErrorBar(angleDiffs, cellfun(@mean, tsneDistances{i}), cellfun(@(x) std(x) / sqrt(length(x)), tsneDistances{i}), 'lineprops',colors{i})
%     shadedErrorBar(angleDiffs, cellfun(@mean, tsneDistances{i}), cellfun(@std, tsneDistances{i}), 'lineprops',colors{i})
end
xticks(angleDiffs)
xlabel('Angle difference (\circ)')
ylabel('Normalized t-SNE distance')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off')
set(gca, 'fontsize', 14)
