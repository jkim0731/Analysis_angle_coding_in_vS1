% clear all; close all; clc;

load('Y:/Whiskernas/JK/suite2p/angle_tuning_summary');
processedResults = struct();

%% params
dataSet = expert;
% dataSet = naive;
% dataSet = expert([1,3:6]);
% dataSet = naive([1,3,4,7,9]);
miceName = {'JK025','JK027','JK030','JK036','JK039','JK052'};
% animalID = 4; % which animal in the dataset to analyse
angles = linspace(45, 135, 7); % how the dataset angle IDs should be classed
% layers = [0 5000]; % choose layer ID limits
layers = [5000 9000]; % choose layer ID limits
% depth = [350 600]; % L4
depth = [136 350]; % L3
% depth = [1 136]; % L2
% depth = [1 350]; % L2/3

angleDiscrim = [45; 60; 75; 90; 105; 120; 135]; % which classes should classifier separate?
reps = 10; % repetitions for classifier training

classificationPerformance = zeros(length(dataSet),3);

% for animalID = 1:length(dataSet)
for animalID = 3:4
    disp(['Starting analysis for animal ' num2str(animalID)]);
    
    %% reshape data
    reshapeData;

    %% t-sne
    runTsne;

    %% classification
%     [~, classificationPerformance(animalID,1)] = angle_tuning_func_reorg_LDA(classificationSetForParameterSetting, angles);
%     [~, classificationPerformance(animalID,2)] = angle_tuning_func_reorg_SVM(classificationSetForParameterSetting, angles);
%     [~, classificationPerformance(animalID,3)] = angle_tuning_func_reorg_KNN(classificationSetForParameterSetting, angles);
%     %% run clasifier
%     runClassifier;
% 
%     %% weight analysis
%     weightAnalysis;
% 
%     %% save out
%     processedResults(animalID).accuracy = accuracy;
%     processedResults(animalID).accuracyShuffled = accuracyShuffled;
%     processedResults(animalID).allCellChoice = allCellChoice;
% %     processedResults(animalID).allClassifiers = allClassifiers;
%     processedResults(animalID).angleDistance = angleDistance;
%     processedResults(animalID).angleDiscrim = angleDiscrim;
%     processedResults(animalID).angles = angles;
%     processedResults(animalID).angleVals = angleVals;
%     processedResults(animalID).bestID = bestID;
%     processedResults(animalID).cmap = cmap;
%     processedResults(animalID).distance = distance;
%     processedResults(animalID).layers = layers;
%     processedResults(animalID).nCells = nCells;
%     processedResults(animalID).properties = properties;
%     processedResults(animalID).tableData = tableData;
%     processedResults(animalID).worstID = worstID;
%     processedResults(animalID).reps = reps;
%     processedResults(animalID).Y = Y;
end

