clear all; close all; clc;

load('Y:/Whiskernas/JK/suite2p/angle_tuning_summary');
processedResults = struct();

%% params
dataSet = expert; 
% animalID = 4; % which animal in the dataset to analyse
angles = linspace(45, 135, 7); % how the dataset angle IDs should be classed
layers = [0 5000]; % choose layer ID limits
angleDiscrim = [45; 60; 75; 90; 105; 120; 135]; % which classes should classifier separate?
reps = 10; % repetitions for classifier training

% for animalID = 1:length(dataSet)
for animalID = 1
    disp(['Starting analysis for animal ' num2str(animalID)]);
    
    %% reshape data
    reshapeData;

    %% t-sne
    runTsne;

    %% run clasifier
    runClassifier;

    %% weight analysis
    weightAnalysis;

    %% save out
    processedResults(animalID).accuracy = accuracy;
    processedResults(animalID).accuracyShuffled = accuracyShuffled;
    processedResults(animalID).allCellChoice = allCellChoice;
%     processedResults(animalID).allClassifiers = allClassifiers;
    processedResults(animalID).angleDistance = angleDistance;
    processedResults(animalID).angleDiscrim = angleDiscrim;
    processedResults(animalID).angles = angles;
    processedResults(animalID).angleVals = angleVals;
    processedResults(animalID).bestID = bestID;
    processedResults(animalID).cmap = cmap;
    processedResults(animalID).distance = distance;
    processedResults(animalID).layers = layers;
    processedResults(animalID).nCells = nCells;
    processedResults(animalID).properties = properties;
    processedResults(animalID).tableData = tableData;
    processedResults(animalID).worstID = worstID;
    processedResults(animalID).reps = reps;
    processedResults(animalID).Y = Y;
end

