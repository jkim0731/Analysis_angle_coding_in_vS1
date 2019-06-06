% data structure to make:
% Whisker kinematic variables. all & divided into planes. Alongside with angles. Better if
% matched with 'vals' of 'angle_tuning_summary'
% averaged in each trial. Multiple touches will have averaged "at touch"
% variables and "during touch" variables

% Dependencies: C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\code\dependentFunctions

    mouse = 25;
    session = 4;

    mouseName = sprintf('JK%03d',mouse);
    sessionName = sprintf('S%02d',session);
    
    %% Load necessary files
    behaviorFolder = 'Y:\Whiskernas\JK\SoloData\';
    load([behaviorFoler, mouseName, filesep, 'behavior_' mouseName '.mat'])
    
    uberFolder = 'Y:\Whiskernas\JK\suite2p\';        
    ufn = sprintf('UberJK%03dS%02d',mouse, session);
    load([uberFolder, sprintf('%03d',mouse), filesep, ufn])
    
    whiskerFolder = ['Y:\Whiskernas\JK\whisker\tracked' filesep mouseName sessionName];
    
    bSessionNums = cellfun(@(x) x.sessionName, b,'uniformoutput',false);

    %find behavioral data matching session and load bMat
    bMatIdx = find(cell2mat(cellfun(@(x) strcmp(x,sessionNumber),bSessionNums,'uniformoutput',false)));
    behavioralStruct = b{bMatIdx};
    wfa = Whisker.WhiskerFinal_2padArray(whiskerFolder);
    
    %% Touch feature builder
    %instantaneous touch features
    it = instantTouchBuilder(behavioralStruct,wfa,'protraction','all');
    %during touch features
    dt = duringTouchBuilder(behavioralStruct,wfa,'protraction','all');
    
    %%
%     it: touch theta, phi, kappaH, kappaV, radialD, touch count
%     dt: max(dTheta), max(dPhi), max(dKappaH), max(dKappaV), max(touchDuration), max(slideDistance), max(abs(dPhi)), max(abs(dKappaV))
    
    touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTFchunks), wfa.trials));
    wkv.theta = cellfun(@(x) mean(x.theta(cellfun(@(y) y(1), x.protractionTFchunks))), wfa.trials(touchTrialInds));
    wkv.phi = cellfun(@(x) mean(x.phi(cellfun(@(y) y(1), x.protractionTFchunks))), wfa.trials(touchTrialInds));
    wkv.kappaH = cellfun(@(x) mean(x.kappaH(cellfun(@(y) y(1), x.protractionTFchunks))), wfa.trials(touchTrialInds));
    wkv.kappaV = cellfun(@(x) mean(x.kappaV(cellfun(@(y) y(1), x.protractionTFchunks))), wfa.trials(touchTrialInds));
    wkv.radialD = cellfun(@(x) mean(x.arcLength(cellfun(@(y) y(1), x.protractionTFchunks))), wfa.trials(touchTrialInds));
    wkv.touchCount = cellfun(@(x) length(x.protractionTFchunks), wfa.trials(touchTrialInds));
    wkv.maxdTheta = 