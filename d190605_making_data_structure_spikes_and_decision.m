% data structure to make:
% Spikes and decision. 
% (1) Spikes from all cells. (spikes (trial num) X (duration))
%     (1)-1 Whether a cell is touch cell, tuned cell, whisking cell.
%     (1)-2 ref cell ID (only from learners; for before learning session and distance coding session)
%     (1)-3 which trials the cell was imaged
%     (1)-4 Touch frames
% (2) Decision
%     (2)-1 Lick patterns in each trial
%     (2)-2 Decision lick
%     (2)-3 outcome (hit, miss, wrong)

% Divide into naive, expert, and nonlearner (naive and expert order
% matched)

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
    


    %% Choice and ttype builder
    outcomes = BMatBuilder(behavioralStruct,wfa);
    outcomes.mouseName = mouseName;
    outcomes.sessionName = sessionName;