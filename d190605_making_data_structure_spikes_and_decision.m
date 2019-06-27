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

% Divide into naive, expert, and nonlearner (naive and expert order matched)

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

naiveInd = find(cellfun(@length, sessions)==2);
nonlearnerInd = find(cellfun(@length, sessions)==1);

naive = struct;
expert = struct;
nonlearner = struct;
for ni = 1 : length(naiveInd)
    mi = naiveInd(ni);
    mouse = mice(mi);
    session = sessions{mi}(1);
    results = make_data_structure_func(mouse, session);
    fields = fieldnames(results);
    for fi = 1 : length(fields)
        naive(ni).(fields{fi}) = results.(fields{fi});
    end
    
    session = sessions{mi}(2);
    results = make_data_structure_func(mouse, session);
    fields = fieldnames(results);
    for fi = 1 : length(fields)
        expert(ni).(fields{fi}) = results.(fields{fi});
    end
end

for ni = 1 : length(nonlearnerInd)
    mi = nonlearnerInd(ni);
    mouse = mice(mi);
    session = sessions{mi};
    results = make_data_structure_func(mouse, session);
    fields = fieldnames(results);
    for fi = 1 : length(fields)
        nonlearner(ni).(fields{fi}) = results.(fields{fi});
    end
end

savefn = 'data_decisionSpk.mat';
save([baseDir, savefn], 'naive', 'expert', 'nonlearner')

function results = make_data_structure_func(mouse, session)
    %% Load uber file
    uberFolder = 'D:\TPM\JK\suite2p\';        
    ufn = sprintf('UberJK%03dS%02d',mouse, session);
    load([uberFolder, sprintf('%03d',mouse), filesep, ufn], 'u')
    
    %%
    choices = zeros(length(u.trials),1); % 1: right, 2: left, 0: miss    
    choices(find(cellfun(@(x) strcmp(x.choice, 'r'), u.trials))) = deal(1);
    choices(find(cellfun(@(x) strcmp(x.choice, 'l'), u.trials))) = deal(2);
    responses = cellfun(@(x) x.response, u.trials); % 1 hit, 0 wrong, -1 miss
    angles = cellfun(@(x) x.angle, u.trials);
%     touchTrialInds = find(cellfun(@(x) ~isempty(x.protractionTouchChunks), u.trials));
%     decisionTrialInds = find(cellfun(@(x) ~isempty(x.answerLickTime), u.trials));
    upperLayerInds = find(cellfun(@(x) ismember(1,x.planes), u.trials));
    lowerLayerInds = find(cellfun(@(x) ismember(5,x.planes), u.trials));    
    
%     testInds = intersect(touchTrialInds, decisionTrialInds);
%     testTrials{1} = u.trials(intersect(testInds, upperLayerInds));
%     testTrials{2} = u.trials(intersect(testInds, lowerLayerInds));
%     testChoices{1} = choices(intersect(testInds, upperLayerInds));
%     testChoices{2} = choices(intersect(testInds, lowerLayerInds));
%     testAngles{1} = angles(intersect(testInds, upperLayerInds));
%     testAngles{2} = angles(intersect(testInds, lowerLayerInds));
    testTrials{1} = u.trials(upperLayerInds);
    testTrials{2} = u.trials(lowerLayerInds);
    testChoices{1} = choices(upperLayerInds);
    testChoices{2} = choices(lowerLayerInds);
    testAngles{1} = angles(upperLayerInds);
    testAngles{2} = angles(lowerLayerInds);
    outcome{1} = responses(upperLayerInds);
    outcome{2} = responses(lowerLayerInds);
    
    decisionFrames = cell(8,1);
    poleUpFrames = cell(8,1);
    for i = 1 : 2
        for j = 1 : 4
            decisionFrames{(i-1)*4+j} = cell(length(testTrials{i}),1);
            poleUpFrames{(i-1)*4+j} = cell(length(testTrials{i}),1);
            decisionTind = find(cellfun(@(x) ~isempty(x.answerLickTime), testTrials{i}));
            nodecisionTind = find(cellfun(@(x) isempty(x.answerLickTime), testTrials{i}));
            decisionFrames{(i-1)*4 + j}(decisionTind) = cellfun(@(x) find(x.tpmTime{j} > x.answerLickTime,1), testTrials{i}(decisionTind), 'uniformoutput', false);
            decisionFrames{(i-1)*4 + j}(nodecisionTind) = cellfun(@(x) find(x.tpmTime{j} > x.poleDownOnsetTime,1), testTrials{i}(nodecisionTind), 'uniformoutput', false);
            poleUpFrames{(i-1)*4 + j} = cellfun(@(x) find(x.tpmTime{j} > x.poleUpOnsetTime,1), testTrials{i}, 'uniformoutput', false);
        end
    end
    
    cellGroup{1} = u.cellNums(u.cellNums < 5000);
    cellGroup{2} = u.cellNums(u.cellNums > 5000);
    spkVal = cell(length(u.cellNums),1);
    for ci = 1 : length(u.cellNums)
        cellID = u.cellNums(ci);
        if cellID < 5000
            currTrials = testTrials{1};
        else
            currTrials = testTrials{2};
        end
        currPlaneNum = floor(cellID/1000);
        spkInd = find(currTrials{1}.neuindSession == cellID);
        spkVal{ci} = cellfun(@(x,y,z) mean(x.spk(spkInd,y:z)), currTrials, poleUpFrames{currPlaneNum}, decisionFrames{currPlaneNum});
    end
    
    %%
    dataTable = cell(1,2);
    for i = 1 : 2
        dataTable{i} = cell2mat(spkVal(find(cellfun(@(x) length(x) == length(testChoices{i}), spkVal)))');
    end
    
    %% load angle tuning and glm results
    atfn = sprintf('JK%03dS%02dangle_tuning',mouse, session);
    load([uberFolder, sprintf('%03d',mouse), filesep, atfn], 'spk')
    touchInd{1} = find(ismember(u.cellNums(find(u.cellNums < 5000)), spk.touchID));
    touchInd{2} = find(ismember(u.cellNums(find(u.cellNums > 5000)), spk.touchID));
    tunedInd{1} = find(ismember(u.cellNums(find(u.cellNums < 5000)), spk.touchID(find(spk.tuned))));
    tunedInd{2} = find(ismember(u.cellNums(find(u.cellNums > 5000)), spk.touchID(find(spk.tuned))));
    
    glmfn = sprintf('JK%03dS%02dglm_cell_function',mouse, session);
    load([uberFolder, sprintf('%03d',mouse), filesep, glmfn], 'glm')
    whiskingInd{1} = find(ismember(u.cellNums(find(u.cellNums < 5000)), glm.whiskingID));
    whiskingInd{2} = find(ismember(u.cellNums(find(u.cellNums > 5000)), glm.whiskingID));
    
    results.spk = dataTable;
    results.choice = testChoices;
    results.angle = testAngles;
    results.outcome = outcome;
    results.touchInd = touchInd;
    results.tunedInd = tunedInd;
    results.whiskingInd = whiskingInd;
    results.cellID = cellGroup;
end
    