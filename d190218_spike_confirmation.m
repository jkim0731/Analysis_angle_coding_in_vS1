upperTinds = find(cellfun(@(x) ismember(1,x.planes),u.trials));
lowerTinds = find(cellfun(@(x) ismember(5,x.planes),u.trials));

spikes = cell(length(u.trials{upperTinds(1)}.neuindSession) + length(u.trials{lowerTinds(1)}.neuindSession), 1);
calcium = cell(length(u.trials{upperTinds(1)}.neuindSession) + length(u.trials{lowerTinds(1)}.neuindSession), 1);
for i = 1 : length(u.trials{upperTinds(1)}.neuindSession)
    spikes{i} = cell2mat(cellfun(@(x) x.spk(i,:), u.trials(upperTinds)', 'uniformoutput', false));
    calcium{i} = cell2mat(cellfun(@(x) x.dF(i,:), u.trials(upperTinds)', 'uniformoutput', false));
end
for i = 1 : length(u.trials{lowerTinds(1)}.neuindSession)
    spikes{length(u.trials{upperTinds(1)}.neuindSession) + i} = cell2mat(cellfun(@(x) x.spk(i,:), u.trials(lowerTinds)', 'uniformoutput', false));
    calcium{length(u.trials{upperTinds(1)}.neuindSession) + i} = cell2mat(cellfun(@(x) x.dF(i,:), u.trials(lowerTinds)', 'uniformoutput', false));
end
%%
cells = [];
for cind = 1:length(spikes)

    % figure, plot(calcium{cind}, 'k-'), hold on, yyaxis right, plot(spikes{cind}, 'r.'), ylim([0 100])


    ca = calcium{cind};
    spk = spikes{cind};

    space = find(u.trials{1}.tpmTime{1} > 1, 1, 'first'); % interval between compared spikes
    spkTimes = find(spk);

    spkWindows = [];
    for i = 1 : length(spkTimes)
        spkWindows = union(spkWindows, [spkTimes(i)-space:spkTimes(i)+space]);
    end
    spkWindows = spkWindows(spkWindows > 0);

    timediff = diff(spkTimes);
    toBeDeleted = union(find(timediff < space)+1, find(timediff < space));
    spkTimes(toBeDeleted) = [];

    if ~isempty(spkTimes)
        ca(spkWindows) = deal(0);    
        caTimes = find(ca >= min(ca(spkTimes)));

        timediff = diff(caTimes);
        toBeDeleted = union(find(timediff < space)+1, find(timediff < space));
        caTimes(toBeDeleted) = [];    

        if ~isempty(caTimes) && ~isempty(spkTimes)
            cells = [cells, cind];
        end
    end
end
   
%%
numrow = 2;
numcol = 2;

%%
close all
figure('units', 'normalized', 'position', [0.1, 0.1, 0.6, 0.45])
for ci = 1:length(cells)
    cind = cells(ci);

    framesBefore = 3;
    framesAfter = 5;

    ca = calcium{cind};
    spk = spikes{cind};


    % figure, plot(ca, 'k-'), hold on, yyaxis right, plot(find(spk), spk(find(spk)), 'r.'), ylim([1 100])


    space = find(u.trials{1}.tpmTime{1} > 1, 1, 'first'); % interval between compared spikes
    spkTimes = find(spk);

    spkWindows = [];
    for i = 1 : length(spkTimes)
        spkWindows = union(spkWindows, [spkTimes(i)-space:spkTimes(i)+space]);
    end
    spkWindows = spkWindows(spkWindows > 0);

    timediff = diff(spkTimes);
    toBeDeletedSpk = union(find(timediff < space)+1, find(timediff < space));
    spkTimes(toBeDeletedSpk) = [];

    caTemp = ca;
    caTemp(spkWindows) = deal(0);
    caTimes = find(caTemp >= min(ca(spkTimes)));

    timediff = diff(caTimes);
    toBeDeleted = union(find(timediff < space)+1, find(timediff < space));
    caTimes(toBeDeleted) = [];
    spkTimes(spkTimes < framesBefore) = [];
    caTimes(caTimes < framesBefore) = [];

    if min(spk(spkTimes)) == 1
        spkTimes = spkTimes(find(spk(spkTimes) == min(spk(spkTimes))));

        caSpikes = zeros(length(spkTimes), framesBefore + framesAfter + 1);
        caNoSpike = zeros(length(caTimes), framesBefore + framesAfter + 1);

        % figure, hold all
        for i = 1 : length(spkTimes)
            caSpikes(i,:) = ca(spkTimes(i)-framesBefore : spkTimes(i) + framesAfter);
        %     plot(-framesBefore:framesAfter, caSpikes(i,:), 'r-');
        end
        for i = 1 : length(caTimes)
            caNoSpike(i,:) =  ca(caTimes(i)-framesBefore : caTimes(i) + framesAfter);
        %     plot(-framesBefore:framesAfter, caNoSpike(i,:), 'k-');
        end
        %
        subplot(numrow,numcol,ci), 
        errorbar(-framesBefore:framesAfter, mean(caSpikes), std(caSpikes)/sqrt(size(caSpikes,1)), 'r-'), hold on
        errorbar(-framesBefore:framesAfter, mean(caNoSpike), std(caNoSpike)/sqrt(size(caNoSpike,1)), 'k-')
        xlim([-framesBefore framesAfter])
        if mod(ci,numcol) == 1
            ylabel('\DeltaF/F_0')
        end
        if ci > length(cells) - numcol
            xlabel('Frames')
        end
    else
        disp('Min spike > 1')
    end
end

%%
close all
figure('units', 'normalized', 'position', [0.1, 0.1, 0.6, 0.45])
for ci = 1:length(cells)
    cind = cells(ci);
    framesBefore = 3;
    framesAfter = 5;

    ca = calcium{cind};
    spk = spikes{cind};


    % figure, plot(ca, 'k-'), hold on, yyaxis right, plot(find(spk), spk(find(spk)), 'r.'), ylim([1 100])


    space = find(u.trials{1}.tpmTime{1} > 1, 1, 'first'); % interval between compared spikes
    spkTimes = find(spk);

    spkWindows = [];
    for i = 1 : length(spkTimes)
        spkWindows = union(spkWindows, [spkTimes(i)-space:spkTimes(i)+space]);
    end
    spkWindows = spkWindows(spkWindows > 0);

    timediff = diff(spkTimes);
    toBeDeletedSpk = union(find(timediff < space)+1, find(timediff < space));
    spkTimes(toBeDeletedSpk) = [];

    caTemp = ca;
    caTemp(spkWindows) = deal(0);
    caTimes = find(caTemp >= min(ca(spkTimes)));

    timediff = diff(caTimes);
    toBeDeleted = union(find(timediff < space)+1, find(timediff < space));
    caTimes(toBeDeleted) = [];
    spkTimes(spkTimes < framesBefore) = [];
    caTimes(caTimes < framesBefore) = [];

    if min(spk(spkTimes)) == 1
        spkTimes = spkTimes(find(spk(spkTimes) == min(spk(spkTimes))));

        caSpikes = zeros(length(spkTimes), framesBefore + framesAfter + 1);
        caNoSpike = zeros(length(caTimes), framesBefore + framesAfter + 1);

        % figure, hold all
        for i = 1 : length(spkTimes)
            caTemp = ca(spkTimes(i)-framesBefore : spkTimes(i) + framesAfter);
            caSpikes(i,:) = (caTemp - min(caTemp)) / (caTemp(framesBefore+1) - min(caTemp));
        %     plot(-framesBefore:framesAfter, caSpikes(i,:), 'r-');
        end
        for i = 1 : length(caTimes)
            caTemp = ca(caTimes(i)-framesBefore : caTimes(i) + framesAfter);
            caNoSpike(i,:) = (caTemp - min(caTemp)) / (caTemp(framesBefore+1) - min(caTemp));
        %     plot(-framesBefore:framesAfter, caNoSpike(i,:), 'k-');
        end
        %
        subplot(numrow,numcol,ci), 
        errorbar(-framesBefore:framesAfter, mean(caSpikes), std(caSpikes)/sqrt(size(caSpikes,1)), 'r-'), hold on
        errorbar(-framesBefore:framesAfter, mean(caNoSpike), std(caNoSpike)/sqrt(size(caNoSpike,1)), 'k-')
        xlim([-framesBefore framesAfter])
        if mod(ci,numcol) == 1
            ylabel('Normalized \DeltaF/F_0')
        end
        if ci > length(cells) - numcol
            xlabel('Frames')
        end
    else
        disp('Min spike > 1')
    end
end


%%
baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]}; 
totalinds = sum(cellfun(@(x) length(x), sessions));
wholeCorrect = zeros(totalinds,1);
pTunedSpk = zeros(totalinds,1);
pTunedCa = zeros(totalinds,1);
pNotuneSpk = zeros(totalinds,1);
pNotuneCa = zeros(totalinds,1);
pTunedGivenCaTuned = zeros(totalinds,1);
pTunedGivenSpkTuned = zeros(totalinds,1);
pNotuneGivenCaNotune = zeros(totalinds,1);
pNotuneGivenSpkNotune = zeros(totalinds,1);
pNoresponseGivenCaNoresponse = zeros(totalinds,1);
pNoresponseGivenSpkNoresponse = zeros(totalinds,1);
pTouchGivenCaTouch = zeros(totalinds,1);
pTouchGivenSpkTouch = zeros(totalinds,1);
pNoresponseGivenSpkNotune = zeros(totalinds,1);
pTunedGivenSpkNotune = zeros(totalinds,1);
pNotuneGivenSpkTuned = zeros(totalinds,1);
pNoresponseGivenSpkTuned = zeros(totalinds,1);
pTunedGivenSpkNoresponse = zeros(totalinds,1);
pNotuneGivenSpkNoresponse = zeros(totalinds,1);

currenti = 0;
for mi = 1 : length(mice)
    mouse = mice(mi);
    cd(sprintf('%s%03d', baseDir, mouse))
    for si = 1 : length(sessions{mi})        
        currenti = currenti + 1;
        
        session = sessions{mi}(si);
        
        ca = load(sprintf('JK%03dS%02dsingleCell_anova_calcium_final',mouse, session), 'cellsTuned', 'cellsNTResponse', 'cellsNTNR', 'tuneAngle', 'tuneSingle', 'cellsTotal');
        spk = load(sprintf('JK%03dS%02dsingleCell_anova_spk_final',mouse, session), 'cellsTuned', 'cellsNTResponse', 'cellsNTNR', 'tuneAngle', 'tuneSingle', 'cellsTotal');

        wholeCorrect(currenti) = (length(intersect(ca.cellsTuned, spk.cellsTuned)) + length(intersect(ca.cellsNTResponse, spk.cellsNTResponse)) + length(intersect(ca.cellsNTNR, spk.cellsNTNR))) / length(ca.cellsTotal);

        pTunedSpk(currenti) = length(spk.cellsTuned)/length(spk.cellsTotal);
        pTunedCa(currenti) = length(ca.cellsTuned)/length(ca.cellsTotal);
        pNotuneSpk(currenti) = length(spk.cellsNTResponse) / length(spk.cellsTotal);
        pNotuneCa(currenti) = length(ca.cellsNTResponse) / length(ca.cellsTotal);
        
        pTunedGivenCaTuned(currenti) = length(intersect(spk.cellsTuned, ca.cellsTuned)) / length(ca.cellsTuned);
        pTunedGivenSpkTuned(currenti) = length(intersect(spk.cellsTuned, ca.cellsTuned)) / length(spk.cellsTuned);

        pNotuneGivenCaNotune(currenti) = length(intersect(spk.cellsNTResponse, ca.cellsNTResponse)) / length(spk.cellsNTResponse);
        pNotuneGivenSpkNotune(currenti) = length(intersect(spk.cellsNTResponse, ca.cellsNTResponse)) / length(ca.cellsNTResponse);

        pNoresponseGivenCaNoresponse(currenti) = length(intersect(spk.cellsNTNR, ca.cellsNTNR)) / length(spk.cellsNTNR);
        pNoresponseGivenSpkNoresponse(currenti) = length(intersect(spk.cellsNTNR, ca.cellsNTNR)) / length(ca.cellsNTNR);

        % touch response classified as touch cells of ca
        pTouchGivenCaTouch(currenti) = length(intersect(union(ca.cellsTuned, ca.cellsNTResponse), union(spk.cellsTuned, spk.cellsNTResponse))) / length(union(ca.cellsTuned, ca.cellsNTResponse));
        % touch response classified as touch cells of spk
        pTouchGivenSpkTouch(currenti) = length(intersect(union(ca.cellsTuned, ca.cellsNTResponse), union(spk.cellsTuned, spk.cellsNTResponse))) / length(union(spk.cellsTuned, spk.cellsNTResponse));
        
        % how much of mismatch in no-tune in ca classified as no-response in spk?
        pNoresponseGivenSpkNotune(currenti) = length(intersect(setdiff(ca.cellsNTResponse, spk.cellsNTResponse), spk.cellsNTNR) ) / length(setdiff(ca.cellsNTResponse, spk.cellsNTResponse));
        % how much of mismatch in no-tune in ca classified as tuned in spk?
        pTunedGivenSpkNotune(currenti) = length(intersect(setdiff(ca.cellsNTResponse, spk.cellsNTResponse), spk.cellsTuned) ) / length(setdiff(ca.cellsNTResponse, spk.cellsNTResponse));

        % how much of mismatch in tuned in ca classified as no-tune in spk?
        pNotuneGivenSpkTuned(currenti) = length(intersect(setdiff(ca.cellsTuned, spk.cellsTuned), spk.cellsNTResponse) ) / length(setdiff(ca.cellsTuned, spk.cellsTuned));
        % how much of mismatch in tuned in ca classified as no-response in spk?
        pNoresponseGivenSpkTuned(currenti) = length(intersect(setdiff(ca.cellsTuned, spk.cellsTuned), spk.cellsNTNR) ) / length(setdiff(ca.cellsTuned, spk.cellsTuned));

        % how much of mismatch in no-response in ca classified as tuned in spk?
        pTunedGivenSpkNoresponse(currenti) = length(intersect(setdiff(ca.cellsNTNR, spk.cellsNTNR), spk.cellsTuned) ) / length(setdiff(ca.cellsNTNR, spk.cellsNTNR));
        % how much of mismatch in no-response in ca classified as no-tune in spk?
        pNotuneGivenSpkNoresponse(currenti) = length(intersect(setdiff(ca.cellsNTNR, spk.cellsNTNR), spk.cellsNTResponse) ) / length(setdiff(ca.cellsNTNR, spk.cellsNTNR));
    end
end

%%
baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]}; 
totalinds = sum(cellfun(@(x) length(x), sessions));

sameTuneAngle = zeros(totalinds,1);
diffTuneDirectionNegative = zeros(totalinds,1);
mismatchTuneDirectionNegativeCa = zeros(totalinds,1);
pNegativeTunedCa = zeros(totalinds,1);
pNegativeTunedSpk = zeros(totalinds,1);

currenti = 0;
for mi = 1 : length(mice)
% for mi = 1    
    mouse = mice(mi);
    cd(sprintf('%s%03d', baseDir, mouse))
    for si = 1 : length(sessions{mi})
%     for si = 1
        currenti = currenti + 1;
        
        session = sessions{mi}(si);
        %%
        ca = load(sprintf('JK%03dS%02dsingleCell_anova_calcium_final',mouse, session), 'cellsTuned', 'tuneAngle', 'tuneSingle', 'tuneDirection', 'cellsTotal');
        spk = load(sprintf('JK%03dS%02dsingleCell_anova_spk_final',mouse, session), 'cellsTuned', 'tuneAngle', 'tuneSingle', 'tuneDirection', 'cellsTotal');
        
        [caind, spkind] = ismember(ca.cellsTuned, spk.cellsTuned);
        
        
        sameinds = find(ca.tuneAngle(caind) == spk.tuneAngle(spkind(find(spkind))));
        
        sameTuneAngle(currenti) = length(sameinds)/sum(caind);
        diffinds = setdiff(find(caind), sameinds);
        diffTuneDirectionNegative(currenti) = length(find(ca.tuneDirection(diffinds) > 1)) / length(diffinds);
        mismatchTuneDirectionNegativeCa(currenti) = length(find(ca.tuneDirection(find(caind == 0)) > 1)) / length(find(caind == 0));
        
        pNegativeTunedCa(currenti) = length(find(ca.tuneDirection > 1)) / length(ca.tuneDirection);
        pNegativeTunedSpk(currenti) = length(find(spk.tuneDirection > 1)) / length(spk.tuneDirection);        
    end
end
        
        