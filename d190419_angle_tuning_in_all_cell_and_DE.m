baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
results = load('glm_results_responseType');
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]};
%% 
i = 1;
mouse = mice(i);
cd(sprintf('%s%03d',baseDir,mouse))
session = sessions{i}(1);
tuningAll = load(sprintf('JK%03dS%02dangle_tuning_all_cell',mouse,session));

range = -0.1:0.02:0.6;

%%
ind = find(tuningAll.spk.tuned);
figure
subplot(121), histogram(results.naive(i).allDE, range)
hold on
histogram(results.naive(i).allDE(ind), range)
xlabel('Deviance Explained')
ylabel('# of cells')
legend({'All cells', 'Angle-tuned'})
subplot(122), plot(results.naive(i).allDE, tuningAll.spk.sharpness, 'k.')
xlabel('Deviance Explained')
ylabel('tuning sharpness')

%% what about modulation?
range = -0.1:0.02:0.6;
ind = find(tuningAll.spk.tuned);
subplot(121), histogram(results.naive(i).allDE, range)
hold on
histogram(results.naive(i).allDE(ind), range)
xlabel('Deviance Explained')
ylabel('# of cells')
legend({'All cells', 'Angle-tuned'})
subplot(122), plot(results.naive(i).allDE, tuningAll.spk.modulation, 'k.')
xlabel('Deviance Explained')
ylabel('tuning sharpness')

%% what is the variability in fitting?
repeat = 10;
eachDE = zeros(length(results.naive(i).allDE), repeat);
for ri = 1 : repeat
    load(sprintf('glmResponseType_JK%03dS%02d_m45_R%02d',mouse, session, ri), 'fitDevExplained')
    eachDE(:,ri) = fitDevExplained;
end

fitVar = std(eachDE,[],2);

figure, scatter(results.naive(i).allDE, tuningAll.spk.sharpness, fitVar*500)


%% is low DE because of sparsity?
lowDEtunedInd = intersect(find(results.naive(i).allDE < 0.1), find(tuningAll.spk.tuned));
highDEtunedInd = intersect(find(results.naive(i).allDE >= 0.1), find(tuningAll.spk.tuned));


range = 0:0.05:1;




lowAllSpkFrames = cellfun(@(x) length(find(cell2mat(x)))/length(cell2mat(x)), tuningAll.spk.val(lowDEtunedInd));

% lowTunedSpk = cellfun(@(x) length(x{})/length(x), tuningAll.spk.val(lowDEtunedInd));
highAllSpkFrames = cellfun(@(x) length(find(cell2mat(x)))/length(cell2mat(x)), tuningAll.spk.val(highDEtunedInd));
% highTunedSpk 
figure, histogram(lowAllSpkFrames, range, 'normalization', 'probability'), hold on, histogram(highAllSpkFrames, range, 'normalization', 'probability')
xlabel('Proportion of active trials')
ylabel('Proportion')
legend({'Dev Exp < 0.1', 'Dev Exp > 0.1'})


%%
range = -0.2:0.05:1;
lowAllSpk = cellfun(@(x) sum(cell2mat(x))/length(cell2mat(x)), tuningAll.spk.val(lowDEtunedInd));
highAllSpk = cellfun(@(x) sum(cell2mat(x))/length(cell2mat(x)), tuningAll.spk.val(highDEtunedInd));
figure, histogram(lowAllSpk, range, 'normalization', 'probability'), hold on, histogram(highAllSpk, range, 'normalization', 'probability')
xlabel('Sum of spikes in each trial / # of trials')
ylabel('Proportion')
legend({'Dev Exp < 0.1', 'Dev Exp > 0.1'})

%% Where is the negative sharpness coming from? highly whisking-tuned?
cd(baseDir)
fit = load('cellFunctionRidgeDE010');
whiskingID = fit.naive(i).whiskingID;
negSharpID = tuningAll.spk.touchID(find(tuningAll.spk.sharpness<0));
length(intersect(whiskingID, negSharpID)) / length(negSharpID)
%% 
%% -> No. less than 10%
%%

%% let's look at individual examples
cd(sprintf('%s%03d',baseDir,mouse))
load(sprintf('UberJK%03dS%02d',mouse,session))
example_angle_tuning_calcium(u, tuningAll.ca, tuningAll.spk, negSharpID)

%%
%% -> Most of negative tuning cells look like having significantly different baseline fluorescence or spike activity
%% -> Baseline ANOVA can be used for a first pass for angle tuning
%%

%% What about low DE positively tuned cells? How do they look like?
lowDEtunedInd = intersect(find(results.naive(i).allDE < 0.1), find(tuningAll.spk.tuned));
highDEtunedInd = intersect(find(results.naive(i).allDE >= 0.1), find(tuningAll.spk.tuned));
negSharpID = tuningAll.spk.touchID(find(tuningAll.spk.sharpness<0));
lowDEtunedID = results.naive(i).cellID(lowDEtunedInd);
lowDEposTunedID = setdiff(lowDEtunedID, negSharpID);

example_angle_tuning_calcium(u, tuningAll.ca, tuningAll.spk, lowDEposTunedID)

%% it seems all about sparsity
%% look at total # of spikes in comparison with dev exp of tuned cells
cd(sprintf('%s%03d',baseDir,mouse))
load(sprintf('UberJK%03dS%02d',mouse,session))
% cellIDs = 


tindCell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
cindSpk = find(u.trials{tindCell(1)}.neuindSession == cID);
planeInd = floor(cID/1000);

testInput = allPredictors{planeInd};
spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cindSpk,:), nan(1,posShift)], u.trials(tindCell)','uniformoutput',false));

figure,
% plot(results.naive(i).allDE(ind), 


%%
%% -> about 2/3 does not look like tuned. Even for the rest 1/3, response amplitude is very low (usually 0.1-0.3 AP / trial)
%%

%% How are these low DE tuned cells distributed, in regards to the depth?

figure,
plot(results.naive(i).allDE(ind), fit.naive(i).cellDepths(ind), 'k.')
hold on
plot([-0.2 1], [350 350], '--', 'color', [0.7 0.7 0.7])
xlabel('Deviance explained')
ylabel('Cell Depth')
xlim([0 0.6])


