% NTR comparison between glm and anova test
% question: why is proportion of tuned/touch so much higher after glm
% fitting, compared to anova test only? (previous analysis)
% hypothesis: most of not-tuned response from anova test only have lower 
% spike-explainability by touch.

% Look at 190524 Object angle coding interim summary.pptx

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};
mi = 3;

baseDir = 'Y:\Whiskernas\JK\suite2p\';
% cd(sprintf('%s%03d',baseDir, mice(mi)))
ac = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_all_cell.mat', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
glm = load(sprintf('%s%03d\\JK%03dS%02dangle_tuning.mat', baseDir, mice(mi), mice(mi), sessions{mi}(1)));
load(sprintf('%sglm_results_responseType.mat', baseDir), 'naive')

%% test with tuned cells
idTunedANOVA = ac.spk.touchID(find(ac.spk.tuned));
idTunedGLM = glm.spk.touchID(find(glm.spk.tuned));

idAnovaSetdiffGlm = setdiff(idTunedANOVA,idTunedGLM);
indAnovaSetdiffGlm = find(ismember(ac.spk.touchID,idAnovaSetdiffGlm));
indTunedGlm = find(ismember(ac.spk.touchID,idTunedGLM));

aDEanovasetdiffglm = naive(mi).averageDE(indAnovaSetdiffGlm);
aDEglm = naive(mi).averageDE(indTunedGlm);

mean(aDEglm)
std(aDEglm)/sqrt(length(aDEglm))
mean(aDEanovasetdiffglm)
std(aDEanovasetdiffglm)/sqrt(length(aDEanovasetdiffglm))

%% results: anova tuned have lower DE than glm tuned
% which is actually by definition...

%% test with not-tuned touch cells
idNTRanova = ac.spk.touchID(find(ac.spk.NTR));
idNTRglm = glm.spk.touchID(find(1-glm.spk.tuned));
idNTRanovaSetdiffGlm = setdiff(idNTRanova, idNTRglm);
indNTRanovaSetdiffGlm = find(ismember(ac.spk.touchID,idNTRanovaSetdiffGlm));
indNTRGlm = find(ismember(ac.spk.touchID,idNTRglm));

aDENTRanovaSetdiffGlm = naive(mi).averageDE(indNTRanovaSetdiffGlm);
aDENTRglm = naive(mi).averageDE(indNTRGlm);

mean(aDENTRglm)
std(aDENTRglm)/sqrt(length(aDENTRglm))
mean(aDENTRanovaSetdiffGlm)
std(aDENTRanovaSetdiffGlm)/sqrt(length(aDENTRanovaSetdiffGlm))

%% results: anova NTR have lower DE than glm NTR
% which is actually by definition...

%% graph
figure, hold on
x = rand(length(aDENTRanovaSetdiffGlm),1);
plot(x, aDENTRanovaSetdiffGlm, 'k.')
x = rand(length(aDENTRglm),1);
plot(2 + x, aDENTRglm, 'k.')

% but there are a few high average DE cells with NTR only at anova only
% test. 

%% check if these high aDE anova NTR cells are not GLM tuned cells
tempind = find(aDENTRanovaSetdiffGlm >= 0.1); % 0.1 is the threshold for average DE
idNTRanovaHighADE = idNTRanovaSetdiffGlm(tempind);
templist = intersect(idTunedGLM, idNTRanovaHighADE);

% %% in case there is (happens very rarely, but they do happen), look at anova p values
% tempid = intersect(idTunedGLM, idNTRanovaHighADE);
% anovaP = zeros(length(tempid),1);
% for i = 1 : length(tempid)
%     tempind = ismember(ac.spk.touchID,tempid(i));
%     tempval = ac.spk.val{tempind};
%     anovaVal = cell2mat(tempval);
%     groupLengths = [0; cumsum(cellfun(@length, tempval))];
%     anovaGroup = zeros(groupLengths(end),1);
%     for j = 1 : length(groupLengths)-1
%         anovaGroup(groupLengths(j)+1:groupLengths(j+1)) = deal(j);
%     end
% 
%     anovaP(i) = anova1(anovaVal, anovaGroup, 'off');
% end
% % compared to all other anova p-values of GLM tuned cells
% anovaPglmTuned = zeros(length(indGlmTuned),1);
% for i = 1 : length(indGlmTuned)
%     tempval = ac.spk.val{indGlmTuned(i)};
%     anovaVal = cell2mat(tempval);
%     groupLengths = [0; cumsum(cellfun(@length, tempval))];
%     anovaGroup = zeros(groupLengths(end),1);
%     for j = 1 : length(groupLengths)-1
%         anovaGroup(groupLengths(j)+1:groupLengths(j+1)) = deal(j);
%     end
% 
%     anovaPglmTuned(i) = anova1(anovaVal, anovaGroup, 'off');
% end
% %
% figure, hold on
% x = rand(length(anovaPglmTuned),1);
% plot(x,anovaPglmTuned,'k.')
% x = rand(length(anovaP),1);
% plot(x,anovaP,'r.')

%% There are some cells with anovaP > 0.05. What's going on?
% look at val from glm.spk
tempinds = find(glm.spk.tuned);
anovaPglm = zeros(length(tempinds),1);
for i = 1 : length(tempinds)
    tempval = glm.spk.val{tempinds(i)};
    anovaVal = cell2mat(tempval);
    groupLengths = [0; cumsum(cellfun(@length, tempval))];
    anovaGroup = zeros(groupLengths(end),1);
    for j = 1 : length(groupLengths)-1
        anovaGroup(groupLengths(j)+1:groupLengths(j+1)) = deal(j);
    end

    anovaPglm(i) = anova1(anovaVal, anovaGroup, 'off');
end
figure, histogram(anovaPglm)


%% Ignore this error for now. (test it later with a new run)
% -> No there is not something like this anymore.
% Fixed on 05/30/2019, confirmed on 05/31/2019.
% It was because of different condition on passing permutation test

%% what is the exclusion method DE diff of NTR(anova excluding glm), compared to NTR(glm)?
tempVal = naive(mi).partialSub;
touchDEdiffASG = cellfun(@(x) x(1), tempVal(indNTRanovaSetdiffGlm));
touchDEdiffGLM = cellfun(@(x) x(1), tempVal(indNTRGlm));

figure, hold on
histogram(touchDEdiffASG)
histogram(touchDEdiffGLM)
xlabel('Touch DE diff'), ylabel('Occurrence'), title(sprintf('JK%03d',mice(mi)))


%% What is the response level of NTR(anova excluding glm), compared to NTR(glm)?
figure, hold on,
histogram(ac.spk.NTRamplitude(indNTRanovaSetdiffGlm), -0.5:0.02:1)
histogram(glm.spk.NTamplitude(find(glm.spk.NTamplitude)), -0.5:0.02:1)
xlabel('Response amplitude'), title(sprintf('Not-tuned touch response JK%03d',mice(mi))), ylabel('Occurrence')
legend({'from ANOVA - GLM', 'from GLM'})
%% results: more spread in GLM compared to ASG (anova setdiff glm), but all negatively modulated cells are gone. How come?

%% check partial DE diff of negatively modulated cells from ASG
% look at the ones lower than -0.2 spikes/touch
indTemp = find(ac.spk.NTRamplitude(indNTRanovaSetdiffGlm));
indNeg = indNTRanovaSetdiffGlm(indTemp);
tempCell = naive(mi).partialSub(indNeg);
DEdiffTouchNeg = cellfun(@(x) x(1), tempCell);
figure, histogram(DEdiffTouchNeg)
%%
figure, plot(DEdiffTouchNeg, ac.spk.NTRamplitude(indNeg), 'k.')

%% Results: all of them has DE diff < 0.1. 
%% have a closer look into them, each cell
% starting from the largest inhibited responding cell, look at the trial
% averaged calcium and spikes (and the errors)
% and also the coefficients

% for this, first load uber array
load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mice(mi), mice(mi), sessions{mi}(1)))

%% from all cell anova NTR - glm NTR

[~, indSorted] = sort(ac.spk.NTRamplitude(indNTRanovaSetdiffGlm));

tempi = 100;
indCell = indNTRanovaSetdiffGlm(indSorted(tempi));

Amp = ac.spk.NTRamplitude(indCell)
DE = naive(mi).averageDE(indCell)

testCellnum = u.cellNums(indCell);

indTrial = find(cellfun(@(x) ~isempty(find(x.neuindSession == testCellnum)), u.trials));
indCellSession = find(u.trials{indTrial(1)}.neuindSession == testCellnum);
numFrames = min(cellfun(@(x) size(x.spk,2), u.trials(indTrial)));
dFAllTrials = cell2mat(cellfun(@(x) x.dF(indCellSession,1:numFrames), u.trials(indTrial), 'uniformoutput', false));
spkAllTrials = cell2mat(cellfun(@(x) x.spk(indCellSession,1:numFrames), u.trials(indTrial), 'uniformoutput', false));

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.5 0.3]), 
subplot(131), errorbar(mean(dFAllTrials),std(dFAllTrials)/sqrt(size(dFAllTrials,2))), title('dF')
subplot(132), errorbar(mean(spkAllTrials),std(spkAllTrials)/sqrt(size(dFAllTrials,2))), title('Spikes')
subplot(133), plot(naive(mi).coeffs{indCell}(2:end)), title('coefficients')

%% In comparison, from glm NTR

indNTRGLM = find(glm.spk.NTamplitude);
[~, indSorted] = sort(glm.spk.NTamplitude(indNTRGLM));

tempi = 2;

indCell = indNTRGlm(indSorted(tempi));

Amp = glm.spk.NTamplitude(indNTRGLM(indSorted(tempi)))
DE = naive(mi).averageDE(indCell)

testCellnum = u.cellNums(indCell);

indTrial = find(cellfun(@(x) ~isempty(find(x.neuindSession == testCellnum)), u.trials));
indCellSession = find(u.trials{indTrial(1)}.neuindSession == testCellnum);
numFrames = min(cellfun(@(x) size(x.spk,2), u.trials(indTrial)));
dFAllTrials = cell2mat(cellfun(@(x) x.dF(indCellSession,1:numFrames), u.trials(indTrial), 'uniformoutput', false));
spkAllTrials = cell2mat(cellfun(@(x) x.spk(indCellSession,1:numFrames), u.trials(indTrial), 'uniformoutput', false));

figure('units', 'normalized', 'outerposition', [0.2 0.2 0.5 0.3]), 
subplot(131), errorbar(mean(dFAllTrials),std(dFAllTrials)/sqrt(size(dFAllTrials,2))), title('dF')
subplot(132), errorbar(mean(spkAllTrials),std(spkAllTrials)/sqrt(size(dFAllTrials,2))), title('Spikes')
subplot(133), plot(naive(mi).coeffs{indCell}(2:end)), title('coefficients')








%% Conclusion
% Most of NTRs from setdiff(ANOVA, GLM) comes from very small touch DE diff. 
% What seemed to be touch response was actually not (just modulated by task engagement), or have very small response to touch, or very small number of spikes at all. 
























