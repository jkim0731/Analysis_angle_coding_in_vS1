% for figure 3a
% representative calcium imaging
% calcium
% spikes
% full model (without whisker touch variables)
% devExp (mean)
% correlation
% predictors
%     touch
%     whisking
%         onset
%         amplitude
%         midpoint
%     licking
%     reward
%     pole up sound

%% basic settings
baseDir = 'C:\Data\suite2p\';
mouse = 25;
session = 4;
repeat = 10;
%% dependent settings
ufn = sprintf('UberJK%03dS%02d',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_m45_R', mouse, session);
%% load files
cd(baseDir)
load('glm_results_responseType', 'naive')
cd(sprintf('%s%03d',baseDir, mouse))
load(ufn, 'u') % loading u
load([glmfnBase, '01'], 'allPredictors', 'posShift')

%%
%% Chose plane 3 from JK025 S04
%% Choose examples of sharp tuning
%%
%% settings
plane = 3;

% Image
% figure, imshow(adapthisteq(mat2gray(u.mimg{plane})))
figure, imshow((mat2gray(u.mimg{plane})))
hold on,
plot(u.c2xpoints, u.c2ypoints, 'w--', 'linewidth', 2)

colors = jet(7);
cIDlist = spk.touchID(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
tunedAngles = spk.tunedAngle(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
indList = find(ismember(u.cellNums, cIDlist));

for ci = 1 : length(cIDlist)
    cID = cIDlist(ci);
    cind = indList(ci);
    if tunedAngles(ci) == 0
        scatter(u.cellx(cind), u.celly(cind), 'markeredgecolor', 'none',  'markerfacecolor', 'w', 'linewidth', 2, 'sizedata', 100, 'markerfacealpha', 0.4)
    else
        angleInd = (tunedAngles(ci)-30)/15;
        scatter(u.cellx(cind), u.celly(cind), 'markeredgecolor', 'none',  'markerfacecolor', colors(angleInd,:), 'linewidth', 2, 'sizedata', 100, 'markerfacealpha', 0.4)
    end
end

% scale bar edge. scale bar length 100 um
scaleBarPix = 100 / u.pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(u.mimg{plane});
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
expectedTextLength = 120;
text(sizes(2)-expectedTextLength, margins, [num2str(-round(mean(mean(u.fovdepth{plane})))), ' \mum'], 'fontsize', 18, 'fontweight', 'bold', 'color', 'w')
set(gcf, 'InvertHardCopy', 'off', 'color', 'w');
print('fig3a_example_map','-depsc2')
%%
figure, histogram(naive(1).allDE)
%%
close all
% 3012        3025        3027        3048        3049        3056       3059        3129
cID = 3012;
ci = find(u.cellNums == cID);
traces = get_traces_per_cell(u, cID, allPredictors, naive(1).coeffs{ci}, posShift);
figure('units', 'normalized', 'outerposition', [0.3 0.3 0.4 0.4]), hold on
plot(min_max_normalization(traces.calcium)*4 + 7 , 'color', [0.1 0.8 0.1], 'linewidth', 2)
plot(min_max_normalization(traces.spikes)*2.5 + 5.2, 'color', [0.8 0.1 0.1], 'linewidth', 2)

plot(min_max_normalization(traces.model)*2.5 + 5.2, 'k-', 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,8)) + 3.5, 'color', [0.7 0.7 0.7], 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,63)) + 2.5, 'color', [0.7 0.7 0.7], 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,70)) + 1.8, 'color', [0.7 0.7 0.7], 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,77)) + 1.1, 'color', [0.7 0.7 0.7], 'linewidth', 2)
lickL = min_max_normalization(traces.predictors(:,83));
% lickL(lickL==0) = nan;
lickR = min_max_normalization(traces.predictors(:,87));
% lickR(lickR==0) = nan;
plot(lickL, 'color', [0.9 0.4 0.4], 'linewidth', 2)
plot(lickR, 'color', [0.4 0.4 0.9], 'linewidth', 2)
reward = min_max_normalization(traces.predictors(:,36));
reward(reward<max(reward)) = 0;
reward(isnan(reward)) = 0;
rewarded = find(reward);
rdiffOne = find(diff(find(reward))==1);
reward(rewarded(rdiffOne+1)) = nan;
reward(reward==0) = nan;
sound = min_max_normalization(traces.predictors(:,25));
sound(sound<max(sound)) = nan;
plot(reward+2.3, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'c')
plot(sound+2.3, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'm')
onedF = 4 / max(traces.calcium);
plot([5100 5100], [9 9+onedF], 'color', [0.1 0.8 0.1], 'linewidth', 2)
plot([5120 5120 + u.frameRate*5], [9 9], 'k-', 'linewidth', 2)
xlim([4475 5177]), ylim([0.01 10])

axis off
set(gcf, 'InvertHardCopy', 'off', 'color', 'w');
print('fig3b_example_traces', '-depsc2')
%%
naive(1).allDE(ci)
naive(1).corrVal(ci)
%% Distribution of deviance explained
nonlearner = [5,6,8,10,11,12];
learner = [1:4,7,9];

range = [-0.02:0.02:0.6];
nonlearnerNaive = zeros(6,length(range)-1);
learnerNaive = zeros(6,length(range)-1);
learnerExpert = zeros(6,length(range)-1);

for i = 1 : length(nonlearner)
    nonlearnerNaive(i,:) = histcounts(naive(nonlearner(i)).allDE, range, 'normalization', 'cdf');
end
for i = 1 : length(learner)
    learnerNaive(i,:) = histcounts(naive(learner(i)).allDE, range, 'normalization', 'cdf');
    learnerExpert(i,:) = histcounts(expert(i).allDE, range, 'normalization', 'cdf');
end

figure, hold on
plot(range(2:end), nonlearnerNaive(1,:), 'c-', 'linewidth', 2)
plot(range(2:end), learnerNaive(1,:), 'b-', 'linewidth', 2)
plot(range(2:end), learnerExpert(1,:), 'r-', 'linewidth', 2)
for i = 1 : 6
    plot(range(2:end), nonlearnerNaive(i,:), 'c-', 'linewidth', 2)
end
for i = 1 : 6
    plot(range(2:end), learnerNaive(i,:), 'b-', 'linewidth', 2)
end
for i = 1 : 6
    plot(range(2:end), learnerExpert(i,:), 'r-', 'linewidth', 2)
end

plot([0.1, 0.1], [0, 1], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
xlabel('Deviance explained')
ylabel('Cumulative proportion')
set(gca, 'linewidth', 2', 'fontweight', 'bold', 'fontsize', 10)
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'southeast')

%%
figure, hold on
plot(range(2:end), mean(nonlearnerNaive), 'c-', 'linewidth', 2)
plot(range(2:end), mean(learnerNaive), 'b-', 'linewidth', 2)
plot(range(2:end), mean(learnerExpert), 'r-', 'linewidth', 2)
x = range(2:end);
y = mean(nonlearnerNaive);
e = std(nonlearnerNaive)/sqrt(6);
boundedline(x, y, e, 'c-')
y = mean(learnerNaive);
e = std(learnerNaive)/sqrt(6);
boundedline(x, y, e, 'b-')
y = mean(learnerExpert);
e = std(learnerExpert)/sqrt(6);
boundedline(x, y, e, 'r-')
plot(range(2:end), mean(nonlearnerNaive), 'c-', 'linewidth', 2)
plot(range(2:end), mean(learnerNaive), 'b-', 'linewidth', 2)
plot(range(2:end), mean(learnerExpert), 'r-', 'linewidth', 2)
plot([0.1, 0.1], [0, 1], '--', 'color', [0.7 0.7 0.7], 'linewidth', 2)
xlabel('Deviance explained')
ylabel('Cumulative proportion')
set(gca, 'linewidth', 2', 'fontweight', 'bold', 'fontsize', 10)
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'southeast')


print('fig3c_deviance_explained', '-depsc2')