% look at the relationship between 
% (1) dKv and dPhi
% (2) dKv and slide distance
% Both in trial-averaged (behavior, Fig2 & 3)
% and in frame-averaged (TPM data, Fig4~6?)

%% Trial-averaged first
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([baseDir, fn], 'groupMdl'); 

dKv = cell(6,1);
dPhi = cell(6,1);
slideDistance = cell(6,1);
poleAngle = cell(6,1);

for mi = 1 : 6
    dKv{mi} = groupMdl{mi}.io.X(:,4);
    dPhi{mi} = groupMdl{mi}.io.X(:,2);
    slideDistance{mi} = groupMdl{mi}.io.X(:,5);
    poleAngle{mi} = groupMdl{mi}.io.Y;
end

figure, 
scatter(cell2mat(dKv), cell2mat(dPhi), 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')
[r,p] = corr(cell2mat(dKv), cell2mat(dPhi))


figure, 
scatter(cell2mat(dKv), cell2mat(slideDistance), 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max(Slide distance)')
[r,p] = corr(cell2mat(dKv), cell2mat(slideDistance))

%% Grouping into pole angles
dkvMat = cell2mat(dKv);
dphiMat = cell2mat(dPhi);
sdMat = cell2mat(slideDistance);
angleMat = cell2mat(poleAngle);
figure, hold on
angleList = 45:15:135;
colorList = jet(7);
legendList = cell(7,1);

for ai = 1 : 7
    aind = find(angleMat == angleList(ai));
    scatter(dkvMat(aind), dphiMat(aind), 10, colorList(ai,:), 'filled')
    legendList{ai} = sprintf('%d^o',angleList(ai));
end
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')
title(sprintf('%s',learned))
axis equal
legend(legendList)

%%



%%
figure, hold on

for ai = 1 : 7
    aind = find(angleMat == angleList(ai));
    scatter(dkvMat(aind), sdMat(aind), 10, colorList(ai,:), 'filled')
end
xlabel('max\Delta\kappa_V')
ylabel('max(slide distance)')
title(sprintf('%s',learned))
axis equal
legend(legendList)
%%



%% going into specific trial (to see what happened during the trial)
figure, 
scatter(cell2mat(dKv), cell2mat(dPhi), 'k.')
[x,y] = ginput(1);
dkvall = cell2mat(dKv);
dphiall = cell2mat(dPhi);
%%
[~,ind] = min((dkvall-x).^2 + (dphiall-y).^2);
groupNums = cumsum(cellfun(@length, dKv));
groupInd = find(groupNums>ind, 1, 'first')-1;
indIngroup = ind-groupNums(groupInd);

%% Can't find which trials were included. Not matching with touch trials from wfa.
%% Select examples from single mouse
%% behavior file is necessary to get the ones before the first lick
%% 

bDir = 'D:\TPM\JK\soloData\';


mi = 2;
mice = [25, 27, 30, 36, 39, 52];
sessions = [4, 3, 3, 1, 1, 3];

mouse = mice(mi);
session = sessions(mi);
wfa = Whisker.WhiskerFinal_2padArray(sprintf('D:\\TPM\\JK\\tracked\\JK%03dS%02d',mouse,session));
touchFrames = cellfun(@(x) x.protractionTFchunksByWhisking, wfa.trials, 'un', 0);

load(sprintf('%sJK%03d\\behavior_JK%03d',bDir, mouse, mouse));
bind = find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',session)), b));
bs = b{bind};
%%
% bsf = 
%%
% dkvCell = cellfun(@(x,y) mean(cellfun(@(z) x.kappaV(z(find(abs(x.kappaV(z)-x.kappaV(z(1))) == max(abs(x.kappaV(z)-x.kappaV(z(1)))) ))) - x.kappaV(z(1)), y)), wfa.trials, touchFrames);
% dphiCell = cellfun(@(x,y) mean(cellfun(@(z) x.phi(z(find(abs(x.phi(z)-x.phi(z(1))) == max(abs(x.phi(z))-x.phi(z(1))) ))), y)), wfa.trials, touchFrames);
% sdCell = cellfun(@(x) mean(cellfun(@max, x.protractionSlideByWhisking)), wfa.trials);

% figure, scatter(dkvCell, dphiCell, 'k.')



% run iterations instead of cellfun because of NaNs
dkv = nan(length(wfa.trials),1);
dphi = nan(length(wfa.trials),1);
slideDistance = nan(length(wfa.trials),1); % slide distance

for wi = 1 : length(wfa.trials)
    wf = wfa.trials{wi};
    if ~isempty(wf.protractionTFchunksByWhisking)
        dkvTemp = nan(length(wf.protractionTFchunksByWhisking),1);
        dphiTemp = nan(length(wf.protractionTFchunksByWhisking),1);
        sdTemp = zeros(length(wf.protractionTFchunksByWhisking),1);
        for i = 1 : length(wf.protractionTFchunksByWhisking)
            if length(wf.protractionTFchunksByWhisking{i}) > 1
                tempInds = wf.protractionTFchunksByWhisking{i};
                [~, maxInd] = max(abs(wf.kappaV(tempInds) - wf.kappaV(tempInds(1))));
                dkvTemp(i) = wf.kappaV(tempInds(maxInd)) - wf.kappaV(tempInds(1));
                
                [~, maxInd] = max(abs(wf.phi(tempInds) - wf.phi(tempInds(1))));
                dphiTemp(i) = wf.phi(tempInds(maxInd)) - wf.phi(tempInds(1));
                
                sdTemp(i) = max(wf.protractionSlideByWhisking{i});
            end
        end
        dkv(wi) = mean(dkvTemp);
        dphi(wi) = mean(dphiTemp);
        slideDistance = mean(sdTemp);
    end    
end
%%
figure, 
ax1 = scatter(dkv, dphi, 'k.');

%%

[x,y] = ginput(1);
[~,ind] = min((dkv-x).^2 + (dphi-y).^2);
wfa.trialNums(ind)
wf = wfa.trials{ind};
wf.poleAngle
%%
figure
subplot(211), hold on
yval = [min(union(wf.theta, wf.phi)), max(union(wf.theta, wf.phi))];
for ti = 1 : length(wf.protractionTFchunksByWhisking)
    tempX = wf.protractionTFchunksByWhisking{ti};    
    patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [0.7 0.7 0.7], 'edgecolor', 'none')
end
plot(wf.theta, 'b')
plot(wf.phi, 'c')

subplot(212), hold on
yval = [min(union(wf.kappaV, wf.kappaH)), max(union(wf.kappaV, wf.kappaH))];
for ti = 1 : length(wf.protractionTFchunksByWhisking)
    tempX = wf.protractionTFchunksByWhisking{ti};    
    patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [0.7 0.7 0.7], 'edgecolor', 'none')
end
plot(wf.kappaH, 'r')
plot(wf.kappaV, 'm')


%% Results: in case of JK025 S04, majority of the 2,4 quadrants were 90 degrees pole or error in video (water stains)


%% plot with angles
figure, hold on
angleList = 45:15:135;
angles = cellfun(@(x) x.poleAngle, wfa.trials);
colorList = jet(7);
for ai = 1 : 7
    aind = find(angles == angleList(ai));
    scatter(dkv(aind), dphi(aind), [], colorList(ai,:));
end




%% for TPM frames
% get features from glm results 
% But these are standardized!
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 
baseFn = 'glmWhisker_lasso_touchCell_NC_'; %JK054S03_R10
experti = [1,2,3,4,7,9];
dkv = cell(length(experti),1);
dphi = cell(length(experti),1);
slideDistance = cell(length(experti),1);
dkh = cell(length(experti),1);
for mi = 1 : length(experti)
    mouse = mice(experti(mi));
    session = sessions{experti(mi)}(1);
    
    dir = sprintf('%s%03d\\',baseDir, mouse);
    fn = sprintf('%sJK%03dS%02d_R10',baseFn, mouse, session);
    load([dir, fn], 'allPredictors')
    
%     whiskerTouchMat = [maxDkappaHMat, maxDkappaVMat, maxDthetaMat, maxDphiMat, maxSlideDistanceMat, maxDurationMat, ...    
%                             thetaAtTouchMat, phiAtTouchMat, kappaHAtTouchMat, kappaVAtTouchMat, arcLengthAtTouchMat, touchCountMat];
% 2: dkv, 4: dphi, 5:slide distance
    dkvTemp = allPredictors{1}(:,(2-1)*3+1);
    dphiTemp = allPredictors{1}(:,(4-1)*3+1);
    slideDistanceTemp = allPredictors{1}(:,(5-1)*3+1);
    dkhTemp = allPredictors{1}(:,1);
    nonaninds = intersect(intersect(find(isfinite(dkvTemp)), find(isfinite(dphiTemp))), intersect(find(isfinite(slideDistanceTemp)),find(isfinite(dkhTemp))) );
    durationTemp = allPredictors{1}(:,(6-1)*3+1);
    duration = durationTemp(find(isfinite(durationTemp)));
    touchinds = find(duration > min(duration));
    inds = intersect(touchinds, nonaninds);
    dkv{mi} = dkvTemp(inds);
    dphi{mi} = dphiTemp(inds);    
    slideDistance{mi} = slideDistanceTemp(inds);    
    dkh{mi} = dkhTemp(inds);
end

figure,
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')
[r,p] = corr(cell2mat(dkv), cell2mat(dphi))

figure, 
scatter(cell2mat(dkv), cell2mat(slideDistance), 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max(Slide distance)')
[r,p] = corr(cell2mat(dkv), cell2mat(slideDistance))


%% figure
i = 7;
figure,
subplot(121),  scatter(dkv{i}, dphi{i}, 'k.')
subplot(122), scatter(dkv{i}, slideDistance{i}, 'k.')


%% group by dKh
figure, hold on
dkvMat = cell2mat(dkv);
dphiMat = cell2mat(dphi);
dkhMat = cell2mat(dkh);
dkhThresh = [2 1 0 -1 -2 -3];
colorList = jet(length(dkhThresh));
legendList = cell(length(dkhThresh),1);
for ti = 1:length(dkhThresh)
    ind = find(dkhMat<dkhThresh(ti));
    scatter(dkvMat(ind), dphiMat(ind), 10, colorList(ti,:), 'filled')
    legendList{ti} = sprintf('dKh < %d', dkhThresh(ti));
end
legend(legendList)

%% group by dKh
figure, hold on
dkvMat = cell2mat(dkv);
dphiMat = cell2mat(dphi);
dkhMat = cell2mat(dkh);
dkhThresh = [-4 -3 -2 -1 0];
colorList = jet(length(dkhThresh));
legendList = cell(length(dkhThresh),1);
for ti = 1:length(dkhThresh)
    ind = find(dkhMat>dkhThresh(ti));
    scatter(dkvMat(ind), dphiMat(ind), 10, colorList(ti,:), 'filled')
    legendList{ti} = sprintf('dKh > %d', dkhThresh(ti));
end
legend(legendList)



%% from uber, calculate Mean response in each frame 
% only take the first plane frames (there are 4 planes in TPM)











%% values are standardized, so hard to interpret. 
%% take frame-mean from uber arrays

plane = 1; % the plane to calculate TPM time frames

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    
    touchOnsetCell = cellfun(@(x) x.protractionTouchOnsetFramesByWhisking{plane}, u.trials, 'un', 0)';
    [~,~,touchOnsetInd] = cellfun(@(x) unique(x), touchOnsetCell, 'un', 0);
    dkvCell{mi} = cell(1,length(u.trials));
    dphiCell{mi} = cell(1,length(u.trials));
    slideDistanceCell{mi} = cell(1,length(u.trials));
    for i = 1 : length(u.trials)
        if ~isempty(touchOnsetInd{i})
            currTrial = u.trials{i};
            touchFrameInds = unique(touchOnsetInd{i});
            tempDkv = zeros(1,length(touchFrameInds));
            tempDphi = zeros(1,length(touchFrameInds));
            tempSd = zeros(1,length(touchFrameInds));
            for j = 1 : length(touchFrameInds)
                tempFrameInds = find(touchOnsetInd{i} == touchFrameInds(j));
                tempDkv(j) = mean(currTrial.protractionTouchDKappaVByWhisking(tempFrameInds));
                tempDphi(j) = mean(currTrial.protractionTouchDPhiByWhisking(tempFrameInds));
                tempSd(j) = mean(currTrial.protractionTouchSlideDistanceByWhisking(tempFrameInds));
            end
            dkvCell{mi}{i} = tempDkv;
            dphiCell{mi}{i} = tempDphi;
            slideDistanceCell{mi}{i} = tempSd;
        end
    end    
end
%%
dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));

inds = intersect(intersect(find(isfinite(dkv)), find(isfinite(dphi))), find(isfinite(slideDistance)));

figure, hold on
scatter(dkv(inds), dphi(inds), 0.001, 'k.')

ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

[r, p] = corr(dkv(inds)', dphi(inds)')

figure, hold on
scatter(dkv(inds), slideDistance(inds), 0.001, 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max(slide distance)')

[r, p] = corr(dkv(inds)', slideDistance(inds)')







%% Look at the scatter and cumulative proportion of slide distance in all mice
%% all naive
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));



inds = intersect(intersect(find(isfinite(dkv)), find(isfinite(dphi))), find(isfinite(slideDistance)));

figure, hold on
scatter(dkv(inds), dphi(inds), 0.001, 'k.')

ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

[r, p] = corr(dkv(inds)', dphi(inds)')

figure, hold on
scatter(dkv(inds), slideDistance(inds), 0.001, 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max(slide distance)')

[r, p] = corr(dkv(inds)', slideDistance(inds)')







%% somethings weird - how about from uber? 
%% (Uber is from ALL touches!!!)
load(sprintf('UberJK%03dS%02d_NC',25,4))
tnums = cellfun(@(x) ones(1,length(x.protractionTouchChunksByWhisking)) * x.trialNum, u.trials, 'uniformoutput', false)';
dkv = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
dphi = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
slideDistance = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';

figure, 
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')
figure,
scatter(cell2mat(dkv), cell2mat(slideDistance), 'k.')

%% how about the mean in each trial from uber?
dkv = cellfun(@(x) mean(x.protractionTouchDKappaVByWhisking), u.trials, 'uniformoutput', false)';
dphi = cellfun(@(x) mean(x.protractionTouchDPhiByWhisking), u.trials, 'uniformoutput', false)';
slideDistance = cellfun(@(x) mean(x.protractionTouchSlideDistanceByWhisking), u.trials, 'uniformoutput', false)';
figure, 
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')
figure,
scatter(cell2mat(dkv), cell2mat(slideDistance), 'k.')

%% look into details of the examples of the quadrants of dphi vs dkv.
% from a single uber array - JK025 S04
load('D:\TPM\JK\suite2p\025\UberJK025S04_NC.mat')
tnums = cellfun(@(x) ones(1,length(x.protractionTouchChunksByWhisking)) * x.trialNum, u.trials, 'uniformoutput', false)';
dkv = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
dphi = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
slideDistance = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';

figure, 
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')

%% pick one point, and get the value

[x,y] = ginput(1)
%%
[~, ind] = min((cell2mat(dkv)-x).^2 + (cell2mat(dphi)-y).^2);
tnumArray = cell2mat(tnums);
tnum = tnumArray(ind);
tnumInd = find(cellfun(@(x) ismember(tnum,x), tnums)); % index of u.trials
allTouchesBefore = sum(cellfun(@length, tnums(1:tnumInd-1)));
touchInd = ind - allTouchesBefore; % index of u.trials.protractionTouchChunksByWhisking

dkv{tnumInd}(touchInd)
dphi{tnumInd}(touchInd)

u.trialNums(tnumInd)
slideDistance{tnumInd}(touchInd)

%% Results: sticking touch - because of friction, not sliding









%% Change thresholds in slide distance and plot on the dphi vs dkv scatter
%% all naive - all touches
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));
angleCell = cell(1,length(mice));
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
    angleCell{mi} = cellfun(@(x) ones(1,length(x.protractionTouchChunksByWhisking)) *x.angle, u.trials, 'un', 0)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));
angle = cell2mat(cellfun(@(x) cell2mat(x), angleCell, 'un', 0));
%%
thList = [1 0.5 0.2 0.1 0.05];

    figure, hold on
colors = jet(length(thList));
legendList = cell(1,length(thList));
legendH = struct;
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    legendH.h(i) = scatter(dkv(sdInd), dphi(sdInd), 0.001, colors(i,:), '.');
    legendList{i} = sprintf('slide distance < %.2f', thList(i));
end
legend([legendH.h], legendList, 'autoupdate', 'off')
    scatter(dkv, dphi, 0.001, 'k.')
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    scatter(dkv(sdInd), dphi(sdInd), 0.001, colors(i,:), '.')
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

%% Applying the threshold
sdThreshold = 1; 
sdInd = find(slideDistance >= sdThreshold);
figure
scatter(dkv(sdInd), dphi(sdInd), 0.01, 'k.')
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')
[r, p] = corr(dkv(sdInd), dphi(sdInd))

%% distribution of slide distance
sdRange = [-1:0.1:5];
figure, histogram(slideDistance, sdRange, 'normalization', 'cdf')
xlabel('Slide distance')
ylabel('Cumulative proportion')



%% Is there a difference in pole angle?
angleList = 45:15:135;
legendList = cell(length(angleList),1);
colorList = jet(7);
figure, hold on
for ai = 1 : 7
    angleInd = find(angle == angleList(ai));
    scatter(dkv(angleInd), dphi(angleInd), 0.0001, colorList(ai,:), '.')
    legendList{ai} = sprintf('%d^o',angleList(ai));
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

legend(legendList)


%% Look at the scatter and cumulative proportion of slide distance in all mice
%% all naive
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));

thList = [1 0.5 0.2 0.1 0.05];

    figure, hold on
colors = jet(length(thList));
legendList = cell(1,length(thList));
legendH = struct;
for i = 1 : length(thList)
%     color = [0.5+(i)*0.1, 0.1*i, 0.1*i];

    sdThreshold = thList(i);

    sdInd = find(slideDistance < sdThreshold);

    legendH.h(i) = scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.');
    legendList{i} = sprintf('slide distance < %.2f', thList(i));
end
legend([legendH.h], legendList, 'autoupdate', 'off')

    scatter(dkv, dphi, 'k.')
for i = 1 : length(thList)


    sdThreshold = thList(i);

    sdInd = find(slideDistance < sdThreshold);

    scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.')
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

%%
sdRange = [0:0.1:3];
figure, histogram(slideDistance, sdRange, 'normalization', 'cdf')
xlabel('Slide distance')
ylabel('Cumulative proportion')

%% from all experts

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));
dkvCell = cell(1,length(indExpert));
dphiCell = cell(1,length(indExpert));
slideDistanceCell = cell(1,length(indExpert));

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    session = sessions{mi}(2);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));

thList = [1 0.5 0.2 0.1 0.05];
%%
    figure, hold on
colors = jet(length(thList));
legendList = cell(1,length(thList));
legendH = struct;
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    legendH.h(i) = scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.');
    legendList{i} = sprintf('slide distance < %.2f', thList(i));
end
legend([legendH.h], legendList, 'autoupdate', 'off')

    scatter(dkv, dphi, 'k.')
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.')
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')
%%
sdRange = [0:0.1:3];
figure, histogram(slideDistance, sdRange, 'normalization', 'cdf')
xlabel('Slide distance')
ylabel('Cumulative proportion')


%% from matching naive

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));
dkvCell = cell(1,length(indExpert));
dphiCell = cell(1,length(indExpert));
slideDistanceCell = cell(1,length(indExpert));

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));

thList = [1 0.5 0.2 0.1 0.05];

    figure, hold on
colors = jet(length(thList));
legendList = cell(1,length(thList));
legendH = struct;
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    legendH.h(i) = scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.');
    legendList{i} = sprintf('slide distance < %.2f', thList(i));
end
legend([legendH.h], legendList, 'autoupdate', 'off')

    scatter(dkv, dphi, 'k.')
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(slideDistance < sdThreshold);

    scatter(dkv(sdInd), dphi(sdInd), [], colors(i,:), '.')
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

sdRange = [0:0.1:3];
figure, histogram(slideDistance, sdRange, 'normalization', 'cdf')
xlabel('Slide distance')
ylabel('Cumulative proportion')


%% Results: Similar thersholds across mice and learning seems to be OK (slide distance < 0.2)
% Their distribution is similar (~ 60% with < 0.2 slide distance)
% But, there seems to be a huge difference in dkv distribution after learning.


%% Compare distribution of dkv and dphi before and after learning
% between expert and matching naives

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));
dkvCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dphiCell = cell(2,length(indExpert));
slideDistanceCell = cell(2,length(indExpert));

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    for i = 1 : 2
        session = sessions{mi}(i);

        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')

        dkvCell{i,ei} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
        dphiCell{i,ei} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
        slideDistanceCell{i,ei} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
    end
end

dkvNaive = cell2mat(cellfun(@(x) cell2mat(x), dkvCell(1,:), 'uniformoutput', false));
dphiNaive = cell2mat(cellfun(@(x) cell2mat(x), dphiCell(1,:), 'uniformoutput', false));
slideDistanceNaive = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell(1,:), 'uniformoutput', false));

dkvExpert = cell2mat(cellfun(@(x) cell2mat(x), dkvCell(2,:), 'uniformoutput', false));
dphiExpert = cell2mat(cellfun(@(x) cell2mat(x), dphiCell(2,:), 'uniformoutput', false));
slideDistanceExpert = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell(2,:), 'uniformoutput', false));

dkvRange = [-0.1:0.01:0.11];
histDkvNaive = histcounts(dkvNaive, dkvRange, 'normalization', 'cdf');
histDkvExpert = histcounts(dkvExpert, dkvRange, 'normalization', 'cdf');

dphiRange = [-10:0.1:10.1];
histDphiNaive = histcounts(dphiNaive, dphiRange, 'normalization', 'cdf');
histDphiExpert = histcounts(dphiExpert, dphiRange, 'normalization', 'cdf');

sdRange = [0:0.1:3.1];
histSdNaive = histcounts(slideDistanceNaive, sdRange, 'normalization', 'cdf');
histSdExpert = histcounts(slideDistanceExpert, sdRange, 'normalization', 'cdf');


figure,
subplot(131), hold on
plot(dkvRange(1:end-1), histDkvNaive)
plot(dkvRange(1:end-1), histDkvExpert)
title('\Delta\kappa_V')
ylabel('Cumulative proportion')

subplot(132), hold on
plot(dphiRange(1:end-1), histDphiNaive)
plot(dphiRange(1:end-1), histDphiExpert)
title('\Delta\phi')

subplot(133), hold on
plot(sdRange(1:end-1), histSdNaive)
plot(sdRange(1:end-1), histSdExpert)
title('Slide distance')

legend({'Matching naive (n = 6)', 'Expert (n = 6)'})
%%
dkv = cell(2,1);
dphi = cell(2,1);
for i = 1 : 2
    dkv{i} = (cellfun(@(x) cell2mat(x), dkvCell(i,:), 'uniformoutput', false));
    dphi{i}= (cellfun(@(x) cell2mat(x), dphiCell(i,:), 'uniformoutput', false));
end
dkvRange = linspace(-0.05, 0.05, 100);
dphiRange = linspace(-5, 5, 100);
dkvDiff = zeros(1,6);
dphiDiff = zeros(1,6);
for mi = 1 : 6
    dkvHistNaive = histcounts(dkv{1}{mi}, dkvRange, 'normalization', 'cdf');
    dkvHistExpert = histcounts(dkv{2}{mi}, dkvRange, 'normalization', 'cdf');
    dkvDiff(mi) = sum(abs(dkvHistNaive - dkvHistExpert));
    
    dphiHistNaive = histcounts(dphi{1}{mi}, dphiRange, 'normalization', 'cdf');
    dphiHistExpert = histcounts(dphi{2}{mi}, dphiRange, 'normalization', 'cdf');
    dphiDiff(mi) = sum(abs(dphiHistNaive - dphiHistExpert));
    
end

figure, scatter(dkvDiff, dphiDiff, 'k', 'filled')
xlabel('Diff in \Delta\kappa_V')
ylabel('Diff in \Delta\phi')
[~,p] = ttest(dkvDiff, dphiDiff)


%% Results: Efficient touch after learning. 
%% The change in distribution is similar between ??V and ??.

%% How about in trial-averaged?

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));
dkvCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dphiCell = cell(2,length(indExpert));
slideDistanceCell = cell(2,length(indExpert));

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    for i = 1 : 2
        session = sessions{mi}(i);

        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')

        dkvCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDKappaVByWhisking), u.trials, 'uniformoutput', false)';
        dphiCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDPhiByWhisking), u.trials, 'uniformoutput', false)';
        slideDistanceCell{i,ei} = cellfun(@(x) mean(x.protractionTouchSlideDistanceByWhisking), u.trials, 'uniformoutput', false)';
    end
end
%%
dkvNaive = cell2mat(cellfun(@(x) cell2mat(x), dkvCell(1,:), 'uniformoutput', false));
dphiNaive = cell2mat(cellfun(@(x) cell2mat(x), dphiCell(1,:), 'uniformoutput', false));
slideDistanceNaive = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell(1,:), 'uniformoutput', false));
nonanInd = intersect(intersect(find(isfinite(dkvNaive)), find(isfinite(dphiNaive))), find(isfinite(slideDistanceNaive)));
dkvNaive = dkvNaive(nonanInd);
dphiNaive = dphiNaive(nonanInd);
slideDistanceNaive = slideDistanceNaive(nonanInd);

dkvExpert = cell2mat(cellfun(@(x) cell2mat(x), dkvCell(2,:), 'uniformoutput', false));
dphiExpert = cell2mat(cellfun(@(x) cell2mat(x), dphiCell(2,:), 'uniformoutput', false));
slideDistanceExpert = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell(2,:), 'uniformoutput', false));
nonanInd = intersect(intersect(find(isfinite(dkvExpert)), find(isfinite(dphiExpert))), find(isfinite(slideDistanceExpert)));
dkvExpert = dkvExpert(nonanInd);
dphiExpert = dphiExpert(nonanInd);
slideDistanceExpert = slideDistanceExpert(nonanInd);


dkvRange = [-0.05:0.001:0.051];
histDkvNaive = histcounts(dkvNaive, dkvRange, 'normalization', 'cdf');
histDkvExpert = histcounts(dkvExpert, dkvRange, 'normalization', 'cdf');

dphiRange = [-5:0.1:5.1];
histDphiNaive = histcounts(dphiNaive, dphiRange, 'normalization', 'cdf');
histDphiExpert = histcounts(dphiExpert, dphiRange, 'normalization', 'cdf');

sdRange = [0:0.1:2.1];
histSdNaive = histcounts(slideDistanceNaive, sdRange, 'normalization', 'cdf');
histSdExpert = histcounts(slideDistanceExpert, sdRange, 'normalization', 'cdf');


figure,
subplot(131), hold on
plot(dkvRange(1:end-1), histDkvNaive)
plot(dkvRange(1:end-1), histDkvExpert)
title('\Delta\kappa_V')
ylabel('Cumulative proportion')

subplot(132), hold on
plot(dphiRange(1:end-1), histDphiNaive)
plot(dphiRange(1:end-1), histDphiExpert)
title('\Delta\phi')

subplot(133), hold on
plot(sdRange(1:end-1), histSdNaive)
plot(sdRange(1:end-1), histSdExpert)
title('Slide distance')

legend({'Matching naive (n = 6)', 'Expert (n = 6)'})

%% distribution difference similarity (Between dkv and dphi)

dkv = cell(2,1);
dphi = cell(2,1);
for i = 1 : 2
    dkv{i} = (cellfun(@(x) cell2mat(x), dkvCell(i,:), 'uniformoutput', false));
    dphi{i}= (cellfun(@(x) cell2mat(x), dphiCell(i,:), 'uniformoutput', false));
end
dkvRange = linspace(-0.05, 0.05, 100);
dphiRange = linspace(-5, 5, 100);
dkvDiff = zeros(1,6);
dphiDiff = zeros(1,6);
for mi = 1 : 6
    dkvHistNaive = histcounts(dkv{1}{mi}, dkvRange, 'normalization', 'cdf');
    dkvHistExpert = histcounts(dkv{2}{mi}, dkvRange, 'normalization', 'cdf');
    dkvDiff(mi) = sum(abs(dkvHistNaive - dkvHistExpert));
    
    dphiHistNaive = histcounts(dphi{1}{mi}, dphiRange, 'normalization', 'cdf');
    dphiHistExpert = histcounts(dphi{2}{mi}, dphiRange, 'normalization', 'cdf');
    dphiDiff(mi) = sum(abs(dphiHistNaive - dphiHistExpert));
    
end

figure, scatter(dkvDiff, dphiDiff, 'k', 'filled')
xlabel('Diff in \Delta\kappa_V')
ylabel('Diff in \Delta\phi')
[~,p] = ttest(dkvDiff, dphiDiff)
%% How does it look when applying the slide distance threshold? - All touches
%% Only from slide distance >= 0.2
sdThreshold = 0.2;

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    slideDistanceCell{mi} = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));


sdInd = find(slideDistance >= sdThreshold);

    figure, hold on
    scatter(dkv(sdInd), dphi(sdInd), 0.001, 'k.')

ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

title(sprintf('Slide distance threshold %.1f', sdThreshold))

nonanind = union(find(isfinite(dkv)), find(isfinite(dphi)));
corrInd = intersect(nonanind, sdInd);
[r, p] = corr(dkv(corrInd)', dphi(corrInd)')



%% How does it look when applying the slide distance threshold? - TPM frames
%% Only from slide distance >= 0.2
sdThreshold = 0.5;
plane = 1; % the plane to calculate TPM time frames

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
slideDistanceCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    
    touchOnsetCell = cellfun(@(x) x.protractionTouchOnsetFramesByWhisking{plane}, u.trials, 'un', 0)';
    [~,~,touchOnsetInd] = cellfun(@(x) unique(x), touchOnsetCell, 'un', 0);
    dkvCell{mi} = cell(1,length(u.trials));
    dphiCell{mi} = cell(1,length(u.trials));
    slideDistanceCell{mi} = cell(1,length(u.trials));
    for i = 1 : length(u.trials)
        if ~isempty(touchOnsetInd{i})
            currTrial = u.trials{i};
            touchFrameInds = unique(touchOnsetInd{i});
            tempDkv = zeros(1,length(touchFrameInds));
            tempDphi = zeros(1,length(touchFrameInds));
            tempSd = zeros(1,length(touchFrameInds));
            for j = 1 : length(touchFrameInds)
                tempFrameInds = find(touchOnsetInd{i} == touchFrameInds(j));
                tempDkv(j) = mean(currTrial.protractionTouchDKappaVByWhisking(tempFrameInds));
                tempDphi(j) = mean(currTrial.protractionTouchDPhiByWhisking(tempFrameInds));
                tempSd(j) = mean(currTrial.protractionTouchSlideDistanceByWhisking(tempFrameInds));
            end
            dkvCell{mi}{i} = tempDkv;
            dphiCell{mi}{i} = tempDphi;
            slideDistanceCell{mi}{i} = tempSd;
        end
    end    
end

%%
dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
slideDistance = cell2mat(cellfun(@(x) cell2mat(x), slideDistanceCell, 'uniformoutput', false));


sdInd = find(slideDistance >= sdThreshold);

    figure, hold on
    scatter(dkv(sdInd), dphi(sdInd), 0.001, 'k.')

ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')

title(sprintf('Slide distance threshold %.1f', sdThreshold))

nonanind = union(find(isfinite(dkv)), find(isfinite(dphi)));
corrInd = intersect(nonanind, sdInd);
[r, p] = corr(dkv(corrInd)', dphi(corrInd)')





%% Results: Still quite many of 2nd and 4th quadrant remains.
%% How about applying touch duration threshold?

%% duration vs slide distance
%% from trial-averaged
baseDir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
sessionGroupName = 'TwoExpert';
whiskDir = 'protraction'; % 'protraction' or 'all' to choose which touches to use
touchOrder = 'all'; % 'first' or 'all' to choose which touches pre-decision
% yOut = 'Ttype'; % can be 'ttype' (45 vs 135), 'discrete' (45:15:135) or 'choice' (lick right probability)
Xhow = 'Mean'; %can be 'mean' or 'individual' so each touch is counted 

learned = 'Naive'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([baseDir, fn], 'groupMdl'); 

% dKv = cell(6,1);
% dPhi = cell(6,1);
slideDistance = cell(6,1);
duration = cell(6,1);

for mi = 1 : 6
%     dKv{mi} = groupMdl{mi}.io.X(:,4);
%     dPhi{mi} = groupMdl{mi}.io.X(:,2);
    slideDistance{mi} = groupMdl{mi}.io.X(:,5);
    duration{mi} = groupMdl{mi}.io.X(:,6);
end

% figure, 
% scatter(cell2mat(dKv), cell2mat(dPhi), 'k.')
% xlabel('max\Delta\kappa_V')
% ylabel('max\Delta\phi')
% [r,p] = corr(cell2mat(dKv), cell2mat(dPhi))
% 
% 
% figure, 
% scatter(cell2mat(dKv), cell2mat(slideDistance), 'k.')
% xlabel('max\Delta\kappa_V')
% ylabel('max(Slide distance)')
% [r,p] = corr(cell2mat(dKv), cell2mat(slideDistance))
% 
figure, 
scatter(cell2mat(slideDistance), cell2mat(duration), 'k.')

xlabel('max(slide distance)')
ylabel('max(duration)')
[r,p] = corr(cell2mat(slideDistance), cell2mat(duration))

%% from each TPM frame

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

slideDistanceCell = cell(1,length(mice));
touchDurationCell = cell(1,length(mice));
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    touchDurationCell{mi} = cell2mat(cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials, 'uniformoutput', false)');
    slideDistanceCell{mi} = cell2mat(cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)');
end
%%
sdTemp = cell2mat(slideDistanceCell')';
tdTemp = cell2mat(touchDurationCell')';
nonanind = intersect(find(isfinite(sdTemp)), find(isfinite(tdTemp)));
sd = sdTemp(nonanind);
td = tdTemp(nonanind);
figure,
scatter(sd, td, 0.001, 'k.')
xlabel('max(slide distance)')
ylabel('max(touch duration)')
[r,p] = corr(sd, td)

%% Results: slide distance and touch duration is well-correlated. 


%% Do the same thing from slide distance
%% all naive
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
durationCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    durationCell{mi} = cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));
duration = cell2mat(cellfun(@(x) cell2mat(x), durationCell, 'uniformoutput', false));

thList = [0.05 0.04 0.03 0.02 0.01];

    figure, hold on
colors = jet(length(thList));
legendList = cell(1,length(thList));
legendH = struct;
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(duration < sdThreshold);

    legendH.h(i) = scatter(dkv(sdInd), dphi(sdInd), 0.01, colors(i,:), '.');
    legendList{i} = sprintf('slide distance < %.2f', thList(i));
end
legend([legendH.h], legendList, 'autoupdate', 'off')

    scatter(dkv, dphi, 0.01, 'k.')
for i = 1 : length(thList)
    sdThreshold = thList(i);
    sdInd = find(duration < sdThreshold);

    scatter(dkv(sdInd), dphi(sdInd), 0.01, colors(i,:), '.')
end
ylim([-15 15])
xlim([-0.1 0.1])
xlabel('max\Delta\kappa_V')
ylabel('max\Delta\phi')


tdRange = [0:0.002:0.152];
figure, histogram(duration, tdRange, 'normalization', 'cdf')
xlabel('Touch duration')
ylabel('Cumulative proportion')

%% Results: using touch duration as the threshold won't work.
%% It's all distributed.

%% Touch duration between naive and expert
%% from trial averaging
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));
durationCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    for i = 1 : 2
        session = sessions{mi}(i);

        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')

        durationCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDurationByWhisking), u.trials, 'uniformoutput', false)';
        
    end
end

duration = cell(2,1);

for i = 1 : 2
     temp = (cellfun(@(x) cell2mat(x), durationCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     duration{i} = temp;
end

durationRange = linspace(0, 0.06, 100);

histDuration = zeros(6,length(durationRange)-1,2);
for mi = 1 : 6
    histDuration(mi,:,1) = histcounts(duration{1}{mi}, durationRange, 'normalization', 'cdf');
    histDuration(mi,:,2) = histcounts(duration{2}{mi}, durationRange, 'normalization', 'cdf');    
end

figure, hold on
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,1))), std(squeeze(histDuration(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,2))), std(squeeze(histDuration(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('Protraction touch duration')
ylabel('Cumulative proportion')
legend({'Matching naive (n = 6)', 'Expert (n = 6)'})



%% Efficiency summary
%% trial averaged

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));

% parameters controlled by movement
dThetaCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dKhCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
durationCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
% parameters resulting from the controlled movement
dPhiCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dKvCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
slideDistanceCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    for i = 1 : 2
        session = sessions{mi}(i);

        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
        
        dThetaCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDThetaByWhisking), u.trials, 'uniformoutput', false)';
        dKhCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDKappaHByWhisking), u.trials, 'uniformoutput', false)';
        durationCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDurationByWhisking), u.trials, 'uniformoutput', false)';
        
        dPhiCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDPhiByWhisking), u.trials, 'uniformoutput', false)';
        dKvCell{i,ei} = cellfun(@(x) mean(x.protractionTouchDKappaVByWhisking), u.trials, 'uniformoutput', false)';
        slideDistanceCell{i,ei} = cellfun(@(x) mean(x.protractionTouchSlideDistanceByWhisking), u.trials, 'uniformoutput', false)';
    end
end

dTheta = cell(2,1);
dKh = cell(2,1);
duration = cell(2,1);

dPhi = cell(2,1);
dKv = cell(2,1);
slideDistance = cell(2,1);

for i = 1 : 2
     temp = (cellfun(@(x) cell2mat(x), dThetaCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dTheta{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dKhCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dKh{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), durationCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     duration{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dPhiCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dPhi{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dKvCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dKv{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), slideDistanceCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     slideDistance{i} = temp;
end
%%
% max(cell2mat(dTheta{1}))
% min(cell2mat(dTheta{1}))
dThetaRange = linspace(0,10,100);
dKhRange = linspace(-0.10, 0, 100);
durationRange = linspace(0, 0.05, 100);
dPhiRange = linspace(-4, 4, 100);
dKvRange = linspace(-0.03, 0.03, 100);
sdRange = linspace(0,2,100); 

histDTheta = zeros(6,length(dThetaRange)-1,2);
histDKh = zeros(6,length(dKhRange)-1,2);
histDuration = zeros(6,length(durationRange)-1,2);

histDPhi = zeros(6,length(dPhiRange)-1,2);
histDKv = zeros(6,length(dKvRange)-1,2);
histSD = zeros(6,length(sdRange)-1,2);
for mi = 1 : 6
    histDTheta(mi,:,1) = histcounts(dTheta{1}{mi}, dThetaRange, 'normalization', 'cdf');
    histDTheta(mi,:,2) = histcounts(dTheta{2}{mi}, dThetaRange, 'normalization', 'cdf');
    histDKh(mi,:,1) = histcounts(dKh{1}{mi}, dKhRange, 'normalization', 'cdf');
    histDKh(mi,:,2) = histcounts(dKh{2}{mi}, dKhRange, 'normalization', 'cdf');    
    histDuration(mi,:,1) = histcounts(duration{1}{mi}, durationRange, 'normalization', 'cdf');
    histDuration(mi,:,2) = histcounts(duration{2}{mi}, durationRange, 'normalization', 'cdf');
    
    histDPhi(mi,:,1) = histcounts(dPhi{1}{mi}, dPhiRange, 'normalization', 'cdf');
    histDPhi(mi,:,2) = histcounts(dPhi{2}{mi}, dPhiRange, 'normalization', 'cdf');
    histDKv(mi,:,1) = histcounts(dKv{1}{mi}, dKvRange, 'normalization', 'cdf');
    histDKv(mi,:,2) = histcounts(dKv{2}{mi}, dKvRange, 'normalization', 'cdf');
    histSD(mi,:,1) = histcounts(slideDistance{1}{mi}, sdRange, 'normalization', 'cdf');
    histSD(mi,:,2) = histcounts(slideDistance{2}{mi}, sdRange, 'normalization', 'cdf');
    
end

figure, 
subplot(231), hold on
shadedErrorBar(dThetaRange(1:end-1), mean(squeeze(histDTheta(:,:,1))), std(squeeze(histDTheta(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dThetaRange(1:end-1), mean(squeeze(histDTheta(:,:,2))), std(squeeze(histDTheta(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(\Delta\theta)')
ylabel('Cumulative proportion')

subplot(232), hold on
shadedErrorBar(dKhRange(1:end-1), mean(squeeze(histDKh(:,:,1))), std(squeeze(histDKh(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dKhRange(1:end-1), mean(squeeze(histDKh(:,:,2))), std(squeeze(histDKh(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(\Delta\kappa_H)')

subplot(233), hold on
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,1))), std(squeeze(histDuration(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,2))), std(squeeze(histDuration(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(protraction touch duration)')
legend({'Matching naive (n = 6)', 'Expert (n = 6)'}, 'location', 'southeast')

subplot(234), hold on
shadedErrorBar(dPhiRange(1:end-1), mean(squeeze(histDPhi(:,:,1))), std(squeeze(histDPhi(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dPhiRange(1:end-1), mean(squeeze(histDPhi(:,:,2))), std(squeeze(histDPhi(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(\Delta\phi)')
ylabel('Cumulative proportion')

subplot(235), hold on
shadedErrorBar(dKvRange(1:end-1), mean(squeeze(histDKv(:,:,1))), std(squeeze(histDKv(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dKvRange(1:end-1), mean(squeeze(histDKv(:,:,2))), std(squeeze(histDKv(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(\Delta\kappa_V)')

subplot(236), hold on
shadedErrorBar(sdRange(1:end-1), mean(squeeze(histSD(:,:,1))), std(squeeze(histSD(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(sdRange(1:end-1), mean(squeeze(histSD(:,:,2))), std(squeeze(histSD(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(slide distance)')




%% all touches
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

indExpert = find(cellfun(@(x) length(x)==2, sessions));

% parameters controlled by movement
dThetaCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dKhCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
durationCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
% parameters resulting from the controlled movement
dPhiCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
dKvCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts
slideDistanceCell = cell(2,length(indExpert)); % 1 for matching naive, 2 for experts

for ei = 1 : length(indExpert)
    mi = indExpert(ei);
    mouse = mice(mi);
    for i = 1 : 2
        session = sessions{mi}(i);

        load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
        
        dThetaCell{i,ei} = cellfun(@(x) (x.protractionTouchDThetaByWhisking), u.trials, 'uniformoutput', false)';
        dKhCell{i,ei} = cellfun(@(x) (x.protractionTouchDKappaHByWhisking), u.trials, 'uniformoutput', false)';
        durationCell{i,ei} = cellfun(@(x) (x.protractionTouchDurationByWhisking), u.trials, 'uniformoutput', false)';
        
        dPhiCell{i,ei} = cellfun(@(x) (x.protractionTouchDPhiByWhisking), u.trials, 'uniformoutput', false)';
        dKvCell{i,ei} = cellfun(@(x) (x.protractionTouchDKappaVByWhisking), u.trials, 'uniformoutput', false)';
        slideDistanceCell{i,ei} = cellfun(@(x) (x.protractionTouchSlideDistanceByWhisking), u.trials, 'uniformoutput', false)';
    end
end

dTheta = cell(2,1);
dKh = cell(2,1);
duration = cell(2,1);

dPhi = cell(2,1);
dKv = cell(2,1);
slideDistance = cell(2,1);

for i = 1 : 2
     temp = (cellfun(@(x) cell2mat(x), dThetaCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dTheta{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dKhCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dKh{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), durationCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     duration{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dPhiCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dPhi{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), dKvCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     dKv{i} = temp;
     
     temp = (cellfun(@(x) cell2mat(x), slideDistanceCell(i,:), 'uniformoutput', false));
     for ti = 1 : length(temp)
         tempInd = find(isfinite(temp{ti}));
         temp{ti} = temp{ti}(tempInd);
     end
     slideDistance{i} = temp;
end
%%
dThetaRange = linspace(0,15,100);
dKhRange = linspace(-0.15, 0.05, 100);
durationRange = linspace(0, 0.1, 100);
dPhiRange = linspace(-5, 5, 100);
dKvRange = linspace(-0.05, 0.05, 100);
sdRange = linspace(0,2.5,100); 

histDTheta = zeros(6,length(dThetaRange)-1,2);
histDKh = zeros(6,length(dKhRange)-1,2);
histDuration = zeros(6,length(durationRange)-1,2);

histDPhi = zeros(6,length(dPhiRange)-1,2);
histDKv = zeros(6,length(dKvRange)-1,2);
histSD = zeros(6,length(sdRange)-1,2);
for mi = 1 : 6
    histDTheta(mi,:,1) = histcounts(dTheta{1}{mi}, dThetaRange, 'normalization', 'cdf');
    histDTheta(mi,:,2) = histcounts(dTheta{2}{mi}, dThetaRange, 'normalization', 'cdf');
    histDKh(mi,:,1) = histcounts(dKh{1}{mi}, dKhRange, 'normalization', 'cdf');
    histDKh(mi,:,2) = histcounts(dKh{2}{mi}, dKhRange, 'normalization', 'cdf');    
    histDuration(mi,:,1) = histcounts(duration{1}{mi}, durationRange, 'normalization', 'cdf');
    histDuration(mi,:,2) = histcounts(duration{2}{mi}, durationRange, 'normalization', 'cdf');
    
    histDPhi(mi,:,1) = histcounts(dPhi{1}{mi}, dPhiRange, 'normalization', 'cdf');
    histDPhi(mi,:,2) = histcounts(dPhi{2}{mi}, dPhiRange, 'normalization', 'cdf');
    histDKv(mi,:,1) = histcounts(dKv{1}{mi}, dKvRange, 'normalization', 'cdf');
    histDKv(mi,:,2) = histcounts(dKv{2}{mi}, dKvRange, 'normalization', 'cdf');
    histSD(mi,:,1) = histcounts(slideDistance{1}{mi}, sdRange, 'normalization', 'cdf');
    histSD(mi,:,2) = histcounts(slideDistance{2}{mi}, sdRange, 'normalization', 'cdf');
end
%%
figure, 
subplot(231), hold on
shadedErrorBar(dThetaRange(1:end-1), mean(squeeze(histDTheta(:,:,1))), std(squeeze(histDTheta(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dThetaRange(1:end-1), mean(squeeze(histDTheta(:,:,2))), std(squeeze(histDTheta(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max\Delta\theta')
ylabel('Cumulative proportion')

subplot(232), hold on
shadedErrorBar(dKhRange(1:end-1), mean(squeeze(histDKh(:,:,1))), std(squeeze(histDKh(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dKhRange(1:end-1), mean(squeeze(histDKh(:,:,2))), std(squeeze(histDKh(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max\Delta\kappa_H')

subplot(233), hold on
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,1))), std(squeeze(histDuration(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(durationRange(1:end-1), mean(squeeze(histDuration(:,:,2))), std(squeeze(histDuration(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(protraction touch duration)')
legend({'Matching naive (n = 6)', 'Expert (n = 6)'}, 'location', 'southeast')

subplot(234), hold on
shadedErrorBar(dPhiRange(1:end-1), mean(squeeze(histDPhi(:,:,1))), std(squeeze(histDPhi(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dPhiRange(1:end-1), mean(squeeze(histDPhi(:,:,2))), std(squeeze(histDPhi(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max\Delta\phi')
ylabel('Cumulative proportion')

subplot(235), hold on
shadedErrorBar(dKvRange(1:end-1), mean(squeeze(histDKv(:,:,1))), std(squeeze(histDKv(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(dKvRange(1:end-1), mean(squeeze(histDKv(:,:,2))), std(squeeze(histDKv(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max\Delta\kappa_V')

subplot(236), hold on
shadedErrorBar(sdRange(1:end-1), mean(squeeze(histSD(:,:,1))), std(squeeze(histSD(:,:,1)))/sqrt(6), 'lineprop', {'color', 'b'})
shadedErrorBar(sdRange(1:end-1), mean(squeeze(histSD(:,:,2))), std(squeeze(histSD(:,:,2)))/sqrt(6), 'lineprop', {'color', 'r'})
title('max(slide distance)')










%% Look at examples of each touch in different quadrants.

%% look into details of the examples of the quadrants of dphi vs dkv.
% from a single uber array - JK030 S03
load('D:\TPM\JK\suite2p\030\UberJK030S03_NC.mat')
tnums = cellfun(@(x) ones(1,length(x.protractionTouchChunksByWhisking)) * x.trialNum, u.trials, 'uniformoutput', false)';
dkv = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
dphi = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
slideDistance = cellfun(@(x) x.protractionTouchSlideDistanceByWhisking, u.trials, 'uniformoutput', false)';
firstLickTimes = cell(1,length(u.trials));
answerLickTimes = cell(1,length(u.trials));
touchOnsetTimes = cell(1,length(u.trials));
for i = 1 : length(u.trials)
    allLicks = union(union(u.trials{i}.leftLickTime, u.trials{i}.rightLickTime), u.trials{i}.answerLickTime);
    if isempty(allLicks)
        firstLickTimes{i} = u.trials{i}.poleMovingTime(end);
    else
        firstLickTimes{i} = allLicks(find(allLicks > u.trials{i}.poleUpTime(1), 1, 'first'));
    end
    
    if isempty(u.trials{i}.answerLickTime)
        answerLickTimes{i} = u.trials{i}.poleMovingTime(end);
    else
        answerLickTimes{i} = u.trials{i}.answerLickTime;
    end
    if isempty(u.trials{i}.protractionTouchChunksByWhisking)
        touchOnsetTimes{i} = [];
    else
        touchOnsetTimes{i} = cellfun(@(x) u.trials{i}.whiskerTime(x(1)), u.trials{i}.protractionTouchChunksByWhisking);
    end
end


figure, 
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')

%% pick one point, and get the value

[x,y] = ginput(1)
%%
[~, ind] = min((cell2mat(dkv)-x).^2 + (cell2mat(dphi)-y).^2);
tnumArray = cell2mat(tnums);
tnum = tnumArray(ind);
tnumInd = find(cellfun(@(x) ismember(tnum,x), tnums)); % index of u.trials
allTouchesBefore = sum(cellfun(@length, tnums(1:tnumInd-1)));
touchInd = ind - allTouchesBefore; % index of u.trials.protractionTouchChunksByWhisking

dkv{tnumInd}(touchInd)
dphi{tnumInd}(touchInd)

u.trialNums(tnumInd)
% slideDistance{tnumInd}(touchInd)



%% Divide by before and after the first licks

dkvBeforeFirstLick = cellfun(@(x,y,z) x(find(y<z)), dkv, touchOnsetTimes, firstLickTimes, 'un', 0);
dkvAfterFirstLick = cellfun(@(x,y,z) x(find(y>=z)), dkv, touchOnsetTimes, firstLickTimes, 'un', 0);
dphiBeforeFirstLick = cellfun(@(x,y,z) x(find(y<z)), dphi, touchOnsetTimes, firstLickTimes, 'un', 0);
dphiAfterFirstLick = cellfun(@(x,y,z) x(find(y>=z)), dphi, touchOnsetTimes, firstLickTimes, 'un', 0);

figure,
subplot(121), hold on
scatter(cell2mat(dkvBeforeFirstLick), cell2mat(dphiBeforeFirstLick), 'k.')
xlim([-0.1 0.1]), ylim([-15 15])
xlabel('max\Delta\kappa_V'), ylabel('max\Delta\phi')
title('Before the first lick')
axis square
subplot(122), hold on
scatter(cell2mat(dkvAfterFirstLick), cell2mat(dphiAfterFirstLick), 'k.')
xlim([-0.1 0.1]), ylim([-15 15])
title('After the first lick')
xlabel('max\Delta\kappa_V'), ylabel('max\Delta\phi')
axis square



%% among all the touches pick one
figure,
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')
xlabel('max\Delta\kappa_V'), ylabel('max\Delta\phi')
axis square

[x,y] = ginput(1);

[~, ind] = min((cell2mat(dkv)-x).^2 + (cell2mat(dphi)-y).^2);
tnumArray = cell2mat(tnums);
tnum = tnumArray(ind);
tnumInd = find(cellfun(@(x) ismember(tnum,x), tnums)); % index of u.trials
allTouchesBefore = sum(cellfun(@length, tnums(1:tnumInd-1)));
touchInd = ind - allTouchesBefore; % index of u.trials.protractionTouchChunksByWhisking

u.trialNums(tnumInd)

close all
utemp = u.trials{tnumInd};
figure('units', 'normalized', 'outerposition', [0 0 1 1])
subplot(211), hold on
plot(utemp.theta, 'b')
plot(utemp.phi, 'c')
yval = [min(union(utemp.theta, utemp.phi)), max(union(utemp.theta, utemp.phi))];
for ti = 1 : length(utemp.protractionTouchChunksByWhisking)
    tempX = utemp.protractionTouchChunksByWhisking{ti};
    if ti == touchInd
        patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [1 1 0], 'edgecolor', 'none')
    else
        patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [0.7 0.7 0.7], 'edgecolor', 'none')
    end
end
xticks([0:50:length(utemp.whiskerTime)])
plot(utemp.theta, 'b')
plot(utemp.phi, 'c')
xlim([0 length(utemp.whiskerTime)])
ylim(yval)
title(sprintf('%s %s trial #%d touch #%d (pole: %d)', u.mouseName, u.sessionName, utemp.trialNum, touchInd, utemp.angle))
legend({'\theta', '\phi'}, 'location', 'northwest')


subplot(212), hold on
plot(utemp.kappaH, 'r')
plot(utemp.kappaV, 'm')
yval = [min(union(utemp.kappaV, utemp.kappaH)), max(union(utemp.kappaV, utemp.kappaH))];
for ti = 1 : length(utemp.protractionTouchChunksByWhisking)
    tempX = utemp.protractionTouchChunksByWhisking{ti};    
    if ti == touchInd
        patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [1 1 0], 'edgecolor', 'none')
    else
        patch(([min(tempX), max(tempX), max(tempX), min(tempX)]), [yval(1) yval(1) yval(2) yval(2)], [0.7 0.7 0.7], 'edgecolor', 'none')
    end
end
plot(utemp.kappaH, 'r')
plot(utemp.kappaV, 'm')
xticks([0:50:length(utemp.whiskerTime)])
xlim([0 length(utemp.whiskerTime)])
ylim(yval)
legend({'\kappa_H', '\kappa_V'}, 'location', 'northwest')

figure, hold on
scatter(cell2mat(dkv), cell2mat(dphi), 'k.')
scatter(dkv{tnumInd}(touchInd), dphi{tnumInd}(touchInd), 'r', 'filled')
xlabel('max\Delta\kappa_V'), ylabel('max\Delta\phi')
axis square
title(sprintf('%s %s', u.mouseName, u.sessionName))



%% how many touches are there in each quadrant? (from all touches, all mice)
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

dkvCell = cell(1,length(mice));
dphiCell = cell(1,length(mice));
durationCell = cell(1,length(mice));

for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(1);
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session),'u')
    
    dkvCell{mi} = cellfun(@(x) x.protractionTouchDKappaVByWhisking, u.trials, 'uniformoutput', false)';
    dphiCell{mi} = cellfun(@(x) x.protractionTouchDPhiByWhisking, u.trials, 'uniformoutput', false)';
    durationCell{mi} = cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials, 'uniformoutput', false)';
end

dkv = cell2mat(cellfun(@(x) cell2mat(x), dkvCell, 'uniformoutput', false));
dphi = cell2mat(cellfun(@(x) cell2mat(x), dphiCell, 'uniformoutput', false));

length(intersect(find(dkv>0), find(dphi>0)))/length(dkv)
length(intersect(find(dkv<0), find(dphi>0)))/length(dkv)
length(intersect(find(dkv<0), find(dphi<0)))/length(dkv)
length(intersect(find(dkv>0), find(dphi<0)))/length(dkv)