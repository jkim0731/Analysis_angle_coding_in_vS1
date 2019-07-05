mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[1:25],[1:6,15:17,23:25],[1:25],[1:21], [1:24], [1:31], [1:28], [1:30], [1:21,26], [1:21], [1:30], [1:13]};
sessions = {[1:18],[1:7,15],[1:20],[1:16], [1:24], [1:31], [1:21], [1:30], [1:20], [1:21], [1:30], [1:13]};

pcbysessions = cell(length(mice),1);
databytrials = cell(length(mice),1);
pcbytrials = cell(length(mice),1);
sessionDurations = cell(length(mice),1);
baseD = 'C:\Data\SoloData\';

for mi = 1 : length(mice)
    mouse = sprintf('JK%03d',mice(mi));
    cd([baseD, mouse])
    fn = sprintf('behavior_%s.mat',mouse);
    load(fn)   
    
    pcbysessions{mi} = zeros(length(sessions{mi}),1);
    databytrials{mi} = [];
    sessionDurations{mi} = zeros(length(sessions{mi}),1);
    for si = 1 : length(sessions{mi})
        bi = find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',sessions{mi}(si))), b));
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{bi}.trials));
        pcbysessions{mi}(si) = sum(b{bi}.hitTrialInds(notooind))/(sum(b{bi}.hitTrialInds(notooind)) + sum(b{bi}.faTrialInds(notooind))) * 100;
        tempCorrects = b{bi}.trialCorrects;
        tempCorrects(tempCorrects==-1) = []; % remove misses
        databytrials{mi} = [databytrials{mi}, tempCorrects];
        sessionDurations{mi}(si) = length(tempCorrects);        
    end
    pcbytrials{mi} = zeros(length(databytrials{mi})-99,1);
    for ti = 1 : length(databytrials{mi})-99
        pcbytrials{mi}(ti) = sum(databytrials{mi}(ti:ti+99));
    end
end
%%
figure, hold all,
for mi = 1 : length(mice)
    if length(find(pcbysessions{mi} > 75 )) > 2
        plot(pcbysessions{mi}, 'k-', 'linewidth', 5)
    else
        plot(pcbysessions{mi}, 'color',[0.7 0.7 0.7], 'linewidth', 3)
    end
end
plot(0:35, ones(36,1)*75, 'r--', 'linewidth', 3)
xlabel('Sessions'), ylabel('% Correct')
set(gca,'fontsize',25)


%%
learned = [25,27,30,36,39,52];
discreteCorrect = zeros(length(learned),5);
for mi = 1 : length(learned)
    mouse = sprintf('JK%03d',learned(mi));
    cd([baseD, mouse])
    fn = sprintf('behavior_%s.mat',mouse);
    load(fn)
    binds = find(cellfun(@(x) strcmpi(x.taskTarget, 'Angle-Discrete'), b));    
    notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(1)}.trials));
    discreteCorrect(mi,1) = sum(b{binds(1)}.hitTrialInds(notooind))/(sum(b{binds(1)}.hitTrialInds(notooind)) + sum(b{binds(1)}.faTrialInds(notooind))) * 100;
    notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(2)}.trials));
    discreteCorrect(mi,3) = sum(b{binds(2)}.hitTrialInds(notooind))/(sum(b{binds(2)}.hitTrialInds(notooind)) + sum(b{binds(2)}.faTrialInds(notooind))) * 100;
    binds = find(cellfun(@(x) strcmpi(x.distractor, 'Discrete'), b));        
    if learned(mi) == 52
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(2)}.trials));
        discreteCorrect(mi,4) = sum(b{binds(2)}.hitTrialInds(notooind))/(sum(b{binds(2)}.hitTrialInds(notooind)) + sum(b{binds(2)}.faTrialInds(notooind))) * 100;
    else
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(1)}.trials));
        discreteCorrect(mi,4) = sum(b{binds(1)}.hitTrialInds(notooind))/(sum(b{binds(1)}.hitTrialInds(notooind)) + sum(b{binds(1)}.faTrialInds(notooind))) * 100;
    end
    if mi == 3
        discreteCorrect(mi,5) = sum(b{end-2}.hitTrialInds)/(sum(b{end-2}.hitTrialInds) + sum(b{end-2}.faTrialInds)) * 100;
    elseif mi == 6
        discreteCorrect(mi,5) = NaN;
    else
        discreteCorrect(mi,5) = sum(b{end-1}.hitTrialInds)/(sum(b{end-1}.hitTrialInds) + sum(b{end-1}.faTrialInds)) * 100;
    end
    mind = find(mice == learned(mi));
    tempPC = sort(pcbysessions{mind}, 'descend');
    discreteCorrect(mi,2) = mean(tempPC(1:3));
%     discreteCorrect(mi,5) = mean(tempPC(1:3))/4 + 50*3/4;
%     discreteCorrect(mi,6) = mean(tempPC(1:3))*3/4 + 50/4;
end
%%
figure, hold all
for mi = 1 : length(learned)
    plot(pcbysessions{find(mice==learned(mi))}, 'linewidth', 3)
end
xlabel('Sessions'), ylabel('% Correct')
legendmice = cell((length(learned)),1);
for i = 1 : length(learned)
    legendmice{i} = sprintf('JK%03d',learned(i));
end
legend(legendmice)

%%
figure, hold all
for mi = 1 : length(learned)
    plot(discreteCorrect(mi,1:5), 'k.', 'markersize',30)    
end

xticks([1:5]), xlim([0.5 5.5]), ylabel('% Correct')
labels = ({'Before Learning', 'Expert', 'After Learning','Radial Distractor', 'No Whisker'});
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
XTickLabel = labels;
set(gca,'fontsize',20, 'xticklabel', labels)
%% 
angles = 45:15:135;
learned = [25,27,30,36,39,52];
plickLBefore = zeros(length(learned),length(angles));
plickLAfter = zeros(length(learned),length(angles));
for mi = 1 : length(learned)
% for mi = 5
    mouse = sprintf('JK%03d',learned(mi));
    cd([baseD, mouse])
    fn = sprintf('behavior_%s.mat',mouse);
    load(fn)
    binds = find(cellfun(@(x) strcmpi(x.taskTarget, 'Angle-Discrete'), b));
        
    lefts = find(cellfun(@(x) strcmp(x.trialType(1), 'l'), b{binds(1)}.trials));
    rights = find(cellfun(@(x) strcmp(x.trialType(1), 'r'), b{binds(1)}.trials));
    lickLonL = intersect(lefts, b{binds(1)}.hitTrialNums);
    lickLonR = intersect(rights, b{binds(1)}.faTrialNums);
    lickL = union(lickLonL, lickLonR);
    lickRonL = intersect(lefts, b{binds(1)}.faTrialNums);
    lickRonR = intersect(rights, b{binds(1)}.hitTrialNums);
    lickR = union(lickRonL, lickRonR);    
    for ai = 1 : length(angles)
        angleTrials = find(cellfun(@(x) x.servoAngle == angles(ai), b{binds(1)}.trials));
        plickLBefore(mi,ai) = length(intersect(lickL,angleTrials))/(length(intersect(lickL,angleTrials)) + length(intersect(lickR,angleTrials)))*100;
    end
    
    lefts = find(cellfun(@(x) strcmp(x.trialType(1), 'l'), b{binds(end)}.trials));
    rights = find(cellfun(@(x) strcmp(x.trialType(1), 'r'), b{binds(end)}.trials));
    lickLonL = intersect(lefts, b{binds(end)}.hitTrialNums);
    lickLonR = intersect(rights, b{binds(end)}.faTrialNums);
    lickL = union(lickLonL, lickLonR);
    lickRonL = intersect(lefts, b{binds(end)}.faTrialNums);
    lickRonR = intersect(rights, b{binds(end)}.hitTrialNums);
    lickR = union(lickRonL, lickRonR);    
    for ai = 1 : length(angles)
        angleTrials = find(cellfun(@(x) x.servoAngle == angles(ai), b{binds(end)}.trials));
        plickLAfter(mi,ai) = length(intersect(lickL,angleTrials))/(length(intersect(lickL,angleTrials)) + length(intersect(lickR,angleTrials)))*100;
    end
end
%%
figure,
subplot(121), hold on
for mi = 1 : length(learned)    
    plot(angles, plickLBefore(mi,:), '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
end
errorbar(angles, mean(plickLBefore),sqrt(std(plickLBefore)), 'r-', 'linewidth', 5)
ylabel('Left lick percentage (%)'), xlabel('Angle (\circ)'), xticks(angles), title('Before learning'), ylim([0 100])
set(gca,'fontsize',15)

subplot(122), hold on
for mi = 1 : length(learned)    
    plot(angles, plickLAfter(mi,:), '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
end
errorbar(angles, mean(plickLAfter),sqrt(std(plickLAfter)), 'r-', 'linewidth', 5)
ylabel('Left lick percentage (%)'), xlabel('Angle (\circ)'), xticks(angles), title('After learning'), ylim([0 100])
set(gca,'fontsize',15)

%%
figure, plot(angles, plickLBefore(5,:))
%%

tn = 7;
wla.trials{tn}.trialNum
wla.trials{tn}.thTouchFrames
wla.trials{tn}.servoAngle

inds = setdiff(1:wla.trials{tn}.nof, wla.trials{tn}.thTouchFrames);
theta = wla.trials{tn}.theta{1};
Ftheta = wla.trials{tn}.theta{2};
kappa = wla.trials{tn}.kappa{1};
% figure, plot3(theta(inds), Ftheta(inds), kappa(inds), 'k.'), hold on,
% plot3(theta(wla.trials{tn}.thTouchFrames), Ftheta(wla.trials{tn}.thTouchFrames),kappa(wla.trials{tn}.thTouchFrames), 'b.')


p3d = [wla.trials{tn}.whiskerEdgeCoord,wla.trials{tn}.apPosition];
p3d = p3d';
p4d = [p3d; ones(1,size(p3d,2))];
A = viewmtx(psi1((wla.trials{tn}.servoAngle-45)/15+1), 90-psi2((wla.trials{tn}.servoAngle-45)/15+1));
projpoints = A*p4d;
% projpoints = round(projpoints*10000)/10000;
% projpoints = unique(projpoints(1:2,:)','rows');

figure, hold all
plot(projpoints(1,:), projpoints(2,:),'k.'), plot(wla.trials{tn}.thPolygon(1,:), wla.trials{tn}.thPolygon(2,:),'r.')
plot(projpoints(1,wla.trials{tn}.thTouchFrames), projpoints(2,wla.trials{tn}.thTouchFrames),'b.'), 

%%
find(theta > 5 & kappa < 0)

%%
thetaInterp = interp1(find(isfinite(theta)), theta(isfinite(theta)), find(isnan(theta)));
figure, plot(1:length(theta), theta, 'k.'), hold on
plot(find(isnan(theta)), thetaInterp, 'r.')
thetaComb = theta;
thetaComb(isnan(theta)) = thetaInterp;
thetaS = smooth(thetaComb);
plot(thetaS, 'b-')
%%
thetaRS = smoothdata(theta, 'movmean',5, 'omitNan');

%%
figure, hold on, plot(thetaRS,'g-'), plot(1:length(theta), theta, 'k.')
%%
kappaInterp = interp1(find(isfinite(kappa)), kappa(isfinite(kappa)), find(isnan(kappa)));
figure, plot(1:length(kappa), kappa, 'k.'), hold on
plot(find(isnan(kappa)), kappaInterp, 'r.')





