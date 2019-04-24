clear
baseDir = 'C:\Data\suite2p\';
cd(baseDir)
load('cellFunctionRidgeDE010.mat')
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];
L4Mi = 15:16;

whiskingInActiveNaive = zeros(length(naiveMi),1);
whiskingInTouchNaive = zeros(length(naiveMi),1);
whiskingInTuningNaive = zeros(length(naiveMi),1);
for mi = 1 : length(naiveMi)
    mouse = mice(naiveMi(mi));
    session = sessions{naiveMi(mi)}(1);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    whiskingInActiveNaive(mi) = length(naive(mi).whiskingID)/length(naive(mi).cellNums);
    whiskingInTouchNaive(mi) = length(intersect(naive(mi).whiskingID, naive(mi).touchID)) / length(naive(mi).touchID);
    whiskingInTuningNaive(mi) = length(intersect(naive(mi).whiskingID, spk.touchID(find(spk.tuned)))) / length(find(spk.tuned));
end
%%
whiskingInActiveExpert = zeros(length(expertMi),1);
whiskingInTouchExpert = zeros(length(expertMi),1);
whiskingInTuningExpert = zeros(length(expertMi),1);
for mi = 1 : length(expertMi)
    mouse = mice(expertMi(mi));
    session = sessions{expertMi(mi)}(2);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    whiskingInActiveExpert(mi) = length(expert(mi).whiskingID)/length(expert(mi).cellNums);
    whiskingInTouchExpert(mi) = length(intersect(expert(mi).whiskingID, expert(mi).touchID)) / length(expert(mi).touchID);
    whiskingInTuningExpert(mi) = length(intersect(expert(mi).whiskingID, spk.touchID(find(spk.tuned)))) / length(find(spk.tuned));
end

%%
whiskingInActiveL4 = zeros(length(expertMi),1);
whiskingInTouchL4 = zeros(length(expertMi),1);
whiskingInTuningL4 = zeros(length(expertMi),1);
for mi = 1 : length(expertMi)
    mouse = mice(L4Mi(mi));
    session = sessions{L4Mi(mi)}(1);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    whiskingInActiveL4(mi) = length(L4(mi).whiskingID)/length(L4(mi).cellNums);
    whiskingInTouchL4(mi) = length(intersect(L4(mi).whiskingID, L4(mi).touchID)) / length(L4(mi).touchID);
    whiskingInTuningL4(mi) = length(intersect(L4(mi).whiskingID, spk.touchID(find(spk.tuned)))) / length(find(spk.tuned));
end

%% whisking in total, touch, and tuning in naive

figure, hold on
bar(0.7,mean(whiskingInActiveNaive), 0.3, 'w', 'linewidth', 3)
bar(1.0,mean(whiskingInTouchNaive), 0.3, 'w', 'linewidth', 3)
bar(1.3,mean(whiskingInTuningNaive), 0.3, 'w', 'linewidth', 3)
errorbar(0.7, mean(whiskingInActiveNaive), std(whiskingInActiveNaive)/sqrt(length(whiskingInActiveNaive)), 'k.', 'linewidth', 3)
errorbar(1.0, mean(whiskingInTouchNaive), std(whiskingInTouchNaive)/sqrt(length(whiskingInTouchNaive)), 'k.', 'linewidth', 3)
errorbar(1.3, mean(whiskingInTuningNaive), std(whiskingInTuningNaive)/sqrt(length(whiskingInTuningNaive)), 'k.', 'linewidth', 3)
ylabel('Proportion of whisking cells')
xticks([0.7 1.0 1.3]), xlim([0.3 1.7])
xticklabels({'Active', 'Touch', 'Angle-tuned'})
xtickangle(45)
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

% in expert
% bar(1.7,mean(whiskingInTotalExpert), 0.28, 'k')
% bar(2.0,mean(whiskingInTouchExpert), 0.28, 'k')
% bar(2.3,mean(whiskingInTuningExpert), 0.28, 'k')
% errorbar(1.7, mean(whiskingInTotalExpert), std(whiskingInTotalExpert)/sqrt(length(whiskingInTotalExpert)), 'k.')
% errorbar(2.0, mean(whiskingInTouchExpert), std(whiskingInTouchExpert)/sqrt(length(whiskingInTouchExpert)), 'k.')
% errorbar(2.3, mean(whiskingInTuningExpert), std(whiskingInTuningExpert)/sqrt(length(whiskingInTuningExpert)), 'k.')
% 


%% L2/3 C2, L2/3 non-C2, L4 C2, and L4 non-C2
whiskingInActiveNaive = zeros(length(naiveMi),4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskingInTouchNaive = zeros(length(naiveMi),4);
whiskingInTuningNaive = zeros(length(naiveMi),4);
for mi = 1 : length(naiveMi)
    mouse = mice(naiveMi(mi));
    session = sessions{naiveMi(mi)}(1);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    
    L23cells = naive(mi).cellNums(find(naive(mi).cellDepths < 350));
    C2cells = naive(mi).cellNums(find(naive(mi).isC2));
    L4cells = naive(mi).cellNums(find(naive(mi).cellDepths >= 350));
    nonC2cells = naive(mi).cellNums(find(naive(mi).isC2==0));
    tunedcells = spk.touchID(find(spk.tuned));
    
    whiskingInActiveNaive(mi,1) = length(intersect(naive(mi).whiskingID, intersect(L23cells, C2cells))) / length(intersect(L23cells, C2cells));
    whiskingInActiveNaive(mi,2) = length(intersect(naive(mi).whiskingID, intersect(L23cells, nonC2cells))) / length(intersect(L23cells, nonC2cells));
    whiskingInActiveNaive(mi,3) = length(intersect(naive(mi).whiskingID, intersect(L4cells, C2cells))) / length(intersect(L4cells, C2cells));
    whiskingInActiveNaive(mi,4) = length(intersect(naive(mi).whiskingID, intersect(L4cells, nonC2cells))) / length(intersect(L4cells, nonC2cells));
    
    whiskingInTouchNaive(mi,1) = length(intersect(naive(mi).whiskingID, intersect(naive(mi).touchID, intersect(L23cells, C2cells)) )) / length(intersect(naive(mi).touchID, intersect(L23cells, C2cells)));
    whiskingInTouchNaive(mi,2) = length(intersect(naive(mi).whiskingID, intersect(naive(mi).touchID, intersect(L23cells, nonC2cells)) )) / length(intersect(naive(mi).touchID, intersect(L23cells, nonC2cells)));
    whiskingInTouchNaive(mi,3) = length(intersect(naive(mi).whiskingID, intersect(naive(mi).touchID, intersect(L4cells, C2cells)) )) / length(intersect(naive(mi).touchID, intersect(L4cells, C2cells)));
    whiskingInTouchNaive(mi,4) = length(intersect(naive(mi).whiskingID, intersect(naive(mi).touchID, intersect(L4cells, nonC2cells)) )) / length(intersect(naive(mi).touchID, intersect(L4cells, nonC2cells)));
    
    whiskingInTuningNaive(mi,1) = length(intersect(naive(mi).whiskingID, intersect(tunedcells, intersect(L23cells, C2cells)) )) / length(intersect(tunedcells, intersect(L23cells, C2cells)));
    whiskingInTuningNaive(mi,2) = length(intersect(naive(mi).whiskingID, intersect(tunedcells, intersect(L23cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L23cells, nonC2cells)));
    whiskingInTuningNaive(mi,3) = length(intersect(naive(mi).whiskingID, intersect(tunedcells, intersect(L4cells, C2cells)) )) / length(intersect(tunedcells, intersect(L4cells, C2cells)));
    whiskingInTuningNaive(mi,4) = length(intersect(naive(mi).whiskingID, intersect(tunedcells, intersect(L4cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L4cells, nonC2cells)));    
end

figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(0.7,mean(whiskingInActiveNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.0,mean(whiskingInTouchNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.3,mean(whiskingInTuningNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    errorbar(0.7, mean(whiskingInActiveNaive(:,i)), std(whiskingInActiveNaive(:,i))/sqrt(length(whiskingInActiveNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.0, mean(whiskingInTouchNaive(:,i)), std(whiskingInTouchNaive(:,i))/sqrt(length(whiskingInTouchNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.3, mean(whiskingInTuningNaive(:,i)), std(whiskingInTuningNaive(:,i))/sqrt(length(whiskingInTuningNaive(:,i))), 'k.', 'linewidth', 3)
    ylabel('Proportion of whisking cells')
    xticks([0.7 1.0 1.3]), xlim([0.3 1.7])
    xticklabels({'Active', 'Touch', 'Angle-tuned'})
    xtickangle(45)
    switch i
        case 1
            title('L2/3 C2')
        case 2
            title('L2/3 non-C2')
        case 3
            title('L4 C2')
        case 4
            title('L4 non-C2')
    end
    set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')    
end

%% L2/3 C2, L2/3 non-C2, L4 C2, and L4 non-C2
%% in experts
whiskingInActiveExpert = zeros(length(expertMi),4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskingInTouchExpert = zeros(length(expertMi),4);
whiskingInTuningExpert = zeros(length(expertMi),4);
for mi = 1 : length(expertMi)
    mouse = mice(expertMi(mi));
    session = sessions{expertMi(mi)}(2);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    
    L23cells = naive(mi).cellNums(find(expert(mi).cellDepths < 350));
    C2cells = naive(mi).cellNums(find(expert(mi).isC2));
    L4cells = naive(mi).cellNums(find(expert(mi).cellDepths >= 350));
    nonC2cells = naive(mi).cellNums(find(expert(mi).isC2==0));
    tunedcells = spk.touchID(find(spk.tuned));
    
    whiskingInActiveExpert(mi,1) = length(intersect(expert(mi).whiskingID, intersect(L23cells, C2cells))) / length(intersect(L23cells, C2cells));
    whiskingInActiveExpert(mi,2) = length(intersect(expert(mi).whiskingID, intersect(L23cells, nonC2cells))) / length(intersect(L23cells, nonC2cells));
    whiskingInActiveExpert(mi,3) = length(intersect(expert(mi).whiskingID, intersect(L4cells, C2cells))) / length(intersect(L4cells, C2cells));
    whiskingInActiveExpert(mi,4) = length(intersect(expert(mi).whiskingID, intersect(L4cells, nonC2cells))) / length(intersect(L4cells, nonC2cells));
    
    whiskingInTouchExpert(mi,1) = length(intersect(expert(mi).whiskingID, intersect(naive(mi).touchID, intersect(L23cells, C2cells)) )) / length(intersect(expert(mi).touchID, intersect(L23cells, C2cells)));
    whiskingInTouchExpert(mi,2) = length(intersect(expert(mi).whiskingID, intersect(naive(mi).touchID, intersect(L23cells, nonC2cells)) )) / length(intersect(expert(mi).touchID, intersect(L23cells, nonC2cells)));
    whiskingInTouchExpert(mi,3) = length(intersect(expert(mi).whiskingID, intersect(naive(mi).touchID, intersect(L4cells, C2cells)) )) / length(intersect(expert(mi).touchID, intersect(L4cells, C2cells)));
    whiskingInTouchExpert(mi,4) = length(intersect(expert(mi).whiskingID, intersect(naive(mi).touchID, intersect(L4cells, nonC2cells)) )) / length(intersect(expert(mi).touchID, intersect(L4cells, nonC2cells)));
    
    whiskingInTuningExpert(mi,1) = length(intersect(expert(mi).whiskingID, intersect(tunedcells, intersect(L23cells, C2cells)) )) / length(intersect(tunedcells, intersect(L23cells, C2cells)));
    whiskingInTuningExpert(mi,2) = length(intersect(expert(mi).whiskingID, intersect(tunedcells, intersect(L23cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L23cells, nonC2cells)));
    whiskingInTuningExpert(mi,3) = length(intersect(expert(mi).whiskingID, intersect(tunedcells, intersect(L4cells, C2cells)) )) / length(intersect(tunedcells, intersect(L4cells, C2cells)));
    whiskingInTuningExpert(mi,4) = length(intersect(expert(mi).whiskingID, intersect(tunedcells, intersect(L4cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L4cells, nonC2cells)));    
end

figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(0.7,mean(whiskingInActiveExpert(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.0,mean(whiskingInTouchExpert(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.3,mean(whiskingInTuningExpert(:,i)), 0.3, 'w', 'linewidth', 3)
    errorbar(0.7, mean(whiskingInActiveExpert(:,i)), std(whiskingInActiveExpert(:,i))/sqrt(length(whiskingInActiveExpert(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.0, mean(whiskingInTouchExpert(:,i)), std(whiskingInTouchExpert(:,i))/sqrt(length(whiskingInTouchExpert(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.3, mean(whiskingInTuningExpert(:,i)), std(whiskingInTuningExpert(:,i))/sqrt(length(whiskingInTuningExpert(:,i))), 'k.', 'linewidth', 3)
    ylabel('Proportion of whisking cells')
    xticks([0.7 1.0 1.3]), xlim([0.3 1.7])
    xticklabels({'Active', 'Touch', 'Angle-tuned'})
    xtickangle(45)
    switch i
        case 1
            title('L2/3 C2')
        case 2
            title('L2/3 non-C2')
        case 3
            title('L4 C2')
        case 4
            title('L4 non-C2')
    end
    set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')    
end


%% L2/3 C2, L2/3 non-C2, L4 C2, and L4 non-C2
%% in matching naive
whiskingInActiveMatchingNaive = zeros(length(expertMi),4); % 1: L2/3 C2, 2: L2/3 non-C2, 3: L4 C2, 4: L4 non-C2
whiskingInTouchMatchingNaive = zeros(length(expertMi),4);
whiskingInTuningMatchingNaive = zeros(length(expertMi),4);
for mi = 1 : length(expertMi)
    mouse = mice(expertMi(mi));
    matchi = expertMi(mi);
    session = sessions{expertMi(mi)}(1);
    cd(sprintf('%s%03d', baseDir, mouse))
    load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
    
    L23cells = naive(matchi).cellNums(find(naive(matchi).cellDepths < 350));
    C2cells = naive(matchi).cellNums(find(naive(matchi).isC2));
    L4cells = naive(matchi).cellNums(find(naive(matchi).cellDepths >= 350));
    nonC2cells = naive(matchi).cellNums(find(naive(matchi).isC2==0));
    tunedcells = spk.touchID(find(spk.tuned));
    
    whiskingInActiveMatchingNaive(mi,1) = length(intersect(naive(matchi).whiskingID, intersect(L23cells, C2cells))) / length(intersect(L23cells, C2cells));
    whiskingInActiveMatchingNaive(mi,2) = length(intersect(naive(matchi).whiskingID, intersect(L23cells, nonC2cells))) / length(intersect(L23cells, nonC2cells));
    whiskingInActiveMatchingNaive(mi,3) = length(intersect(naive(matchi).whiskingID, intersect(L4cells, C2cells))) / length(intersect(L4cells, C2cells));
    whiskingInActiveMatchingNaive(mi,4) = length(intersect(naive(matchi).whiskingID, intersect(L4cells, nonC2cells))) / length(intersect(L4cells, nonC2cells));
    
    whiskingInTouchMatchingNaive(mi,1) = length(intersect(naive(matchi).whiskingID, intersect(naive(matchi).touchID, intersect(L23cells, C2cells)) )) / length(intersect(naive(matchi).touchID, intersect(L23cells, C2cells)));
    whiskingInTouchMatchingNaive(mi,2) = length(intersect(naive(matchi).whiskingID, intersect(naive(matchi).touchID, intersect(L23cells, nonC2cells)) )) / length(intersect(naive(matchi).touchID, intersect(L23cells, nonC2cells)));
    whiskingInTouchMatchingNaive(mi,3) = length(intersect(naive(matchi).whiskingID, intersect(naive(matchi).touchID, intersect(L4cells, C2cells)) )) / length(intersect(naive(matchi).touchID, intersect(L4cells, C2cells)));
    whiskingInTouchMatchingNaive(mi,4) = length(intersect(naive(matchi).whiskingID, intersect(naive(matchi).touchID, intersect(L4cells, nonC2cells)) )) / length(intersect(naive(matchi).touchID, intersect(L4cells, nonC2cells)));
    
    whiskingInTuningMatchingNaive(mi,1) = length(intersect(naive(matchi).whiskingID, intersect(tunedcells, intersect(L23cells, C2cells)) )) / length(intersect(tunedcells, intersect(L23cells, C2cells)));
    whiskingInTuningMatchingNaive(mi,2) = length(intersect(naive(matchi).whiskingID, intersect(tunedcells, intersect(L23cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L23cells, nonC2cells)));
    whiskingInTuningMatchingNaive(mi,3) = length(intersect(naive(matchi).whiskingID, intersect(tunedcells, intersect(L4cells, C2cells)) )) / length(intersect(tunedcells, intersect(L4cells, C2cells)));
    whiskingInTuningMatchingNaive(mi,4) = length(intersect(naive(matchi).whiskingID, intersect(tunedcells, intersect(L4cells, nonC2cells)) )) / length(intersect(tunedcells, intersect(L4cells, nonC2cells)));    
end

figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(0.7,mean(whiskingInActiveMatchingNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.0,mean(whiskingInTouchMatchingNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    bar(1.3,mean(whiskingInTuningMatchingNaive(:,i)), 0.3, 'w', 'linewidth', 3)
    errorbar(0.7, mean(whiskingInActiveMatchingNaive(:,i)), std(whiskingInActiveMatchingNaive(:,i))/sqrt(length(whiskingInActiveMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.0, mean(whiskingInTouchMatchingNaive(:,i)), std(whiskingInTouchMatchingNaive(:,i))/sqrt(length(whiskingInTouchMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.3, mean(whiskingInTuningMatchingNaive(:,i)), std(whiskingInTuningMatchingNaive(:,i))/sqrt(length(whiskingInTuningMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    ylabel('Proportion of whisking cells')
    xticks([0.7 1.0 1.3]), xlim([0.3 1.7])
    xticklabels({'Active', 'Touch', 'Angle-tuned'})
    xtickangle(45)
    switch i
        case 1
            title('L2/3 C2')
        case 2
            title('L2/3 non-C2')
        case 3
            title('L4 C2')
        case 4
            title('L4 non-C2')
    end
    set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')    
end

%%

figure, 
for i = 1 : 4
    subplot(2,2,i), hold on
    bar(0.8,mean(whiskingInActiveMatchingNaive(:,i)), 0.4, 'w', 'linewidth', 3)
    bar(1.2,mean(whiskingInActiveExpert(:,i)), 0.4, 'k', 'linewidth', 3)
    
    bar(1.8,mean(whiskingInTouchMatchingNaive(:,i)), 0.4, 'w', 'linewidth', 3)
    bar(2.2,mean(whiskingInTouchExpert(:,i)), 0.4, 'k', 'linewidth', 3)
    
    bar(2.8,mean(whiskingInTuningMatchingNaive(:,i)), 0.4, 'w', 'linewidth', 3)
    bar(3.2,mean(whiskingInTuningExpert(:,i)), 0.4, 'k', 'linewidth', 3)
    
    errorbar(0.8, mean(whiskingInActiveMatchingNaive(:,i)), [], std(whiskingInActiveMatchingNaive(:,i))/sqrt(length(whiskingInActiveMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(1.2, mean(whiskingInActiveExpert(:,i)), [], std(whiskingInActiveExpert(:,i))/sqrt(length(whiskingInActiveExpert(:,i))), 'k.', 'linewidth', 3)
    
    errorbar(1.8, mean(whiskingInTouchMatchingNaive(:,i)), [], std(whiskingInTouchMatchingNaive(:,i))/sqrt(length(whiskingInTouchMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(2.2, mean(whiskingInTouchExpert(:,i)), [], std(whiskingInTouchExpert(:,i))/sqrt(length(whiskingInTouchExpert(:,i))), 'k.', 'linewidth', 3)
    
    errorbar(2.8, mean(whiskingInTuningMatchingNaive(:,i)), [], std(whiskingInTuningMatchingNaive(:,i))/sqrt(length(whiskingInTuningMatchingNaive(:,i))), 'k.', 'linewidth', 3)
    errorbar(3.2, mean(whiskingInTuningExpert(:,i)), [], std(whiskingInTuningExpert(:,i))/sqrt(length(whiskingInTuningExpert(:,i))), 'k.', 'linewidth', 3)
    
    switch i
        case 1
            title('L2/3 C2')
            ylabel('Proportion of whisking cells')
            xticks([]), xlim([0.4 3.6])    

        case 2
            title('L2/3 non-C2')
            xticks([]), xlim([0.4 3.6])
            legend({'Naive', 'Learned'}, 'box', 'off', 'location', 'northwest')

        case 3
            title('L4 C2')
            ylabel('Proportion of whisking cells')
            xticks([1 2 3]), xlim([0.4 3.6])
            xticklabels({'Active', 'Touch', 'Angle-tuned'}), xtickangle(45)

        case 4
            title('L4 non-C2')            
            xticks([1 2 3]), xlim([0.4 3.6])
            xticklabels({'Active', 'Touch', 'Angle-tuned'}), xtickangle(45)
            
    end
    set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')
end

%% There is increase in whisking cells among touch and tuned cells AFTER learning
%% Is this because of larger variance in whisking after learning?
%% Check whisking parameter distribution before and after learning
%% theta (or onset), amplitude, midpoint
%% their variances

%% (1) Using predictors
baseDir = 'C:\Data\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};
naiveMi = 1:12;
expertMi = [1,2,3,4,7,9];
nonlearnerMi = setdiff(naiveMi, expertMi);

onsetNotLearnedTpmtime = cell(length(nonlearnerMi), 1);
onsetBeforeLearningTpmtime = cell(length(expertMi), 1);
onsetAfterLearningTpmtime = cell(length(expertMi), 1);

amplitudeNotLearnedTpmtime = cell(length(nonlearnerMi), 1);
amplitudeBeforeLearningTpmtime = cell(length(expertMi), 1);
amplitudeAfterLearningTpmtime = cell(length(expertMi), 1);

midpointNotLearnedTpmtime = cell(length(nonlearnerMi), 1);
midpointBeforeLearningTpmtime = cell(length(expertMi), 1);
midpointAfterLearningTpmtime = cell(length(expertMi), 1);

thetaNotLearned = cell(length(nonlearnerMi), 1);
thetaBeforeLearning = cell(length(expertMi), 1);
thetaAfterLearning = cell(length(expertMi), 1);

amplitudeNotLearned = cell(length(nonlearnerMi), 1);
amplitudeBeforeLearning = cell(length(expertMi), 1);
amplitudeAfterLearning = cell(length(expertMi), 1);

midpointNotLearned = cell(length(nonlearnerMi), 1);
midpointBeforeLearning = cell(length(expertMi), 1);
midpointAfterLearning = cell(length(expertMi), 1);

for i = 1 : length(nonlearnerMi)
    mouse = mice(nonlearnerMi(i));
    session = sessions{nonlearnerMi(i)};
    cd(sprintf('%s%03d',baseDir, mouse))
      
    ufn = sprintf('UberJK%03dS%02d',mouse,session);
    load(ufn)
    theta = cell2mat(cellfun(@(x) x.theta, u.trials, 'uniformoutput', false));
    [~, ~, thetaS, ~, amplitudeS, midpointS] = jkWhiskerDecomposition(theta);
    
    thetaNotLearned{i} = theta;
    amplitudeNotLearned{i} = amplitudeS;
    midpointNotLearned{i} = midpointS;
    
    whiskingOnsetCell = cell(1, length(u.trials));
    whiskingAmplitudeCell = cell(1, length(u.trials));
    whiskingMidpointCell = cell(1, length(u.trials));
    for ti = 1 : length(u.trials)
        currTrial = u.trials{ti};
        time = [0, currTrial.tpmTime{1}];
        wtimes = currTrial.whiskerTime;
        if iscell(wtimes)
            wtimes = wtimes{1};
        end
        [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
        whiskerVideoFrameDuration = u.trials{ti}.frameDuration; % in s
        onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
        whiskingOnsetCell{ti} = histcounts(onsetTimes, time);

        tempMid = zeros(1,length(time)-1);
        tempAmp = zeros(1,length(time)-1);
        for timei = 1 : length(tempMid)
            startInd = find(wtimes >= time(timei), 1, 'first');
            endInd = find(wtimes < time(timei+1), 1, 'last');
            tempMid(timei) = mean(midpoint(startInd:endInd));
            tempAmp(timei) = max(amplitude(startInd:endInd));
        end
        tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
        tempAmp(isnan(tempAmp)) = deal(mode(tempAmp(isfinite(tempAmp))));
        whiskingMidpointCell{ti} = tempMid;
        whiskingAmplitudeCell{ti} = tempAmp;
    end
    onsetNotLearnedTpmtime{i} = cell2mat(whiskingOnsetCell);
    amplitudeNotLearnedTpmtime{i} = cell2mat(whiskingAmplitudeCell);
    midpointNotLearnedTpmtime{i} = cell2mat(whiskingMidpointCell);
end

for i = 1 : length(expertMi)
    mouse = mice(expertMi(i));
    session = sessions{expertMi(i)}(1);
    cd(sprintf('%s%03d',baseDir, mouse))
      
    ufn = sprintf('UberJK%03dS%02d',mouse,session);
    load(ufn)
    theta = cell2mat(cellfun(@(x) x.theta, u.trials, 'uniformoutput', false));
    [~, ~, thetaS, ~, amplitudeS, midpointS] = jkWhiskerDecomposition(theta);
    
    thetaBeforeLearning{i} = theta;
    amplitudeBeforeLearning{i} = amplitudeS;
    midpointBeforeLearning{i} = midpointS;
    
    whiskingOnsetCell = cell(1, length(u.trials));
    whiskingAmplitudeCell = cell(1, length(u.trials));
    whiskingMidpointCell = cell(1, length(u.trials));
    for ti = 1 : length(u.trials)
        currTrial = u.trials{ti};
        time = [0, currTrial.tpmTime{1}];
        wtimes = currTrial.whiskerTime;
        if iscell(wtimes)
            wtimes = wtimes{1};
        end
        [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
        whiskerVideoFrameDuration = u.trials{ti}.frameDuration; % in s
        onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
        whiskingOnsetCell{ti} = histcounts(onsetTimes, time);

        tempMid = zeros(1,length(time)-1);
        tempAmp = zeros(1,length(time)-1);
        for timei = 1 : length(tempMid)
            startInd = find(wtimes >= time(timei), 1, 'first');
            endInd = find(wtimes < time(timei+1), 1, 'last');
            tempMid(timei) = mean(midpoint(startInd:endInd));
            tempAmp(timei) = max(amplitude(startInd:endInd));
        end
        tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
        tempAmp(isnan(tempAmp)) = deal(mode(tempAmp(isfinite(tempAmp))));
        whiskingMidpointCell{ti} = tempMid;
        whiskingAmplitudeCell{ti} = tempAmp;
    end
    onsetBeforeLearningTpmtime{i} = cell2mat(whiskingOnsetCell);
    amplitudeBeforeLearningTpmtime{i} = cell2mat(whiskingAmplitudeCell);
    midpointBeforeLearningTpmtime{i} = cell2mat(whiskingMidpointCell);
end

for i = 1 : length(expertMi)
    mouse = mice(expertMi(i));
    session = sessions{expertMi(i)}(2);
    cd(sprintf('%s%03d',baseDir, mouse))
      
    ufn = sprintf('UberJK%03dS%02d',mouse,session);
    load(ufn)
    theta = cell2mat(cellfun(@(x) x.theta, u.trials, 'uniformoutput', false));
    [~, ~, thetaS, ~, amplitudeS, midpointS] = jkWhiskerDecomposition(theta);
    
    thetaAfterLearning{i} = theta;
    amplitudeAfterLearning{i} = amplitudeS;
    midpointAfterLearning{i} = midpointS;
    
    whiskingOnsetCell = cell(1, length(u.trials));
    whiskingAmplitudeCell = cell(1, length(u.trials));
    whiskingMidpointCell = cell(1, length(u.trials));
    for ti = 1 : length(u.trials)
        currTrial = u.trials{ti};
        time = [0, currTrial.tpmTime{1}];
        wtimes = currTrial.whiskerTime;
        if iscell(wtimes)
            wtimes = wtimes{1};
        end
        [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
        whiskerVideoFrameDuration = u.trials{ti}.frameDuration; % in s
        onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
        whiskingOnsetCell{ti} = histcounts(onsetTimes, time);

        tempMid = zeros(1,length(time)-1);
        tempAmp = zeros(1,length(time)-1);
        for timei = 1 : length(tempMid)
            startInd = find(wtimes >= time(timei), 1, 'first');
            endInd = find(wtimes < time(timei+1), 1, 'last');
            tempMid(timei) = mean(midpoint(startInd:endInd));
            tempAmp(timei) = max(amplitude(startInd:endInd));
        end
        tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
        tempAmp(isnan(tempAmp)) = deal(mode(tempAmp(isfinite(tempAmp))));
        whiskingMidpointCell{ti} = tempMid;
        whiskingAmplitudeCell{ti} = tempAmp;
    end
    onsetAfterLearningTpmtime{i} = cell2mat(whiskingOnsetCell);
    amplitudeAfterLearningTpmtime{i} = cell2mat(whiskingAmplitudeCell);
    midpointAfterLearningTpmtime{i} = cell2mat(whiskingMidpointCell);
end

%%
% gathering information from individual mice
thetaRange = -40:20;
amplitudeRange = 0:20;
midpointRange = -40:20;
onsetRange = 0:6;

histTheta = cell(6,3); % rows: mice, columns: 1 - nonlearner, 2 - before learning, 3 - after learning
histAmplitude = cell(6,3);
histMidpoint = cell(6,3);
stdTheta = zeros(6,3);
stdAmplitude = zeros(6,3);
stdMidpoint = zeros(6,3);

histOnset = cell(6,3);
histAmpTpm = cell(6,3);
histMidTpm = cell(6,3);
stdOnset = zeros(6,3);
stdAmpTpm = zeros(6,3);
stdMidTpm = zeros(6,3);

for mi = 1 : 6
    histTheta{mi,1} = histcounts(thetaNotLearned{mi}, thetaRange, 'normalization', 'probability');
    histAmplitude{mi,1} = histcounts(amplitudeNotLearned{mi}, amplitudeRange, 'normalization', 'probability');
    histMidpoint{mi,1} = histcounts(midpointNotLearned{mi}, midpointRange, 'normalization', 'probability');
    stdTheta(mi,1) = nanstd(thetaNotLearned{mi});
    stdAmplitude(mi,1) = nanstd(amplitudeNotLearned{mi});
    stdMidpoint(mi,1) = nanstd(midpointNotLearned{mi});
    
    histTheta{mi,2} = histcounts(thetaBeforeLearning{mi}, thetaRange, 'normalization', 'probability');
    histAmplitude{mi,2} = histcounts(amplitudeBeforeLearning{mi}, amplitudeRange, 'normalization', 'probability');
    histMidpoint{mi,2} = histcounts(midpointBeforeLearning{mi}, midpointRange, 'normalization', 'probability');
    stdTheta(mi,2) = nanstd(thetaBeforeLearning{mi});
    stdAmplitude(mi,2) = nanstd(amplitudeBeforeLearning{mi});
    stdMidpoint(mi,2) = nanstd(midpointBeforeLearning{mi});
    
    histTheta{mi,3} = histcounts(thetaAfterLearning{mi}, thetaRange, 'normalization', 'probability');
    histAmplitude{mi,3} = histcounts(amplitudeAfterLearning{mi}, amplitudeRange, 'normalization', 'probability');
    histMidpoint{mi,3} = histcounts(midpointAfterLearning{mi}, midpointRange, 'normalization', 'probability');
    stdTheta(mi,3) = nanstd(thetaAfterLearning{mi});
    stdAmplitude(mi, 3) = nanstd(amplitudeAfterLearning{mi});
    stdMidpoint(mi,3) = nanstd(midpointAfterLearning{mi});
    
    
    
    
    histOnset{mi,1} = histcounts(onsetNotLearnedTpmtime{mi}, onsetRange, 'normalization', 'probability');
    histAmpTpm{mi,1} = histcounts(amplitudeNotLearnedTpmtime{mi}, amplitudeRange, 'normalization', 'probability');
    histMidTpm{mi,1} = histcounts(midpointNotLearnedTpmtime{mi}, midpointRange, 'normalization', 'probability');
    stdOnset(mi,1) = std(onsetNotLearnedTpmtime{mi});
    stdAmpTpm(mi,1) = std(amplitudeNotLearnedTpmtime{mi});
    stdMidTpm(mi,1) = std(midpointNotLearnedTpmtime{mi});
    
    histOnset{mi,2} = histcounts(onsetBeforeLearningTpmtime{mi}, onsetRange, 'normalization', 'probability');
    histAmpTpm{mi,2} = histcounts(amplitudeBeforeLearningTpmtime{mi}, amplitudeRange, 'normalization', 'probability');
    histMidTpm{mi,2} = histcounts(midpointBeforeLearningTpmtime{mi}, midpointRange, 'normalization', 'probability');
    stdOnset(mi,2) = nanstd(onsetBeforeLearningTpmtime{mi});
    stdAmpTpm(mi,2) = nanstd(amplitudeBeforeLearningTpmtime{mi});
    stdMidTpm(mi,2) = nanstd(midpointBeforeLearningTpmtime{mi});
    
    histOnset{mi,3} = histcounts(onsetAfterLearningTpmtime{mi}, onsetRange, 'normalization', 'probability');
    histAmpTpm{mi,3} = histcounts(amplitudeAfterLearningTpmtime{mi}, amplitudeRange, 'normalization', 'probability');
    histMidTpm{mi,3} = histcounts(midpointAfterLearningTpmtime{mi}, midpointRange, 'normalization', 'probability');
    stdOnset(mi,3) = nanstd(onsetAfterLearningTpmtime{mi});
    stdAmpTpm(mi,3) = nanstd(amplitudeAfterLearningTpmtime{mi});
    stdMidTpm(mi,3) = nanstd(midpointAfterLearningTpmtime{mi});
end

figure,
subplot(2,4,1), hold on
boundedline(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,1))), std(cell2mat(histTheta(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,2))), std(cell2mat(histTheta(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,3))), std(cell2mat(histTheta(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,1))), 'c-', 'linewidth', 3)
plot(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,2))), 'b-', 'linewidth', 3)
plot(thetaRange(1:end-1)+0.5, mean(cell2mat(histTheta(:,3))), 'r-', 'linewidth', 3)
xlabel('Angle (\circ)')
ylabel('Proportion')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,2), hold on
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,1))), std(cell2mat(histAmplitude(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,2))), std(cell2mat(histAmplitude(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,3))), std(cell2mat(histAmplitude(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,1))), 'c-', 'linewidth', 3)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,2))), 'b-', 'linewidth', 3)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmplitude(:,3))), 'r-', 'linewidth', 3)
xlabel('Amplitude (\circ)')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,3), hold on
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,1))), std(cell2mat(histMidpoint(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,2))), std(cell2mat(histMidpoint(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,3))), std(cell2mat(histMidpoint(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,1))), 'c-', 'linewidth', 3)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,2))), 'b-', 'linewidth', 3)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidpoint(:,3))), 'r-', 'linewidth', 3)
xlabel('Midpoint (\circ)')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,4), hold on
bar(0.7, mean(stdTheta(:,1)), 0.3, 'c')
bar(1.0, mean(stdTheta(:,2)), 0.3, 'b')
bar(1.3, mean(stdTheta(:,3)), 0.3, 'r')
errorbar(0.7, mean(stdTheta(:,1)), std(stdTheta(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(1.0, mean(stdTheta(:,2)), std(stdTheta(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(1.3, mean(stdTheta(:,3)), std(stdTheta(:,3))/sqrt(6), 'r', 'linewidth', 3)
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

bar(1.7, mean(stdAmplitude(:,1)), 0.3, 'c')
bar(2.0, mean(stdAmplitude(:,2)), 0.3, 'b')
bar(2.3, mean(stdAmplitude(:,3)), 0.3, 'r')
errorbar(1.7, mean(stdAmplitude(:,1)), std(stdAmplitude(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(2.0, mean(stdAmplitude(:,2)), std(stdAmplitude(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(2.3, mean(stdAmplitude(:,3)), std(stdAmplitude(:,3))/sqrt(6), 'r', 'linewidth', 3)
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

bar(2.7, mean(stdMidpoint(:,1)), 0.3, 'c')
bar(3.0, mean(stdMidpoint(:,2)), 0.3, 'b')
bar(3.3, mean(stdMidpoint(:,3)), 0.3, 'r')
errorbar(2.7, mean(stdMidpoint(:,1)), std(stdMidpoint(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(3.0, mean(stdMidpoint(:,2)), std(stdMidpoint(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(3.3, mean(stdMidpoint(:,3)), std(stdMidpoint(:,3))/sqrt(6), 'r', 'linewidth', 3)
xticks([1:3]), xticklabels({'Theta', 'Amplitude', 'Midpoint'}), xtickangle(45)
ylabel('Standard deviation')
legend({'Nonlearner', 'Before learning', 'After learning'}, 'box', 'off', 'location', 'northeast')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')




subplot(2,4,5), hold on
boundedline(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,1))), std(cell2mat(histOnset(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,2))), std(cell2mat(histOnset(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,3))), std(cell2mat(histOnset(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,1))), 'c-', 'linewidth', 3)
plot(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,2))), 'b-', 'linewidth', 3)
plot(onsetRange(1:end-1)+0.5, mean(cell2mat(histOnset(:,3))), 'r-', 'linewidth', 3)
xlabel('# of whisking onsets')
ylabel('Proportion')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,6), hold on
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,1))), std(cell2mat(histAmpTpm(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,2))), std(cell2mat(histAmpTpm(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,3))), std(cell2mat(histAmpTpm(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,1))), 'c-', 'linewidth', 3)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,2))), 'b-', 'linewidth', 3)
plot(amplitudeRange(1:end-1)+0.5, mean(cell2mat(histAmpTpm(:,3))), 'r-', 'linewidth', 3)
xlabel('Amplitude (\circ)')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,7), hold on
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,1))), std(cell2mat(histMidTpm(:,1))) / sqrt(6), 'c-', 'transparency', 0.2)        
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,2))), std(cell2mat(histMidTpm(:,2))) / sqrt(6), 'b-', 'transparency', 0.2)        
boundedline(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,3))), std(cell2mat(histMidTpm(:,3))) / sqrt(6), 'r-', 'transparency', 0.2)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,1))), 'c-', 'linewidth', 3)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,2))), 'b-', 'linewidth', 3)
plot(midpointRange(1:end-1)+0.5, mean(cell2mat(histMidTpm(:,3))), 'r-', 'linewidth', 3)
xlabel('Midpoint (\circ)')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

subplot(2,4,8), hold on
bar(0.7, mean(stdOnset(:,1)), 0.3, 'c')
bar(1.0, mean(stdOnset(:,2)), 0.3, 'b')
bar(1.3, mean(stdOnset(:,3)), 0.3, 'r')
errorbar(0.7, mean(stdOnset(:,1)), std(stdOnset(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(1.0, mean(stdOnset(:,2)), std(stdOnset(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(1.3, mean(stdOnset(:,3)), std(stdOnset(:,3))/sqrt(6), 'r', 'linewidth', 3)
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

bar(1.7, mean(stdAmpTpm(:,1)), 0.3, 'c')
bar(2.0, mean(stdAmpTpm(:,2)), 0.3, 'b')
bar(2.3, mean(stdAmpTpm(:,3)), 0.3, 'r')
errorbar(1.7, mean(stdAmpTpm(:,1)), std(stdAmpTpm(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(2.0, mean(stdAmpTpm(:,2)), std(stdAmpTpm(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(2.3, mean(stdAmpTpm(:,3)), std(stdAmpTpm(:,3))/sqrt(6), 'r', 'linewidth', 3)
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')

bar(2.7, mean(stdMidTpm(:,1)), 0.3, 'c')
bar(3.0, mean(stdMidTpm(:,2)), 0.3, 'b')
bar(3.3, mean(stdMidTpm(:,3)), 0.3, 'r')
errorbar(2.7, mean(stdMidTpm(:,1)), std(stdMidTpm(:,1))/sqrt(6), 'c', 'linewidth', 3)
errorbar(3.0, mean(stdMidTpm(:,2)), std(stdMidTpm(:,2))/sqrt(6), 'b', 'linewidth', 3)
errorbar(3.3, mean(stdMidTpm(:,3)), std(stdMidTpm(:,3))/sqrt(6), 'r', 'linewidth', 3)
xticks([1:3]), xticklabels({'Onset', 'Amplitude', 'Midpoint'}), xtickangle(45)
ylabel('Standard deviation')
set(gca, 'linewidth', 3, 'fontsize', 12, 'fontweight', 'bold')