%% Check for similarity between touch and choice
% because they look suspiciously similar (not exactly the same, though)
dir = 'C:\Users\shires\Documents\GitHub\AngleDiscrimBehavior\matlab\datastructs\';
Xhow = 'Mean'; %'Individual' or 'Mean'
learned = 'Expert'; % 'Naive', 'Expert', or ''
timing = 'lick'; % 'lick' or 'answer'
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'

Yout = 'Touch'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
touchData = load([dir, fn]);

Yout = 'Choice'; % 'Touch' or 'Choice'
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
choiceData = load([dir, fn]);

numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);
histRange = -3:0.2:3.1;
%% Check if the data structure is correct

for mi = 1:numMice
    compare(touchData.groupMdl{mi}.io.X, choiceData.groupMdl{mi}.io.X)
end
performance = zeros(numMice,1);
for mi = 1 : numMice
    tempY = choiceData.groupMdl{mi}.io.Y;
    angleY = touchData.groupMdl{mi}.io.Y;
    numDiff = length(find(abs(tempY-angleY)));
    performance(mi) = 1-mean(abs(tempY-angleY));
end
performance

%% results:
% data structure is simple. The only difference is in Y, and the performance matches





%%
choice = sort(unique(choiceData.groupMdl{1}.io.Y), 'descend');
choiceFeatures = cell(length(choiceData.groupMdl), length(choice), size(choiceData.groupMdl{1}.io.X,2));
for mi = 1 : length(choiceData.groupMdl)
    for ci = 1 : length(choice)
        tempInd = find(choiceData.groupMdl{mi}.io.Y == choice(ci));
        for i = 1 : size(choiceFeatures,3)
            choiceFeatures{mi,ci,i} = choiceData.groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

angles = sort(unique(touchData.groupMdl{1}.io.Y), 'descend');
touchFeatures = cell(length(touchData.groupMdl), length(angles), size(touchData.groupMdl{1}.io.X,2));
for mi = 1 : length(touchData.groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(touchData.groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(touchFeatures,3)
            touchFeatures{mi,ai,i} = touchData.groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

%%

fi = 2;
colorList = {'b','r';'c','m'};
figure
for mi = 1 : 6
    subplot(2,3,mi), hold on
    for aiorci = 1 : 2        
        choiceTrace = choiceFeatures{mi,aiorci,fi};
        touchTrace = touchFeatures{mi,aiorci,fi};
        choiceHist = histcounts(choiceTrace, histRange, 'normalization', 'probability');
        touchHist = histcounts(touchTrace, histRange, 'normalization', 'probability');
        plot(histRange(1:end-1), choiceHist, 'color', colorList{2,aiorci}, 'linewidth', 2);
        plot(histRange(1:end-1), touchHist, 'color', colorList{1,aiorci}, 'linewidth', 2);        
    end
    if mi == 2
        title({touchData.groupMdl{1}.fitCoeffsFields{fi};sprintf('Performance %.3f', performance(mi))})
    else
        title(sprintf('Performance %.3f', performance(mi)))
    end
    ylim([0 0.4])
end



%%

fi = 2;
colorList = {'b','r';'c','m'};
figure
for mi = 1 : 6
    subplot(2,3,mi), hold on
    for aiorci = 1 : 2      
        totalChoiceLength = sum(cellfun(@length, choiceFeatures(mi,:,fi)));
        totalTouchLength = sum(cellfun(@length, touchFeatures(mi,:,fi)));
        choiceTrace = choiceFeatures{mi,aiorci,fi};
        touchTrace = touchFeatures{mi,aiorci,fi};
        choiceHist = histcounts(choiceTrace, histRange, 'normalization', 'probability');
        touchHist = histcounts(touchTrace, histRange, 'normalization', 'probability');
        plot(histRange(1:end-1), choiceHist*length(choiceTrace)/totalChoiceLength, 'color', colorList{2,aiorci}, 'linewidth', 2);
        plot(histRange(1:end-1), touchHist*length(touchTrace)/totalTouchLength, 'color', colorList{1,aiorci}, 'linewidth', 2);        
    end
    if mi == 2
        title({touchData.groupMdl{1}.fitCoeffsFields{fi};sprintf('Performance %.3f', performance(mi))})
    else
        title(sprintf('Performance %.3f', performance(mi)))
    end
    ylim([0 0.2])
end









%% Simulation of the difference in angle and choice
% 1. Match the number of errors, and randomly assign their choice (multiple times), and see how the distribution difference changes
% 2. Gradually increase the mismatch in choice, by 5 % each, and see how the distribution difference changes

%% Method 1
histRange = -3:0.2:3.1;
simNum = 1000;
numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);
choiceFeaturesSim = cell(numMice, length(choice), numFeat, simNum);
performance = zeros(numMice,1);
for mi = 1 : numMice
    tempY = choiceData.groupMdl{mi}.io.Y;
    angleY = touchData.groupMdl{mi}.io.Y;
    numDiff = length(find(abs(tempY-angleY)));
    performance(mi) = 1-mean(abs(tempY-angleY));
    
    choice = sort(unique(tempY), 'descend');
    
    
    for si = 1 : simNum   
        tempY = angleY;
        ind2flip = randperm(length(angleY), numDiff);
        tempY(ind2flip) = 1-tempY(ind2flip);

        for ci = 1 : length(choice)
            tempInd = find(tempY == choice(ci));
            for fi = 1 : numFeat
                choiceFeaturesSim{mi,ci,fi,si} = choiceData.groupMdl{mi}.io.X(tempInd,fi);
            end
        end
    end
end
%%
figure,
fi = 2;
for mi = 1:numMice
    subplot(2,3,mi), hold on    
    for aiorci = 1 : 2
%         for si = 1 : simNum
%             choiceTraceSim = choiceFeaturesSim{mi,aiorci,fi,si};
%             choiceHistSim = histcounts(choiceTraceSim, histRange, 'normalization', 'probability');
%             plot(histRange(1:end-1), choiceHistSim, 'color', colorList{2,aiorci});
%         end
        choiceHistSim = zeros(simNum,length(histRange)-1);
        coeffProportion = length(choiceFeatures{mi,aiorci,fi}) / (length(choiceFeatures{mi,1,fi}) + length(choiceFeatures{mi,2,fi}));
        for si = 1 : simNum
            choiceTraceSim = choiceFeaturesSim{mi,aiorci,fi,si};
            choiceHistSim(si,:) = histcounts(choiceTraceSim, histRange, 'normalization', 'probability') * coeffProportion;
        end
        yMean = mean(choiceHistSim);
        ySEM = std(choiceHistSim)/sqrt(simNum);
        CI95 = tinv([0.005 0.995], simNum-1);
        yCI95 = bsxfun(@times, ySEM, CI95(:));
        plot(histRange(1:end-1), yCI95+yMean, 'color', colorList{2,aiorci})
        
        choiceTrace = choiceFeatures{mi,aiorci,fi};
        choiceHist = histcounts(choiceTrace, histRange, 'normalization', 'probability') * coeffProportion;
        plot(histRange(1:end-1), choiceHist, 'color', colorList{1,aiorci}, 'linewidth', 2);
        if mi == 2
            title({touchData.groupMdl{1}.fitCoeffsFields{fi};sprintf('Performance %.3f', performance(mi))})
        else
            title(sprintf('Performance %.3f', performance(mi)))
        end
        ylim([0 0.25])
%         touchTrace = touchFeatures{mi,aiorci,fi};
%         touchHist = histcounts(touchTrace, histRange, 'normalization', 'probability');
%         plot(histRange(1:end-1), touchHist, 'color', colorList{1,aiorci}, 'linewidth', 2);
        
    end
    
end

% legend({'Right 99% CI', 'Right', 'Left 99% CI', 'Left'})


%% Method 2
histRange = -3:0.2:3.1;
errorGroupNum = 11; % 100:-5:50 %
simNum = 1000;
numMice = length(choiceData.groupMdl);
numFeat = size(choiceData.groupMdl{1}.io.X,2);
choiceFeaturesSim = cell(numMice, length(choice), numFeat, errorGroupNum,simNum);
performance = zeros(numMice,1);
for mi = 1 : numMice
    tempY = choiceData.groupMdl{mi}.io.Y;
    angleY = touchData.groupMdl{mi}.io.Y;
    performance(mi) = 1-mean(abs(tempY-angleY));
    
    for ei = 1 : errorGroupNum
        simPerf = 1-0.05*(ei-1);
        numDiff = round(length(angleY)*(1-simPerf));
        for si = 1 : simNum   
            tempY = angleY;
            ind2flip = randperm(length(angleY), numDiff);
            tempY(ind2flip) = 1-tempY(ind2flip);

            for ci = 1 : length(choice)
                tempInd = find(tempY == choice(ci));
                for fi = 1 : numFeat
                    choiceFeaturesSim{mi,ci,fi,ei,si} = choiceData.groupMdl{mi}.io.X(tempInd,fi);
                end
            end
        end
    end
end

%%
figure,
fi = 4;
for mi = 1:numMice
    subplot(2,3,mi), hold on
    for aiorci = 1 : 2
        coeffProportion = length(choiceFeatures{mi,aiorci,fi}) / (length(choiceFeatures{mi,1,fi}) + length(choiceFeatures{mi,2,fi}));
        for ei = 1 : errorGroupNum
            choiceHistSim = zeros(simNum,length(histRange)-1);
            for si = 1 : simNum
                choiceTraceSim = choiceFeaturesSim{mi,aiorci,fi,ei,si};
                choiceHistSim(si,:) = histcounts(choiceTraceSim, histRange, 'normalization', 'probability');            
            end
            yMean = mean(choiceHistSim)*coeffProportion;
            if aiorci == 1
                plot(histRange(1:end-1), yMean, 'color', [(ei-1)*0.08 (ei-1)*0.08 1])
            else
                plot(histRange(1:end-1), yMean, 'color', [1 (ei-1)*0.08 (ei-1)*0.08])
            end
            
        end
        
        choiceTrace = choiceFeatures{mi,aiorci,fi};
        choiceHist = histcounts(choiceTrace, histRange, 'normalization', 'probability') * coeffProportion;
        plot(histRange(1:end-1), choiceHist, 'color', colorList{2,aiorci}, 'linewidth', 2);
        if mi == 2
            title({touchData.groupMdl{1}.fitCoeffsFields{fi};sprintf('Performance %.3f', performance(mi))})
        else
            title(sprintf('Performance %.3f', performance(mi)))
        end
    end
end