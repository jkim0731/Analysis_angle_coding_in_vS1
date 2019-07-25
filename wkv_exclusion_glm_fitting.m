% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'D:\TPM\JK\suite2p\';
loadfn = 'glmResults_WKV_touchCell_exclusion';
load(sprintf('%s%s',baseDir,loadfn))

% whiskerTouchMat = [maxDkappaHMat, maxDkappaVMat, maxDthetaMat, maxDphiMat, maxSlideDistanceMat, maxDurationMat, ...    
%                             thetaAtTouchMat, phiAtTouchMat, kappaHAtTouchMat, kappaVAtTouchMat, arcLengthAtTouchMat, touchCountMat];
presentationOrder = [3,4,1,2,5,6, ...
                        7,8,9,10,11,12];
barWidth = 0.2;
                    
                    
expertInd = [1:4,7,9];
nonlearnerInd = setdiff(1:12,expertInd);

numCellNonlearner = zeros(length(nonlearnerInd),1);
for i = 1 : length(nonlearnerInd)
    numCellNonlearner(i) = length(naive(nonlearnerInd(i)).cID);
end
numCellNaive = zeros(length(expertInd),1);
for i = 1 : length(expertInd)
    numCellNaive(i) = length(naive(expertInd(i)).cID);
end
numCellExpert = zeros(length(expertInd),1);
for i = 1 : length(expertInd)
    numCellExpert(i) = length(expert(i).cID);
end

cumsumNonlearner = [0;cumsum(numCellNonlearner)];
cumsumNaive = [0;cumsum(numCellNaive)];
cumsumExpert = [0;cumsum(numCellExpert)];

deDiffNonlearner = zeros(cumsumNonlearner(end),12);
deDiffNaive = zeros(cumsumNaive(end),12);
deDiffExpert = zeros(cumsumExpert(end),12);

for i = 1 : length(nonlearnerInd)    
    deDiffNonlearner(cumsumNonlearner(i)+1:cumsumNonlearner(i+1),:) = naive(nonlearnerInd(i)).whiskerVariableDEdiff(:,presentationOrder);
end
for i = 1 : length(expertInd)    
    deDiffNaive(cumsumNaive(i)+1:cumsumNaive(i+1),:) = naive(expertInd(i)).whiskerVariableDEdiff(:,presentationOrder);
    deDiffExpert(cumsumExpert(i)+1:cumsumExpert(i+1),:) = expert(i).whiskerVariableDEdiff(:,presentationOrder);
end

%%
figure('units', 'normalized', 'outerposition', [0.1 0.1 0.6 0.5]), hold on
bar([[1:6, 8:13]]-barWidth, mean(deDiffNonlearner), barWidth, 'c', 'edgecolor', 'c')
bar([[1:6, 8:13]], mean(deDiffNaive), barWidth, 'b', 'edgecolor', 'b')
bar([[1:6, 8:13]]+barWidth, mean(deDiffExpert), barWidth, 'r', 'edgecolor', 'r')
errorbar([[1:6, 8:13]]-barWidth, mean(deDiffNonlearner), std(deDiffNonlearner)/sqrt(size(deDiffNonlearner,1)),'c.')
errorbar([[1:6, 8:13]], mean(deDiffNaive), std(deDiffNaive)/sqrt(size(deDiffNaive,1)), 'b.')
errorbar([[1:6, 8:13]]+barWidth, mean(deDiffExpert), std(deDiffExpert)/sqrt(size(deDiffExpert,1)), 'r.')
xticks([1:6, 8:13])
xticklabels({'max\Delta\theta', 'max\Delta\phi', 'max\Delta\kappa_H', 'max\Delta\kappa_V', 'Slide distance', 'Duration', ...
    '\theta', '\phi', '\kappa_H', '\kappa_V', 'Arc length', 'Touch count'})
xtickangle(45)
set(gca, 'fontsize', 16)
ylabel('Deviance explained difference')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'northeast')