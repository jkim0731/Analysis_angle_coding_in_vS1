% baseDir = 'Y:\Whiskernas\JK\suite2p\';
baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
touch = load('glmResults_devExp_touch_NC');
wkv = load('glmResults_devExp_WKV_touchCell_NC');

expertInd = [1:4,7,9];
%% all naive
histRange = [-0.1:0.01:0.2];
diffDE = cell(length(wkv.naive),1);
for i = 1 : length(wkv.naive)
    tempCID = wkv.naive(i).cellID;
%     diffDE{i} = zeros(length(tempCID),1);
    pfTemp = zeros(length(tempCID),1); % parfor temp
    parfor ci = 1 : length(tempCID)
        wkvInd = find(wkv.naive(i).cellID == tempCID(ci));
        touchInd = find(touch.naive(i).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.naive(i).allDE(touchInd) - wkv.naive(i).allDE(wkvInd);
    end
    diffDE{i} = pfTemp;
end

% figure,
% histogram(cell2mat(diffDE), range, 'normalization', 'cdf')


%% expert
histRange = [-0.1:0.01:0.2];
diffDE = cell(length(wkv.expert),1);
for i = 1 : length(wkv.expert)
    tempCID = wkv.expert(i).cellID;
%     diffDE{i} = zeros(length(tempCID),1);
    pfTemp = zeros(length(tempCID),1); % parfor temp
    parfor ci = 1 : length(tempCID)
        wkvInd = find(wkv.expert(i).cellID == tempCID(ci));
        touchInd = find(touch.expert(i).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.expert(i).allDE(touchInd) - wkv.expert(i).allDE(wkvInd);
    end
    diffDE{i} = pfTemp;
end

figure,
histogram(cell2mat(diffDE), histRange, 'normalization', 'cdf')

%% cdf of nonlearner, naive, and expert

histRange = [-0.11:0.01:0.2];
nonlearnerInd = setdiff([1:1:12], expertInd);
nonlearnerCDF = zeros(length(nonlearnerInd),length(histRange)-1);
naiveCDF = zeros(length(expertInd),length(histRange)-1);
expertCDF = zeros(length(expertInd),length(histRange)-1);

for i = 1 : length(nonlearnerInd)
    tempCID = wkv.naive(nonlearnerInd(i)).cellID;
    pfTemp = zeros(length(tempCID),1); % parfor temp
    parfor ci = 1 : length(tempCID)
        wkvInd = find(wkv.naive(nonlearnerInd(i)).cellID == tempCID(ci));
        touchInd = find(touch.naive(nonlearnerInd(i)).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.naive(nonlearnerInd(i)).allDE(touchInd) - wkv.naive(nonlearnerInd(i)).allDE(wkvInd);
    end
    nonlearnerCDF(i,:) = histcounts(pfTemp, histRange, 'normalization', 'cdf');
end

for i = 1 : length(expertInd)
    tempCID = wkv.naive(expertInd(i)).cellID;
    pfTemp = zeros(length(tempCID),1); % parfor temp
    parfor ci = 1 : length(tempCID)
        wkvInd = find(wkv.naive(expertInd(i)).cellID == tempCID(ci));
        touchInd = find(touch.naive(expertInd(i)).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.naive(expertInd(i)).allDE(touchInd) - wkv.naive(expertInd(i)).allDE(wkvInd);
    end
    naiveCDF(i,:) = histcounts(pfTemp, histRange, 'normalization', 'cdf');
end
%%
for i = 1 : length(expertInd)
    tempCID = wkv.expert(i).cellID;
    pfTemp = zeros(length(tempCID),1); % parfor temp
    parfor ci = 1 : length(tempCID)
        wkvInd = find(wkv.expert(i).cellID == tempCID(ci));
        touchInd = find(touch.expert(i).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.expert(i).allDE(touchInd) - wkv.expert(i).allDE(wkvInd);
    end
    expertCDF(i,:) = histcounts(pfTemp, histRange, 'normalization', 'cdf');
end

figure, hold on,
shadedErrorBar(histRange(2:end), mean(nonlearnerCDF), std(nonlearnerCDF)/sqrt(length(nonlearnerInd)), 'lineprop', 'c-')
shadedErrorBar(histRange(2:end), mean(naiveCDF), std(naiveCDF)/sqrt(length(expertInd)), 'lineprop', 'b-')
shadedErrorBar(histRange(2:end), mean(expertCDF), std(expertCDF)/sqrt(length(expertInd)), 'lineprop', 'r-')
ylabel('Cumulative proportion')
xlabel('DE(touch) - DE(whisker)')
legend({'Nonlearner', 'Naive', 'Expert'}, 'box', 'off', 'location', 'southeast')
set(gca, 'fontsize', 16)