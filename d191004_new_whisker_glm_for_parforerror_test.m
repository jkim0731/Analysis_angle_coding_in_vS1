dirnew = 'Y:\Whiskernas\JK\suite2p\025\';
dirold = 'Y:\Whiskernas\JK\suite2p\025\whisker_allCell_OG\';

fnBase = 'glmWhisker_lasso_allCell_NC_JK025S04_R';

%% load one for setting

load(sprintf('%s%s%02d',dirold, fnBase, 1), 'done');
numCell = length(done);
%% load and average
DEoldAll = zeros(numCell, 10);
DEnewAll = zeros(numCell, 10);

coeffOldAll = cell(numCell, 10);
coeffNewAll = cell(numCell, 10);
for i = 1 : 10
    load(sprintf('%s%s%02d',dirold, fnBase, i), 'fitCoeffs', 'fitDevExplained')
    DEoldAll(:,i) = fitDevExplained;
    coeffOldAll(:,i) = fitCoeffs;
    
    load(sprintf('%s%s%02d',dirnew, fnBase, i), 'fitCoeffs', 'fitDevExplained')
    DEnewAll(:,i) = fitDevExplained;
    coeffNewAll(:,i) = fitCoeffs;    
end

DEold = mean(DEoldAll,2);
DEnew = mean(DEnewAll,2);
coeffOld = cell(numCell,1);
coeffNew = cell(numCell,1);
%%
for i = 1 : numCell
    coeffOld{i} = mean(cell2mat(coeffOldAll(i,:)),2);
    coeffNew{i} = mean(cell2mat(coeffNewAll(i,:)),2);
end

%% Compare average DE (scatter)

figure
scatter(DEold, DEnew)


%% Compare max Coeff > 0.05 identity
maxNew = zeros(numCell,1);
maxOld = zeros(numCell,1);

for i = 1 : numCell
    if ~isempty(coeffNew{i})
        [~, maxNew(i)] = max(coeffNew{i}(2:end));
    end    
    if ~isempty(coeffOld{i})
        [~, maxOld(i)] = max(coeffOld{i}(2:end));
    end
end
oldFit = find(DEold >= 0.1);
newFit = find(DEnew >= 0.1);
bothFit = intersect(oldFit, newFit);

figure, 
subplot(131), scatter(maxOld(oldFit), maxNew(oldFit), '.')
subplot(132), scatter(maxOld(newFit), maxNew(newFit), '.')
subplot(133), scatter(maxOld(bothFit), maxNew(bothFit), '.')


