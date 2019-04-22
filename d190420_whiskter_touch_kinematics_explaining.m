%% The influence of whisker kinematics
%% Compare between drop-out fitting and drop-out prediction.
%% DE diff comparison, DE comparison, ratio of error ratio comparison, etc.
%% Try finding out anything that makes both of the method be correlated with each other.
%% If not, try picking out the ones that are correlated, saying the others are hard to assign.
%% The question is, how many of them can be effectively explained by whisker touch kinematics?
%% Which touch kinematics best describe those cells?

%% first, look at difference in DE diff in each one
clear
baseDir = 'C:\JK\';
cd(baseDir)
fullModel = load('glm_cell_function_error_ratio_withWTV_shuffling', 'naive', 'expert');
wtvModel = load('glm_cell_function_error_ratio_WTV_ONLY', 'naive', 'expert');
touchModel = load('glm_results_responseType', 'naive', 'expert');

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

%% Are the two methods correlated?

figure, 
subplot(121), hold on
for i = 1 : 12

    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) - DEdiffFromFull(:,2), DEdiffFromPartial(:,1) - DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch - wtv)')
ylabel('DE diff from partial models (touch - wtv)')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) - DEdiffFromFull(:,2), DEdiffFromPartial(:,1) - DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch - wtv)')
ylabel('DE diff from partial models (touch - wtv)')
title('Expert (n=6)')

%% Results: both naive and expert have correlated relationship between two

%% ratio of DE diff

figure, 
subplot(121), hold on
% for i = 1 : 12
for i = 1 : 8
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) ./ DEdiffFromFull(:,2), DEdiffFromPartial(:,1) ./ DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch / wtv)')
ylabel('DE diff from partial models (touch / wtv)')
title('All naive (n=8*)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromFull(:,1) ./ DEdiffFromFull(:,2), DEdiffFromPartial(:,1) ./ DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff from full model (touch / wtv)')
ylabel('DE diff from partial models (touch / wtv)')
title('Expert (n=6)')


%% Results: There are some negative values. Mostly (large ones) from partial wvt model (touch ONLY model)
% Which means that in these cases adding touch was making the model worse
% There are some (very few) very high values in partial model, which means partial wvt model did not explain any.
% But majority of the data are hard to explain.

%% Compare between touch and wtv in each method
% is there negative correlation? how are the values distributed?
% (1) partial fitting method

figure, 
subplot(121), hold on
for i = 1 : 12
% for i = 1 : 8
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    cinds = find(ismember(touchModel.expert(i).cellID, cID));
    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?

    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')


%% (2) partial prediction method

figure, 
subplot(121), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('All naive (n=12)')

subplot(122), hold on
for i = 1 : 6

    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 

    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert (n=6)')

%% Results:
% They are not negatively correlated in both ways. 
% partial fitting method is more difficult to interpret because of negative values.


%% What are the values in partial prediction method for those of negative values is partial fittin method?
touchNegvalInd = cell(12,1);
wtvNegvalInd = cell(12,1);
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    touchNegvalInd{i} = find(DEdiffFromPartial(:,1) < 0);
    wtvNegvalInd{i} = find(DEdiffFromPartial(:,2) < 0);
end

figure,
subplot(221), hold on
for i = 1 : 12
    
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(touchNegvalInd{i},1), DEdiffFromFull(touchNegvalInd{i},2), 'r.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Naive touch only better than full model')


subplot(222), hold on
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(wtvNegvalInd{i},1), DEdiffFromFull(wtvNegvalInd{i},2), 'b.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Naive WTV only better than full model')
    

% now for expert
touchNegvalInd = cell(6,1);
wtvNegvalInd = cell(6,1);
for i = 1 : 6
    cID = wtvModel.expert(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.expert(i).cellID, cID));

    DEdiffFromPartial(:,1) = fullModel.expert(i).devExp - wtvModel.expert(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.expert(i).devExp - touchModel.expert(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    touchNegvalInd{i} = find(DEdiffFromPartial(:,1) < 0);
    wtvNegvalInd{i} = find(DEdiffFromPartial(:,2) < 0);
end


subplot(223), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(touchNegvalInd{i},1), DEdiffFromFull(touchNegvalInd{i},2), 'r.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert touch only better than full model')


subplot(224), hold on
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(:,1), DEdiffFromFull(:,2), 'k.')
end
for i = 1 : 6
    fullDEdiff = cell2mat(fullModel.expert(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);

    plot(DEdiffFromFull(wtvNegvalInd{i},1), DEdiffFromFull(wtvNegvalInd{i},2), 'b.')
end
xlabel('DE diff touch')
ylabel('DE diff WTV')
title('Expert WTV only better than full model')

%% Results: scattered a lot. Ther eis a tendency, but very weak.

%% How does DE diff < 0.1 from partial prediction look like in partial fitting?
% naive
touchLowvalInd = cell(12,1);
wtvLowvalInd = cell(12,1);
for i = 1 : 12
    fullDEdiff = cell2mat(fullModel.naive(i).DEdiff);
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    touchLowvalInd{i} = find(DEdiffFromFull(:,1) < 0.1);
    wtvLowvalInd{i} = find(DEdiffFromFull(:,2) < 0.1);
end

figure,
subplot(221), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(touchLowvalInd{i},1), DEdiffFromPartial(touchLowvalInd{i},2), 'r.')
end

xlabel('DE diff touch from partial fitting')
ylabel('DE diff WTV from partial fitting')
title('Naive without touch similar to full model')

subplot(222), hold on
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(:,1), DEdiffFromPartial(:,2), 'k.')
end
for i = 1 : 12
    cID = wtvModel.naive(i).cID; % touch cells
    % DEdiffFromFull = zeros(length(cID),2); % 1 for touch, 2 for wtv
    DEdiffFromPartial = zeros(length(cID),2); 
    cinds = find(ismember(touchModel.naive(i).cellID, cID));
    
    DEdiffFromPartial(:,1) = fullModel.naive(i).devExp - wtvModel.naive(i).devExp; % how much does touch predictors affect the prediction?
    DEdiffFromPartial(:,2) = fullModel.naive(i).devExp - touchModel.naive(i).allDE(cinds); % how much does wtv predictors affect the prediction?
    plot(DEdiffFromPartial(wtvLowvalInd{i},1), DEdiffFromPartial(wtvLowvalInd{i},2), 'b.')
end

xlabel('DE diff touch from partial fitting')
ylabel('DE diff WTV from partial fitting')
title('Naive without WTV similar to full model')

%% Results: again, there is a tendency, but very weak and scattered or covers most of the points.
%% Negative values from partial fitting are not necessarily low values in partial prediction (although there is a tendency).
%% It makes it difficult to connect both of them, so just focus on partial prediction method. I have shuffling for it too.
%% It is good enough to know that there is a correlation.

%% Try finding the right threshold from shuffling
% std of each DE diff, calculated back from error ratio (stupid...)
% plot each partial DE diff value to each std (one plot for touch another
% for WTV)
clear
baseDir = 'C:\JK\';
cd(baseDir)
load('glm_cell_function_error_ratio_withWTV_shuffling', 'naive', 'expert');

touchSTD = cell(12,1);
WTVSTD = cell(12,1);
for i = 1 : 12
    fullDEdiff = cell2mat(naive(i).DEdiff);
    
    DEdiffFromFull = fullDEdiff(:,[1,6]);
    
    touchSTD{i} = find(DEdiffFromFull(:,1) < 0.1);
    WTVSTD{i} = find(DEdiffFromFull(:,2) < 0.1);
end




