tuneNew = load('angle_tuning_summary_preAnswer_perTouch_NC');
tuneOld = load('angle_tuning_summary_predecision_NC');

%%
% compare(tuneOld.naive(1).tuned, tuneNew.naive(1).tuned)

tuned = zeros(length(tuneNew.naive),2);
diffs = zeros(length(tuneNew.naive),2); % 1: abs #, 2: proportion
for i = 1 : length(tuneNew.naive)
    tuned(i,1) = sum(tuneOld.naive(i).tuned);
    tuned(i,2) = sum(tuneNew.naive(i).tuned);
    diffs(i,1) = length(find(abs(tuneOld.naive(i).tuned - tuneNew.naive(i).tuned)==1));
    diffs(i,2) = length(find(abs(tuneOld.naive(i).tuned - tuneNew.naive(i).tuned)==1))/length(tuneOld.naive(i).touchID);
end

mean(tuned(:,1) - tuned(:,2))
std(tuned(:,1) - tuned(:,2))

mean(diffs(:,2))
std(diffs(:,2))

%%

tuned = zeros(length(tuneNew.expert),2);
diffs = zeros(length(tuneNew.expert),2); % 1: abs #, 2: proportion
for i = 1 : length(tuneNew.expert)
    tuned(i,1) = sum(tuneOld.expert(i).tuned);
    tuned(i,2) = sum(tuneNew.expert(i).tuned);
    diffs(i,1) = length(find(abs(tuneOld.expert(i).tuned - tuneNew.expert(i).tuned)==1));
    diffs(i,2) = length(find(abs(tuneOld.expert(i).tuned - tuneNew.expert(i).tuned)==1))/length(tuneOld.expert(i).touchID);
end

mean(tuned(:,1) - tuned(:,2))
std(tuned(:,1) - tuned(:,2))

mean(diffs(:,2))
std(diffs(:,2))


%%
% tunedAngle
% unimodalSingle
% unimodalBroad
% multimodal

diffs = zeros(length(tuneNew.naive),4); % 1: tunedAngle, 2: unimodalSingle, 3: unimodalBroad, 4: multimodal
for i = 1 : length(tuneNew.naive)
    matchInds = intersect(find(tuneNew.naive(i).tuned), find(tuneOld.naive(i).tuned));
    diffs(i,1) = length(find(abs(tuneOld.naive(i).tunedAngle(matchInds) - tuneNew.naive(i).tunedAngle(matchInds))==1))/length(tuneOld.naive(i).touchID);
    diffs(i,2) = length(find(abs(tuneOld.naive(i).unimodalSingle(matchInds) - tuneNew.naive(i).unimodalSingle(matchInds))==1))/length(tuneOld.naive(i).touchID);
    diffs(i,3) = length(find(abs(tuneOld.naive(i).unimodalBroad(matchInds) - tuneNew.naive(i).unimodalBroad(matchInds))==1))/length(tuneOld.naive(i).touchID);
    diffs(i,4) = length( intersect( setdiff( find(tuneOld.naive(i).multimodal(matchInds)), find(tuneOld.naive(i).unimodalBroad(matchInds)) ), find(tuneNew.naive(i).multimodal(matchInds)) ) )...
        /length(tuneOld.naive(i).touchID);
end

mean(diffs)
std(diffs)

%%
diffs = zeros(length(tuneNew.expert),4); % 1: tunedAngle, 2: unimodalSingle, 3: unimodalBroad, 4: multimodal
for i = 1 : length(tuneNew.expert)
    matchInds = intersect(find(tuneNew.expert(i).tuned), find(tuneOld.expert(i).tuned));
    diffs(i,1) = length(find(abs(tuneOld.expert(i).tunedAngle(matchInds) - tuneNew.expert(i).tunedAngle(matchInds))==1))/length(tuneOld.expert(i).touchID);
    diffs(i,2) = length(find(abs(tuneOld.expert(i).unimodalSingle(matchInds) - tuneNew.expert(i).unimodalSingle(matchInds))==1))/length(tuneOld.expert(i).touchID);
    diffs(i,3) = length(find(abs(tuneOld.expert(i).unimodalBroad(matchInds) - tuneNew.expert(i).unimodalBroad(matchInds))==1))/length(tuneOld.expert(i).touchID);
    diffs(i,4) = length( intersect( setdiff( find(tuneOld.expert(i).multimodal(matchInds)), find(tuneOld.expert(i).unimodalBroad(matchInds)) ), find(tuneNew.expert(i).multimodal(matchInds)) ) )...
        /length(tuneOld.expert(i).touchID);
end

mean(diffs)
std(diffs)