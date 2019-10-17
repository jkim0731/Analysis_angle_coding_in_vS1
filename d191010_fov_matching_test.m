baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,39,52];

mi = 2;
mouse = mice(mi);

load(sprintf('%s%03d\\cellIDmatchTest_JK%03d',baseDir,mouse,mouse))

for sessionInd = [1,3]
    figure, 
    numPlane = length(us.sessions(1).mimg);
    for planeInd = 1:numPlane
        if ~isempty(us.sessions(2).mimg{planeInd})
            subplot(3,3,planeInd)
            refIm = us.sessions(2).mimg{planeInd};
            movingIm = us.sessions(sessionInd).mimg{planeInd};
            movedIm = imwarp(movingIm, us.sessions(sessionInd).tform{planeInd}, 'OutputView', imref2d(size(refIm)));

            imshowpair(movedIm, refIm)
        end
    end
    
    figure, 
    for planeInd = 1:numPlane
        if ~isempty(us.sessions(2).mimg{planeInd})
            subplot(3,3,planeInd)
            refIm = us.sessions(2).cellmap{planeInd};
            movingIm = us.sessions(sessionInd).cellmap{planeInd};
            movedIm = imwarp(movingIm, us.sessions(sessionInd).tform{planeInd}, 'OutputView', imref2d(size(refIm)));

            imshowpair(movedIm, refIm)
        end
    end
end

%%


%%
figure
subplot(121), imshow(mat2gray(us.sessions(2).mimg{8}))
subplot(122), imshow(mat2gray(us.sessions(1).mimg{8}))

%%
planeI = 1;
sessionI = 3;
ref = mat2gray(us.sessions(2).mimg{planeI});
moving = mat2gray(us.sessions(sessionI).mimg{planeI});

%     ref2 = adapthisteq(ref);
%     moving2 = imhistmatch(moving, ref2);

    ref2 = adapthisteq(ref);
    moving2 = adapthisteq(moving);

tformTest = imregtform(moving2, ref2, 'rigid', optimizer, metric);
moved = imwarp(moving2, tformTest, 'OutputView', imref2d(size(ref2)));
figure, 
subplot(131), imshow(mat2gray(ref2))
subplot(132), imshow(mat2gray(moved))
subplot(133), imshowpair(ref2, moved)



%% Quantification of matched and non-matched cells (appearred, disappearred)

clear
baseDir = 'Y:\Whiskernas\JK\suite2p\';
mice = [25,27,30,36,39,52];
numMice = length(mice);

refi = 2;
moviList = [1,3];
    
similarity = cell(numMice,1);
app = cell(numMice,1);
disapp = cell(numMice,1);

for mi = 1 : numMice
    mouse = mice(mi);

    load(sprintf('%s%03d\\cellIDmatch_JK%03d',baseDir,mouse,mouse))

%     for i = 1 : length(us.sessions)
%         try
%             uber(i) = load(sprintf('%s%03d\\UberJK%03d%s_NC',baseDir,mouse,mouse,us.sessions(i).sessionName));
%         catch
%             uber(i) = load(sprintf('%s%03d\\UberJK%03d%s',baseDir,mouse,mouse,us.sessions(i).sessionName));
%         end
%     end
    numPlane = length(us.sessions(1).cellmap);
    
    similarity{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    app{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    disapp{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    
    for planei = 1 : numPlane        
        for movii = 1 : length(moviList)            
            movi = moviList(movii);
            cm = us.sessions(refi).cellmap{planei};
            if ~isempty(cm)
                refx(1) = min(find(sum(cm)));
                refx(2) = max(find(sum(cm)));
                refy(1) = min(find(sum(cm,2)));
                refy(2) = max(find(sum(cm,2)));
                refTemplate = zeros(size(us.sessions(refi).cellmap{planei}));
                refTemplate(refy(1):refy(2), refx(1):refx(2)) =1;

                cm = us.sessions(movi).cellmap{planei};
                movx(1) = min(find(sum(cm)));
                movx(2) = max(find(sum(cm)));
                movy(1) = min(find(sum(cm,2)));
                movy(2) = max(find(sum(cm,2)));

                movTemplate = zeros(size(us.sessions(movi).cellmap{planei}));
                movTemplate(movy(1):movy(2), movx(1):movx(2)) =1;

                % figure, imshowpair(refTemplate,movTemplate)
                movedTemplate = imwarp(movTemplate, us.sessions(movi).tform{planei}, 'outputview', imref2d(size(refTemplate)));
                movedTemplate(find(movedTemplate)) = 1;
                movedCellMapTemp = imwarp(us.sessions(movi).cellmap{planei}, us.sessions(movi).tform{planei}, 'outputview', imref2d(size(refTemplate)));
                movedCellMap = zeros(size(movedCellMapTemp));
                integerInds = find(movedCellMapTemp - floor(round(movedCellMapTemp)) == 0);
                nonZeroInds = find(movedCellMapTemp);
                finalInds = intersect(integerInds, nonZeroInds);
                movedCellMap(finalInds) = movedCellMapTemp(finalInds);
                % figure, imshowpair(refTemplate,movedTemplate)
                jointTerritory = refTemplate .* movedTemplate;
                refCells = setdiff(unique(us.sessions(refi).cellmap{planei} .* jointTerritory), 0);
                movCellTemp = setdiff(unique(movedCellMap .* jointTerritory), 0);
                id = unique(round(movCellTemp/1000));
                freq = zeros(length(id),1);
                for idi = 1 : length(id)
                    freq(idi) = length(find(floor(movCellTemp/1000) == id(idi)));
                end
                [~, maxind] = max(freq);
                maxid = id(maxind);

                inds = find(floor(movCellTemp/1000) == id(idi));
                movCells = movCellTemp(inds);

                movCi = find(ismember(us.sessions(movi).cellID, movCells));
                candCi = us.sessions(movi).matchedRefCellID(movCi);
                trackedCimov = movCi(find(candCi));
                trackedCIDmov = us.sessions(movi).cellID(trackedCimov);
                trackedCIDref = candCi(find(candCi));

                similarity{mi}(planei,movii) = length(find(candCi)) / (length(refCells) + length(movCells) - length(find(candCi)));
                app{mi}(planei,movii) = length(setdiff(refCells,trackedCIDref)) / length(refCells);
                disapp{mi}(planei,movii) = length(setdiff(movCells, trackedCIDmov)) / length(movCells);
            end
        end
    end
end


%% Distribution of the similarity

figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), similarity, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), similarity, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), similarity, 'un', 0)));

bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between before & after learning')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Similarity')
ylim([0 0.6])


figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), app, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), app, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), app, 'un', 0)));
bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between before & after learning')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Appearred rate')
ylim([0 0.6])


figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), disapp, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), disapp, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), disapp, 'un', 0)));
bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between before & after learning')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Disappearred rate')
ylim([0 0.6])


%% Between discrete angles and radial distances

figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,2)), similarity, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,2)), similarity, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,2)), similarity, 'un', 0)));

bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between discrete angles & radial distances')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Similarity')
ylim([0 0.6])


figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,2)), disapp, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,2)), disapp, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,2)), disapp, 'un', 0)));
bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between discrete angles & radial distances')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Appearred rate')
ylim([0 0.6])


figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,2)), app, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,2)), app, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,2)), app, 'un', 0)));
bar(1, nanmean(mat1), 'k')
bar(2, nanmean(mat2), 'k')
bar(3, nanmean(mat3), 'k')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
title('Between discrete angles & radial distances')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Disappearred rate')
ylim([0 0.6])

