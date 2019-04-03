clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
baseDir = 'D:\TPM\JK\suite2p\';
positiveOnly = 1;
negativeOnly = 0;
c2only = 1;
singleOnly = 1;
% range = [1, 136, 243, 350, 700];
range = [1, 136, 350, 700]; % for L2, L3, L4
layerGroup = {[2],[3],[4],[2:3],[2,4],[3,4],[2:4]};
% % layers = 2:4; % layers = 2; layers = 3; layers = 2:3; layers = 4;


permLength = 10000;

for mi = 1 : length(mice)
% for mi = 1
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        
        resultsfn = sprintf('JK%03dS%02dsingleCell_anova_FT_1secB_absthresh02_perm_sharp.mat', mouse, session);
        savefn = sprintf('JK%03dS%02dsingleCell_spatial_clustering_index_c2only.mat', mouse, session);


        targetDir = sprintf('%s%03d',baseDir,mouse);
%         ufn = sprintf('UberJK%03dS%02d.mat',mouse,session);

        cd(targetDir)
        load(resultsfn, 'cell*', 'tune*', 'angles')
%         load(ufn)

        clear cell

        if positiveOnly
            tuneIndsTotal = find(ismember(cellsTotal, cellsTuned(tuneDirection == 1)));
        else
            if negativeOnly
            else
                tuneIndsTotal = find(ismember(cellsTotal, cellsTuned));        
            end
        end

        
        sci = zeros(length(layerGroup),1);
        sciperm = cell(length(layerGroup),1);
        scipval = zeros(length(layerGroup),1);
        
        for li = 1 : length(layerGroup)
            layers = layerGroup{li};
        
            if length(layers) < 3
                currDepthInds = [];
                for i = 1 : length(layers)
                    currDepthInds = [currDepthInds, find(cellDepths > range(layers(i)-1) & cellDepths <= range(layers(i)))];
                end
                tuneInds = intersect(tuneIndsTotal, currDepthInds);                
            else 
                tuneInds = tuneIndsTotal;
            end
            
            if c2only                    
                tuneInds = intersect(tuneInds, find(cellisC2));
            end
            
            cellsTunedInds = find(ismember(cellsTuned, cellsTotal(tuneInds)));
            tuneAngles = tuneAngle(cellsTunedInds);
            withinGroupDist = cell(length(angles),1);
            betweenGroupDist = cell(length(angles)-1,1);
            for i = 1 : length(angles)
                currAngleInd = find(tuneAngles == angles(i));
                if length(currAngleInd) > 1
                    currAngleComb = combnk(currAngleInd,2);
                    withinGroupDist{i} = zeros(size(currAngleComb,1),1);
                    for j = 1 : length(withinGroupDist{i})
                        ind1 = tuneInds(currAngleComb(j,1));
                        ind2 = tuneInds(currAngleComb(j,2));
                        withinGroupDist{i}(j) = sqrt((celly(ind1) - celly(ind2))^2 + (cellx(ind1) - cellx(ind2))^2);
                    end
                else
                    withinGroupDist{i} = [];
                end
            end
            for i = 1 : length(angles)-1
                currAngleInd = find(tuneAngles == angles(i));
                if isempty(currAngleInd)
                    betweenGroupDist{i} = [];
                else
                    currOtherInd = [];
                    for j = i+1 : length(angles)
                        currOtherInd = [currOtherInd; find(tuneAngles == angles(j))];
                    end
                    if isempty(currOtherInd)
                        betweenGroupDist{i} = [];
                    else
                        betweenGroupDist{i} = zeros(length(currAngleInd) * length(currOtherInd),1);
                        for j = 1 : length(currAngleInd)
                            ind1 = tuneInds(currAngleInd(j));
                            for k = 1 : length(currOtherInd)            
                                ind2 = tuneInds(currOtherInd(k));
                                betweenGroupDist{i}((j-1)*length(currOtherInd) + k) = sqrt((celly(ind1) - celly(ind2))^2 + (cellx(ind1) - cellx(ind2))^2);
                            end
                        end
                    end
                end
            end

            withinGroupDistMean = mean(cell2mat(withinGroupDist));
            betweenGroupDistMean = mean(cell2mat(betweenGroupDist));

            sci(li) = betweenGroupDistMean / withinGroupDistMean;

            
            tempsciperm = zeros(permLength,1);
            parfor pi = 1 : permLength
                permTuneAngles = tuneAngles(randperm(length(tuneAngles)));

                pwithinGroupDist = cell(length(angles),1);
                pbetweenGroupDist = cell(length(angles)-1,1);
                for i = 1 : length(angles)
                    currAngleInd = find(permTuneAngles == angles(i));
                    if length(currAngleInd) > 1
                        currAngleComb = combnk(currAngleInd,2);
                        pwithinGroupDist{i} = zeros(size(currAngleComb,1),1);
                        for j = 1 : length(pwithinGroupDist{i})
                            ind1 = tuneInds(currAngleComb(j,1));
                            ind2 = tuneInds(currAngleComb(j,2));
                            pwithinGroupDist{i}(j) = sqrt((celly(ind1) - celly(ind2))^2 + (cellx(ind1) - cellx(ind2))^2);
                        end
                    else
                        pwithinGroupDist{i} = [];
                    end
                end
                for i = 1 : length(angles)-1
                    currAngleInd = find(permTuneAngles == angles(i));
                    if isempty(currAngleInd)
                        pbetweenGroupDist{i} = [];
                    else
                        currOtherInd = [];
                        for j = i+1 : length(angles)
                            currOtherInd = [currOtherInd; find(permTuneAngles == angles(j))];
                        end
                        if isempty(currOtherInd)
                            pbetweenGroupDist{i} = [];
                        else
                            pbetweenGroupDist{i} = zeros(length(currAngleInd) * length(currOtherInd),1);
                            for j = 1 : length(currAngleInd)
                                ind1 = tuneInds(currAngleInd(j));
                                for k = 1 : length(currOtherInd)            
                                    ind2 = tuneInds(currOtherInd(k));
                                    pbetweenGroupDist{i}((j-1)*length(currOtherInd) + k) = sqrt((celly(ind1) - celly(ind2))^2 + (cellx(ind1) - cellx(ind2))^2);
                                end
                            end
                        end
                    end
                end

                pwithinGroupDistMean = mean(cell2mat(pwithinGroupDist));
                pbetweenGroupDistMean = mean(cell2mat(pbetweenGroupDist));

                tempsciperm(pi) = pbetweenGroupDistMean / pwithinGroupDistMean;
            end
            sciperm{li} = tempsciperm;
            scipval(li) = length(find(sciperm{li} > sci(li))) / permLength;
        end

        save(savefn, 'sci*', 'range', 'layerGroup')
    end
end