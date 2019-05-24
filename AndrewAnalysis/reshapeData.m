%% reshape data
tableData = [];
properties = [];
for roiID = 1:length(dataSet(animalID).touchID)
    roiName = dataSet(animalID).touchID(roiID);
    if (roiName > layers(1) && roiName < layers(2))
        responses = dataSet(animalID).val{roiID};

        angleVals = [];
        responseVals = [];
        for a = 1:length(responses)
            angleVals = [angleVals; repmat(angles(a), length(responses{a}), 1)];
            responseVals = [responseVals; responses{a}];
        end
        properties = [properties; [roiName, dataSet(animalID).tuned(roiID),...
                                   dataSet(animalID).tunedAngle(roiID),...
                                   dataSet(animalID).tuneDirection(roiID),...
                                   dataSet(animalID).sharpness(roiID),...
                                   dataSet(animalID).modulation(roiID)]];
        tableData = [tableData, responseVals];
    end
end