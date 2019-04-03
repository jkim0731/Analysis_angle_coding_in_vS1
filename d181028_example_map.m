%%
% load u and anova results first
%%
P = [0,0,143;0,47,255;0,223,255;143,255,111;255,207,0;255,31,0;127,0,0]; % for 7 angles
% P = [0,0,143;0,143,255;143,255,111;255,111,0;127,0,0]; % for 5 distances
%% (2) linking lines of each distance across depth

scalebar = 50; % in micron
marginPix = 20;

figure, 
for i = 1 : length(u.cellmap)

%     if i < 3
%         position = i;
%     elseif i < 5
%         position = i + 2;
%     elseif i < 7
%         position = i - 2;
%     else
%         position = i;
%     end
%     subplot(2, 4, position), imshow(mat2gray(u.mimg{i})+0.5), hold on

    subplot(2, 2, i), imshow(mat2gray(u.mimg{i})+0.5), hold on

    
    
    for j = 1 : length(angles)
        currColor = P(j,:)/255;
        currCells = cellsTuned(tuneAngle == angles(j));
        currCells = currCells(find(currCells > i * 1000 & currCells <= (i+1) * 1000));
        for ci = 1 : length(currCells)
            cellInd = find(u.cellNums == currCells(ci));
            tuneInd = find(cellsTuned == currCells(ci));
            currSharp = tuneSharpness(tuneInd);
            if currSharp == 0
                currSharp = 1;
            end
%             plot(u.cellx(cellInd), u.celly(cellInd), 'marker', 'o', 'markersize', min(5 + tuneModulationMaxmin(currInds(ci))/0.3, 10), 'markerfacecolor', currColor, 'markeredgecolor', 'none')
            scatter(u.cellx(cellInd), u.celly(cellInd), min(50 + tuneModulationMaxmin(tuneInd)*30, 120), currColor, 'filled', 'markerfacealpha', 0.2 + (length(angles) + 1 - currSharp)/length(angles) * 0.8)
        end
    end
    
    currCells = cellsNTResponse;
    currCells = currCells(find(currCells > i * 1000 & currCells <= (i+1) * 1000));
    for ci = 1 : length(currCells)
        cellInd = find(u.cellNums == currCells(ci));
        scatter(u.cellx(cellInd), u.celly(cellInd), 30, [0 0 0], 'filled')
    end
    
    plot(u.c2xpoints, u.c2ypoints, 'k--', 'linewidth', 5)
%     sbendpoint = flip(size(u.mimg{i})) - marginPix;
%     sbstartpoint = [sbendpoint(1) - scalebar/u.pixResolution, sbendpoint(2)];    
%     plot([sbendpoint(1), sbstartpoint(1)], [sbendpoint(2), sbstartpoint(2)], 'w-', 'linewidth', 5)
%     depth = [num2str(-round(mean(mean(u.fovdepth{i})))), ' \mum'];
%     text(size(u.mimg{i},2) - 100, 30, depth, 'color', 'black', 'fontsize', 15);
end
