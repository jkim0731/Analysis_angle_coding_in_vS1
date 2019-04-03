% clear('*Array','b')

clim = [];

h = figure('Units','normalized', 'Position', [0 0.2 0.5 0.3]);
% figure

% % Select one among these
% showCellNums = u.cellNums;
% showCellNums = cellsTuned(tuneSharpness==6);
showCellNums = cellsNTResponse;
% showCellNums = tuneBroad;

ci = 1;
while true
% for i = 1 : length(u.cellNums)
    
    cellid = ci;
%     cellNum = u.cellNums(cellid);
    if cellid > length(showCellNums)
        cellid = mod(cellid,length(showCellNums));
    end        
    cellNum = showCellNums(cellid);
        
    dF = dFtotal{(cellsTotal == cellNum)};
    dfmean = cell2mat(cellfun(@(x) nanmean(x), dF, 'uniformoutput', false));
    %%
    % ANOVA
    touchNumGroups = cellfun(@(x) size(x,1), dF);
    timeAveragedF = zeros(sum(touchNumGroups),1);
    anovaGroupsF = zeros(sum(cellfun(@(x) size(x,1), dF)),1);
    for ai = 1 : length(angles)
        timeAveragedF( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = mean(dF{ai}(:,baseFrameNum+1:end),2);
        anovaGroupsF( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = deal(ai);        
    end
    
    %% ANOVA
    tTestH = cellfun(@(x) ttest(nanmean(x(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2)), dF);
    tTestInd = find(tTestH);
    means = cellfun(@(x) mean(nanmean(x(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2)), dF);
    sems = cell2mat(cellfun(@(x) std(nanmean(x(:,baseFrameNum+1:baseFrameNum+afterFrameNum),2))/sqrt(size(x,1)), dF, 'uniformoutput',false));

    %%    
    
    if isempty(clim)
        subplot(121), imagesc(dfmean);
    else
        subplot(121), imagesc(dfmean, clim);
    end
    hold on
    for i = 1 : size(dF,1)
        message = ['n=',num2str(size(dF{i},1))];
        text(5.3,i, message,'FontWeight', 'bold','FontSize',17,'Color','white');
    end
    c = colorbar;
%     c.Label.String = '\Delta(\Delta F / F_0)';
%     c.Label.FontSize = 25;
%     c.Label.FontWeight = 'bold';
%     set(gca,'XTick',1:2:baseFrameNum + afterFrameNum)
%     set(gca,'XTickLabel',round([-baseFrameNum+1:2:afterFrameNum]*1/frameRate*100)/100)
    set(gca,'XTick',2:2:baseFrameNum + afterFrameNum)
    set(gca,'XTickLabel',round([-baseFrameNum+2:2:afterFrameNum]*1/frameRate*100)/100)
    xlim([size(dfmean,2)-afterFrameNum-3, size(dfmean,2)])
    xlabel('Time after touch (s)')
    set(gca,'YTickLabel', angles, 'fontsize',17, 'fontweight', 'bold')
    yticks(1:7)
    ylabel('Pole angle (\circ)')
    title(['cell ID = ', num2str(cellid), ', cell Number = ', num2str(cellNum)])
    set(gca, 'fontweight', 'bold', 'fontsize', 17, 'linewidth', 5, 'box', 'off')
    hold off
    
    subplot(122)
    errorbar(angles, means, sems, 'k-', 'linewidth', 5, 'capsize', 10), hold on    
%     plot(angles(tTestInd), means(tTestInd) + sems(tTestInd)*1.5, 'b*')
%     plot(angles, ones(length(angles),1) * cellThresholds(cellsTotal == cellNum,1), 'b--', 'linewidth', 1)
%     plot(angles, ones(length(angles),1) * cellThresholds(cellsTotal == cellNum,2), 'b--', 'linewidth', 1)
    
%     [~, maxResInd] = max(abs(anovaStat.means));
%     ind__1 = find(pairComp(:,1) == maxResInd);
%     ind__2 = find(pairComp(:,2) == maxResInd);
%     testInd = union(ind__1, ind__2);
%     insigInd = find(pairComp(testInd,6) >= categoryThreshold);
%     sigInd = find(pairComp(testInd,6) < categoryThreshold);
%     temp = pairComp(testInd(insigInd),1:2);
%     insigIndGroup = unique(temp(:)); % sorted. Include tunedAngleInd.
%     temp = pairComp(testInd(sigInd),1:2);
%     sigIndGroup = setdiff(temp(:), maxResInd);
%     errorbar(angles(maxResInd), meanNsem(maxResInd,1), meanNsem(maxResInd,2), 'r.')
%     errorbar(angles(sigIndGroup), meanNsem(sigIndGroup,1), meanNsem(sigIndGroup,2), 'y.')
    
    set(gca,'XTick', angles, 'fontsize',17), xlabel('Pole angle (\circ)'), ylabel('\Delta(\Delta F / F_0)')
%     title({['ANOVA p = ', num2str(anovaP(cellsTotal == cellNum))]})
    set(gca, 'fontweight', 'bold', 'fontsize', 17, 'linewidth', 5, 'box', 'off')
    hold off
    

%%
% %     subplot(223), 
%     if ismember(cellNum, cellsTuned)
%         aind = find(angles == tuneAngle(find(cellsTuned==cellNum)));
%         temp = dF{aind};
%         temp = temp - mean(temp(:,1:baseFrameNum),2);
%         temp = temp(:,baseFrameNum-3:end);
%         for i = 1 : size(temp,1)
%             if mean(temp(i,3:end)) > 0.2
%                 plot(1:size(temp,2), temp(i,:), 'color', [0.7 0.7 0.7], 'linewidth', 3), hold on 
%             end
%         end
%         temp = temp(find(mean(temp(:,3:end))) > 0.2,:);
%         errorbar(1:size(temp,2), mean(temp), std(temp)/sqrt(size(temp,1)), '-', 'color', [0 114 189]/255, 'linewidth', 7, 'capsize', 12)
%     else
%         temp = cell2mat(cellfun(@(x) x, dF, 'uniformoutput', false));
%         temp = temp - mean(temp(:,1:baseFrameNum),2);
%         temp = temp(:,baseFrameNum-3:end);
%         c = jet(7);        
%         for i = 1 : size(temp,1)
%             if mean(temp(i,3:end)) > 0.2
%                 plot(1:size(temp,2), temp(i,:), 'color', [0.7 0.7 0.7]), hold on 
% %             plot(1:size(temp,2), temp(i,:), 'color', c(anovaGroupsF(i),:), 'linewidth', 3), hold on 
%             end
%         end
%         temp = temp(find(mean(temp(:,3:end)) > 0.2),:);
%         errorbar(1:size(temp,2), mean(temp), std(temp)/sqrt(size(temp,1)), '-', 'color', [217 83 25]/255, 'linewidth', 7, 'capsize', 12)
%     end
%     hold off
%     set(gca,'XTick',1:2:baseFrameNum + afterFrameNum)
% %     set(gca,'XTickLabel',round([-baseFrameNum+1:afterFrameNum]*1/frameRate*100)/100)
%     set(gca,'XTickLabel',round([-baseFrameNum+1+5:2:afterFrameNum]*1/frameRate*100)/100)
%     set(gca, 'linewidth', 5, 'fontweight', 'bold', 'fontsize', 15, 'box', 'off')
%     xlabel('Time after touch (s)')    
%     ylabel('\Delta F / F_0')  
%     ylim([-2 8])
% %     ylim([min(mean(temp)) - max(std(temp))/sqrt(size(temp,1)) * 3, max(mean(temp)) + max(std(temp))/sqrt(size(temp,1)) * 3])
    
    while true
        if waitforbuttonpress
            value = double(get(gcf,'CurrentCharacter'));
            switch value
                case 28 % <-
                    if ci == 1
                        ci = length(showCellNums);
                    else
                        ci = ci - 1;
                    end
                case 29 % ->
                    if ci == length(showCellNums)
                        ci = 1;
                    else
                        ci = ci + 1;
                    end
                case 30 % up arrow
                    if ci > length(showCellNums) - 10
                        ci = ci + 10 - length(showCellNums);
                    else
                        ci = ci + 10;
                    end
                case 31 % down arrow
                    if ci < 10
                        ci = length(showCellNums) - 10 + ci;
                    else
                        ci = ci - 10;
                    end
                case 93 % ]
                    if ci > length(showCellNums) - 100
                        ci = ci + 100 - length(showCellNums);
                    else
                        ci = ci + 100;
                    end
%                     disp('Tuned')
%                     cdfthresucc = [cdfthresucc, sum(signalCDF)];
                case 91 % [
                    if ci < 100
                        ci = length(showCellNums) - 100 + ci;
                    else
                        ci = ci - 100;
                    end
%                     disp('Not tuned')
%                     cdfthrefail = [cdfthrefail, sum(signalCDF)];                    
                    
                    
                case 27 % when pressing esc                    
                    ci =0;
                
                    
            end
            break
        end
    end
    if ci == 0
        break
    end
end

