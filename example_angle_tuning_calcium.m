% For figure 3b.
% Drawing example calcium responses from u.dF
% dF aligned with pole up onset, and 2 sec after
% Separated in different angles
% Summarizing heatmap (subtracting baseline)
% Summarizing time series in different angles (subtracting baseline)
% Can navigate using arrow keys
% exit with esc button
%

% dependency:
% boundedline.m (and its package. kakearney-boundedline-pkg)
%
% % input:
% u
% cellID. a single number, or an array of numbers. cell ID's.
% (if no cell ID is given, starts with 1001)
% No navigation option when only single number is given.
% 
% 
% % output:
% all calcium dF heatmap in different angles
% angle-averaged heatmap
% angle-averaged time series with sem

function example_angle_tuning_calcium(u, ca, spk, varargin)

afterDur = 2; % in sec

navigate = 1;
if nargin > 3
    cellIDList = varargin{1};
    if length(varargin{1}) == 1
        navigate = 0;
    end
else
    cellIDList = ca.touchID;
end

ci = 1;
while true 
%     close all
    cellID = cellIDList(ci);
    caspkInd = find(ca.touchID == cellID);
    angles = unique(cellfun(@(x) x.angle, u.trials));
    tind = find(cellfun(@(x) ismember(cellID, x.neuindSession), u.trials));
    cind = find(u.trials{tind(1)}.neuindSession == cellID);
    tpmInd = mod(floor(cellID/1000)-1,4)+1;
    baselineFrames = min(cellfun(@(x) min(cellfun(@(y) find(y>x.poleUpOnsetTime,1), x.tpmTime)), u.trials));
    afterFrames = round(afterDur * u.frameRate);

    tindGroups = cell(length(angles),1);
    egHeatMapCa = cell(length(angles),1);
    sumHeatMapCa = cell(length(angles),1);
    sumTimeSeriesCa = cell(length(angles),1);

    egHeatMapSpk = cell(length(angles),1);
    sumHeatMapSpk = cell(length(angles),1);
    sumTimeSeriesSpk = cell(length(angles),1);
    tempHMS = cell(length(angles),1);
    for ai = 1 : length(angles)
        tindGroups{ai} = intersect(tind, find(cellfun(@(x) x.angle == angles(ai), u.trials)));    
        egHeatMapCa{ai} = cellfun(@(x) x.dF(cind, find(x.tpmTime{tpmInd} > x.poleUpOnsetTime,1) - baselineFrames + 1 : find(x.tpmTime{tpmInd} > x.poleUpOnsetTime,1) + afterFrames), u.trials(tindGroups{ai}), 'uniformoutput', false);
        sumTimeSeriesCa{ai} = cell2mat(cellfun(@(x) x-mean(x(1:baselineFrames)), egHeatMapCa{ai}, 'uniformoutput', false));
        sumHeatMapCa{ai} = mean(sumTimeSeriesCa{ai});

        egHeatMapSpk{ai} = cellfun(@(x) x.spk(cind, find(x.tpmTime{tpmInd} > x.poleUpOnsetTime,1) - baselineFrames + 1 : find(x.tpmTime{tpmInd} > x.poleUpOnsetTime,1) + afterFrames), u.trials(tindGroups{ai}), 'uniformoutput', false);
        sumTimeSeriesSpk{ai} = cell2mat(cellfun(@(x) x-mean(x(1:baselineFrames)), egHeatMapSpk{ai}, 'uniformoutput', false));
        sumHeatMapSpk{ai} = mean(sumTimeSeriesSpk{ai});
        tempHMS{ai} = cell2mat(egHeatMapSpk{ai});
    end
    sumHeatMapCa = cell2mat(sumHeatMapCa);
    sumHeatMapSpk = cell2mat(sumHeatMapSpk);

    
        colors = jet(length(angles));
        
        
    
    % allValHeatMap = cell2mat(cellfun(@(x) cell2mat(x), egHeatMap, 'uniformoutput', false));
    climCa = [min(cellfun(@(x) min(min(cell2mat(x))), egHeatMapCa)), max(cellfun(@(x) max(max(cell2mat(x))), egHeatMapCa))];
    climSpk = [min(cellfun(@(x) min(min(cell2mat(x))), egHeatMapSpk)), max(cellfun(@(x) max(max(cell2mat(x))), egHeatMapSpk))];
    % clim = [prctile(allValHeatMap(:), 10), prctile(allValHeatMap(:), 90)];
    h1 = figure(1);
    h1.Units = 'normalized';
    h1.OuterPosition = [0.1, 0.1, 0.2, 0.6];
    for i = 1 : length(angles)
        subplot(length(angles),2,(i-1)*2+1)
        imagesc(cell2mat(egHeatMapCa{i}), climCa), colormap gray
        if i == length(angles)
            xticks([baselineFrames]);
            xticklabels({'^'})
        else
            xticks(0);
        end
        yticks([])
        ylabel([num2str(angles(i)), ' \circ'])
        if i == 1
            title('Calcium (Z-score)')
        end
        set(gca, 'fontweight', 'bold')
        subplot(length(angles),2,i*2)
        imagesc(cell2mat(egHeatMapSpk{i}), climSpk), colormap gray
        if i == length(angles)
            xticks([baselineFrames]);
            xticklabels({'^'})
        else
            xticks(0);
        end
        yticks([])
        if i == 1
            title('Spike (#)')
        end
        set(gca, 'fontweight', 'bold')
    end
    % imagesc(allValHeatMap, clim), colormap gray
    %
    h2 = figure(2);
    h2.Units = 'normalized';
    h2.OuterPosition = [0.3, 0.1, 0.3, 0.6];

    subplot(321)
    imagesc(sumHeatMapCa)
    xticks([]);
    yticks([1:7]);
    yticklabels({'45\circ', '60\circ', '75\circ', '90\circ', '105\circ', '120\circ', '135\circ'})
    title('Calcium during pole in')
    set(gca, 'fontweight', 'bold')
%     colorbar
    
    subplot(323)
    hold off
    plot(-baselineFrames+1:afterFrames, sumHeatMapCa(1,:))
    hold on
    for i = 1 : length(angles)
        boundedline(-baselineFrames+1:afterFrames, sumHeatMapCa(i,:), std(sumTimeSeriesCa{i}) / sqrt(size(sumTimeSeriesCa{i},1)), 'cmap', colors(i,:), 'transparency', 0.2)        
    end
    for i = 1 : length(angles)
        plot(-baselineFrames+1:afterFrames, sumHeatMapCa(i,:), 'linewidth', 3, 'color', colors(i,:))
    end    
    xlim([-baselineFrames+1 afterFrames])
%     xticks([round(-baselineFrames + u.frameRate/2), baselineFrames, round(baselineFrames + u.frameRate/2), round(baselineFrames + u.frameRate), round(baselineFrames + u.frameRate*3/2)])
    xticks([-3, 0, 3, 6, 9, 12])
    xticklabels({'-0.5', '0', '0.5', '1', '1.5', '2'})
    xlabel('Time after pole onset (s)')
    ylabel('\DeltaF/F_0 Z-score')
    set(gca, 'fontweight', 'bold')
    
    subplot(325)
    errorbar(angles, cellfun(@(x) mean(x), ca.val{caspkInd}), cellfun(@(x) std(x)/sqrt(length(x)), ca.val{caspkInd}), 'k-', 'linewidth', 3)
    ylabel('\DeltaF/F_0 Z-score')
    xticks([45:15:135])
    xlim([45 135])
    xlabel('Pole angle (\circ)')
    set(gca, 'fontweight', 'bold')
    
    subplot(322)
    imagesc(sumHeatMapSpk)
    colors = jet(length(angles));
    yticks([]);
    xticks([]);
    title('Spikes during pole in')
%     colorbar
    
    subplot(324)
    
    hold off
    % to remove previous plots, and for legend
    p1 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(1,:), 'linewidth', 3, 'color', colors(1,:));
    hold on
    p2 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(2,:), 'linewidth', 3, 'color', colors(2,:));
    p3 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(3,:), 'linewidth', 3, 'color', colors(3,:));
    p4 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(4,:), 'linewidth', 3, 'color', colors(4,:));
    p5 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(5,:), 'linewidth', 3, 'color', colors(5,:));
    p6 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(6,:), 'linewidth', 3, 'color', colors(6,:));
    p7 = plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(7,:), 'linewidth', 3, 'color', colors(7,:));


%     hold on
    for i = 1 : length(angles)        
        boundedline(-baselineFrames+1:afterFrames, sumHeatMapSpk(i,:), std(sumTimeSeriesSpk{i}) / sqrt(size(sumTimeSeriesSpk{i},1)), 'cmap', colors(i,:), 'transparency', 0.2)
    end
    for i = 1 : length(angles)
        plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(i,:), 'linewidth', 3, 'color', colors(i,:))
    end
    legend([p1,p2,p3,p4,p5,p6,p7], {'45\circ','60\circ','75\circ','90\circ','105\circ','120\circ','135\circ'})
    ylabel('\Delta #Spk')
    xlim([-baselineFrames+1 afterFrames])
    xticks([-3, 0, 3, 6, 9, 12])
    xticklabels({'-0.5', '0', '0.5', '1', '1.5', '2'})
    xlabel('Time after pole onset (s)')
    
    set(gca, 'fontweight', 'bold')
    
    subplot(326)
    errorbar(angles, cellfun(@(x) mean(x), spk.val{caspkInd}), cellfun(@(x) std(x)/sqrt(length(x)), spk.val{caspkInd}), 'k-', 'linewidth', 3)
    xticks([45:15:135])
    xlim([45 135])
    xlabel('Pole angle (\circ)')
    ylabel('\Delta #Spk')
    set(gca, 'fontweight', 'bold')
%     title('# of spikes only during touch frames')
    
%     subplot(428)
%     tempcell = cellfun(@(x) mean(x(:,baselineFrames+1:end),2) - mean(x(:,1:baselineFrames),2), tempHMS, 'uniformoutput', false);
%     errorbar(angles, cellfun(@(x) mean(x), tempcell), cellfun(@(x) std(x)/sqrt(length(x)), tempcell), 'r-', 'linewidth', 3)
%     title('Time-averaged # of spikes')
%     xticks([45:15:135])
%     xlim([45 135])
%     xlabel('Pole Angle \o')
    
    while true && navigate
        if waitforbuttonpress
            value = double(get(gcf,'CurrentCharacter'));
            switch value
                case 28 % <-
                    if ci == 1
                        ci = length(cellIDList);
                    else
                        ci = ci - 1;
                    end
                case 29 % ->
                    if ci == length(cellIDList)
                        ci = 1;
                    else
                        ci = ci + 1;
                    end
                case 30 % up arrow
                    if ci > length(cellIDList) - 10
                        ci = ci + 10 - length(showCellNums);
                    else
                        ci = ci + 10;
                    end
                case 31 % down arrow
                    if ci < 10
                        ci = length(cellIDList) - 10 + ci;
                    else
                        ci = ci - 10;
                    end
                case 93 % ]
                    if ci > length(cellIDList) - 100
                        ci = ci + 100 - length(cellIndList);
                    else
                        ci = ci + 100;
                    end
%                     disp('Tuned')
%                     cdfthresucc = [cdfthresucc, sum(signalCDF)];
                case 91 % [
                    if ci < 100
                        ci = length(cellIDList) - 100 + ci;
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
    if ci == 0 || navigate == 0
        break
    end
end