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
% cellID. a single number, or an array of numbers. cell indices.
% (if no cell ID is given, starts with 1)
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
    cellIDList = ca.touchID(varargin);
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

    % allValHeatMap = cell2mat(cellfun(@(x) cell2mat(x), egHeatMap, 'uniformoutput', false));
    climCa = [min(cellfun(@(x) min(min(cell2mat(x))), egHeatMapCa)), max(cellfun(@(x) max(max(cell2mat(x))), egHeatMapCa))];
    climSpk = [min(cellfun(@(x) min(min(cell2mat(x))), egHeatMapSpk)), max(cellfun(@(x) max(max(cell2mat(x))), egHeatMapSpk))];
    % clim = [prctile(allValHeatMap(:), 10), prctile(allValHeatMap(:), 90)];
    h1 = figure(1);
    h1.Units = 'normalized';
    h1.OuterPosition = [0.1, 0.1, 0.2, 0.8];
    for i = 1 : length(angles)
        subplot(length(angles),2,(i-1)*2+1)
        imagesc(cell2mat(egHeatMapCa{i}), climCa), colormap gray
        if i == length(angles)
            xticks([baselineFrames]);
            xticklabels({'0'})
        else
            xticks(0);
        end
        if i == 1
            title('Calcium (Z-score)')
        end
        subplot(length(angles),2,i*2)
        imagesc(cell2mat(egHeatMapSpk{i}), climSpk), colormap gray
        if i == length(angles)
            xticks([baselineFrames]);
            xticklabels({'0'})
        else
            xticks(0);
        end
        if i == 1
            title('Spike (#)')
        end
    end
    % imagesc(allValHeatMap, clim), colormap gray
    %
    h2 = figure(2);
    h2.Units = 'normalized';
    h2.OuterPosition = [0.3, 0.3, 0.6, 0.5];

    subplot(331)
    imagesc(sumHeatMapCa)
    xticks([]);
    yticks([]);
    ylabel('Calcium during pole up')
    colors = jet(length(angles));    
    
    subplot(332)
    hold off
    plot(-baselineFrames+1:afterFrames, sumHeatMapCa(1,:))
    hold on
    for i = 1 : length(angles)
        boundedline(-baselineFrames+1:afterFrames, sumHeatMapCa(i,:), std(sumTimeSeriesCa{i}) / sqrt(size(sumTimeSeriesCa{i},1)), 'cmap', colors(i,:), 'transparency', 0.2)        
    end
    for i = 1 : length(angles)
        plot(-baselineFrames+1:afterFrames, sumHeatMapCa(i,:), 'linewidth', 3, 'color', colors(i,:))
    end    
    xticks([]);
    ylabel('\DeltaF/F_0 Z-score')
    
    subplot(333)
    errorbar(angles, cellfun(@(x) mean(x), ca.val{caspkInd}), cellfun(@(x) std(x)/sqrt(length(x)), ca.val{caspkInd}), 'k-', 'linewidth', 3)
    xticks([]);
    
    subplot(334)
    imagesc(sumHeatMapSpk)
    colors = jet(length(angles));
    yticks([]);
    ylabel('Spikes during pole up')
    xticks([round(baselineFrames - u.frameRate/2), baselineFrames, round(baselineFrames + u.frameRate/2), round(baselineFrames + u.frameRate), round(baselineFrames + u.frameRate*3/2)])
    xticklabels({'-0.5', '0', '0.5', '1', '1.5'})
    xlabel('Time after pole onset (s)')
    
    subplot(335)
    hold off
    plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(1,:))
    hold on
    for i = 1 : length(angles)        
        boundedline(-baselineFrames+1:afterFrames, sumHeatMapSpk(i,:), std(sumTimeSeriesSpk{i}) / sqrt(size(sumTimeSeriesSpk{i},1)), 'cmap', colors(i,:), 'transparency', 0.2)
    end
    for i = 1 : length(angles)
        plot(-baselineFrames+1:afterFrames, sumHeatMapSpk(i,:), 'linewidth', 3, 'color', colors(i,:))
    end
    ylabel('\Delta #Spk')
    xticks([-3, 0, 3, 6, 9])
    xticklabels({'-0.5', '0', '0.5', '1', '1.5'})
    xlabel('Time after pole onset (s)')
    
    subplot(336)
    errorbar(angles, cellfun(@(x) mean(x), spk.val{caspkInd}), cellfun(@(x) std(x)/sqrt(length(x)), spk.val{caspkInd}), 'k-', 'linewidth', 3)
    xticks([]);
    title('# of spikes only during touch frames')
    
    subplot(339)
    tempcell = cellfun(@(x) mean(x(:,baselineFrames+1:end),2) - mean(x(:,1:baselineFrames),2), tempHMS, 'uniformoutput', false);
    errorbar(angles, cellfun(@(x) mean(x), tempcell), cellfun(@(x) std(x)/sqrt(length(x)), tempcell), 'r-', 'linewidth', 3)
    title('Time-averaged # of spikes')
    xticks([45:15:135])
    xlabel('Pole Angle \o')
    
    while true
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