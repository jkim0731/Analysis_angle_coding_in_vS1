% basic settings
tempMaskZ = polyval(ws.polyFitsMask{1}{1}, linspace(-0.3, 1.3));
zInterval = max(tempMaskZ) - min(tempMaskZ); % in pixels
numPoints = round(zInterval/ws.pxPerMm*33)+1; % total num points in defining the 3d mask is going to be numPoints^2
% 30 um interval.

% find the apex point of top-view mask (in y direction)
maskTopY = polyval(ws.polyFitsMask{1}{2}, linspace(-0.3, 1.3, numPoints));
% [~,maxIndTopView] = max(maskTopY);
maskTopX = polyval(ws.polyFitsMask{1}{1}, linspace(-0.3, 1.3, numPoints));
% apexTopView = [maskTopY(maxIndTopView), maskTopX(maxIndTopView)];

% maskTopY and maskTopX is a reference line at the top-down view (phi = 0).

% find the apex point of front-view mask (in y direction)

maskFrontY = polyval(ws.polyFitsMask{2}{2}, linspace(-0.3, 1.3, numPoints));
[~,maxIndFrontView] = max(maskFrontY);
maskFrontZ = polyval(ws.polyFitsMask{2}{1}, linspace(-0.3, 1.3, numPoints));
apexFrontView = [maskFrontY(maxIndFrontView), maskFrontZ(maxIndFrontView)];

% along the z axis, move the reference line following the relationship in
% front mask (maskFrontY and maskFrontz)
%%
mask3d = zeros(numPoints^2,3);
for i = 1 : numPoints
    mask3d((i-1)*numPoints+1 : i*numPoints,1) = maskTopX;
    mask3d((i-1)*numPoints+1 : i*numPoints,2) = maskTopY - (apexFrontView(1)-maskFrontY(i));
    mask3d((i-1)*numPoints+1 : i*numPoints,3) = maskFrontZ(i);
end


figure, plot3(mask3d(:,1), mask3d(:,2), mask3d(:,3), 'k.')
xlabel('x'), ylabel('y'), zlabel('z'), axis equal





%% finding the intersection point (whisker base)

testFrame = 500;
wpo = ws.whiskerPadOrigin;
vwidth = ws.imagePixelDimsXY(1);

R = [cosd(ws.mirrorAngle) -sind(ws.mirrorAngle); sind(ws.mirrorAngle) cosd(ws.mirrorAngle)]; % rotation matrix in top view
            
% find intersection points from both views
z = polyval(ws.polyFits{2}{1}(testFrame,:), linspace(0,1))';
w = polyval(ws.polyFits{2}{2}(testFrame,:), linspace(0,1))';

if sqrt(sum((wpo-[z(end) w(end)]).^2)) < sqrt(sum((wpo-[z(1) w(1)]).^2))
    % c(q_max) is closest to whisker pad origin, so reverse the (z,w) sequence
    z = z(end:-1:1);
    w = w(end:-1:1);
end
whiskerFront = [z'; w'];
maskFront = [polyval(ws.polyFitsMask{2}{1},linspace(-0.3,1.3)); polyval(ws.polyFitsMask{2}{2},linspace(-0.3,1.3))];
Pfront = Whisker.InterX(whiskerFront, maskFront);
if isempty(Pfront)
    z = polyval(ws.polyFits{2}{1}(testFrame,:), linspace(-0.1, 1))';
    w = polyval(ws.polyFits{2}{2}(testFrame,:), linspace(-0.1, 1))';
    if sqrt(sum((wpo-[z(end) w(end)]).^2)) < sqrt(sum((wpo-[z(1) w(1)]).^2))
        % c(q_max) is closest to whisker pad origin, so reverse the (z,w) sequence
        z = z(end:-1:1);
        w = w(end:-1:1);
    end
    whiskerFront = [z';w'];
    Pfront = Whisker.InterX(whiskerFront, maskFront);
end

if ~isempty(Pfront) % further calculation carried out only when there's a base point intersection in the front view
    Pfront = Pfront(:,end); % just in case where there is 2 intersection points.
    % top-view mask is picked from 3D mask, based on the z-axis intersection value.
    whiskerBaseZ = Pfront(1);
    [~,maskZind] = min(abs(maskFrontZ - whiskerBaseZ));
    maskZ = maskFrontZ(maskZind);
    mask3dInd = find(mask3d(:,3) == maskZ);
    maskTop = mask3d(mask3dInd,1:2);
    x = polyval(ws.polyFits{1}{1}(testFrame,:), linspace(0,1))';
    y = polyval(ws.polyFits{1}{2}(testFrame,:), linspace(0,1))';
    if sqrt(sum((wpo-[x(end) y(end)]).^2)) < sqrt(sum((wpo-[x(1) y(1)]).^2))
        % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
        x = x(end:-1:1);
        y = y(end:-1:1);
    end
    whiskerTop = [x'; y'];
    Ptop = Whisker.InterX(whiskerTop, maskTop');


    if isempty(Ptop)
        x = polyval(ws.polyFits{1}{1}(testFrame,:), linspace(-0.1, 1.1))';
        y = polyval(ws.polyFits{1}{2}(testFrame,:), linspace(-0.1, 1.1))';
        if sqrt(sum((wpo-[x(end) y(end)]).^2)) < sqrt(sum((wpo-[x(1) y(1)]).^2))
            % c(q_max) is closest to whisker pad origin, so reverse the (x,y) sequence
            x = x(end:-1:1);
            y = y(end:-1:1);
        end
        whiskerTop = [x'; y'];
        Ptop = Whisker.InterX(whiskerTop, maskTop');
    end

    if ~isempty(Ptop)
        Ptop = Ptop(:,end);

        tempData = NaN(length(x),3);
        distFromBase = zeros(length(x),1);
        tempData(:,1:2) = (R * [x' - Ptop(1); y' - Ptop(2)] + Ptop)'; % rotate in regard to the base (intersection between whisker and mask)
        for j = 1 : length(x)
            distFromBase(j) = Ptop(2) - y(j); % lateral distance from the mask, calculated from the top-view
            line = [1, vwidth; Pfront(2) - distFromBase(j), Pfront(2) - distFromBase(j)]; % corresponding line for front-view
            P = Whisker.InterX(line, whiskerFront);
            if size(P,2) > 0
                % projection to the axis orthogonal to the body axis. Height (j,3) does not have to change
                tempData(j,3) = P(1);
            end
        end

        [~,baseInd] = nanmin(abs(distFromBase));
        finiteInds = find(isfinite(sum(tempData,2)));
    
    %     [~, obj.baseInd(i)] = min(abs(finiteInds - baseInd));
    %     tempData = tempData(finiteInds,:);
    %     obj.base(i,:) = tempData(obj.baseInd(i),:);
    %     s = cumsum(sqrt([0; diff(tempData(:,1))].^2 + [0; diff(tempData(:,2))].^2) + [0; diff(tempData(:,3))].^2);
    %     q = s ./ max(s);
    % 
    %     if max(s) > obj.rInMm * obj.pxPerMm % 3D tracking should be at least rInMm long
    %         obj.trackerData{i} = tempData;
    %         px = polyfit(q',tempData(:,1)',obj.fitorder);
    %         py = polyfit(q',tempData(:,2)',obj.fitorder);
    %         pz = polyfit(q',tempData(:,3)',obj.fitorder);
    %         obj.fit3Data{i} = [px', py', pz'];
    %     end
    %     frameNum = round(ws.time{1}(tdtopind(i))/ws.framePeriodInSec + 1);
    %     if ~isempty(ws.whiskerPoleIntersection{frameNum,1}) && ~isempty(ws.whiskerPoleIntersection{frameNum,2})
    %         obj.intersectPoint(i,1:2) = (R * [ws.whiskerPoleIntersection{frameNum,1}(1) - Ptop(1); ws.whiskerPoleIntersection{frameNum,1}(2) - Ptop(2)] + Ptop)';
    %         obj.intersectPoint(i,3) = ws.whiskerPoleIntersection{frameNum,2}(1);
    %     end

        [~, tempInd] = min(abs(finiteInds - baseInd));
        tempData = tempData(finiteInds,:);
        objBase(1,:) = tempData(tempInd,:);
        s = cumsum(sqrt([0; diff(tempData(:,1))].^2 + [0; diff(tempData(:,2))].^2) + [0; diff(tempData(:,3))].^2);
        q = s ./ max(s);

        if max(s) > 3 * ws.pxPerMm % 3D tracking should be at least rInMm long        
            px = polyfit(q',tempData(:,1)',5);
            py = polyfit(q',tempData(:,2)',5);
            pz = polyfit(q',tempData(:,3)',5);        
        end
    end
end

%%

figure, plot3(mask3d(:,1), mask3d(:,2), mask3d(:,3), 'k.'), hold on
plot3(tempData(:,1), tempData(:,2), tempData(:,3), 'r.')
plot3(objBase(1,1), objBase(1,2), objBase(1,3), 'c.', 'markersize', 20)
xlabel('x'), ylabel('y'), zlabel('z'), axis equal


%% implemented in Whisker.Whisker3D_2pad.m