%% first, figure out the tangent plane.
frame = 916;
x = w3.trackerData{frame}(w3.prePoint(frame):w3.postPoint(frame),1);
y = w3.trackerData{frame}(w3.prePoint(frame):w3.postPoint(frame),2);
z = w3.trackerData{frame}(w3.prePoint(frame):w3.postPoint(frame),3);

q = linspace(0,1, w3.postPoint(frame)-w3.prePoint(frame)+1);
px = polyfit(q',x,2);
py = polyfit(q',y,2);

pxDot = polyder(px);
pyDot = polyder(py);
xDot = polyval(pxDot,q);
yDot = polyval(pyDot,q);


midpoint = round((w3.postPoint(frame) - w3.prePoint(frame))/2);

%%
figure, hold on
plot3(x, y, z, 'k.')
for i = 1 : length(z)
    planeX = x(midpoint) + xDot(midpoint)*linspace(-1, 1);
    planeY = y(midpoint) + yDot(midpoint)*linspace(-1, 1);    
    planeZ = ones(length(planeX),1) * z(i);
    plot3(planeX, planeY, planeZ, 'r.');
end
axis equal
xlabel('x'), ylabel('y'), zlabel('z')

%% calculate the angle (alpha)
if strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'rightward')
    alpha = atand(xDot(midpoint) / yDot(midpoint));
    
elseif strcmp(ws.faceSideInImage,'top') && strcmp(ws.protractionDirection,'leftward')
    alpha = -atand(xDot(midpoint) / yDot(midpoint));
    
elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'downward')
    alpha = atand(yDot(midpoint) / xDot(midpoint));
    
elseif strcmp(ws.faceSideInImage,'left') && strcmp(ws.protractionDirection,'upward')
    alpha = -atand(yDot(midpoint) / xDot(midpoint));
    
elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'upward')
    alpha = atand(yDot(midpoint) / xDot(midpoint));

elseif strcmp(ws.faceSideInImage,'right') && strcmp(ws.protractionDirection,'downward')
    alpha = -atand(yDot(midpoint) / xDot(midpoint));
    
elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'rightward')
    alpha = -atand(xDot(midpoint) / yDot(midpoint));

elseif strcmp(ws.faceSideInImage,'bottom') && strcmp(ws.protractionDirection,'leftward')
    alpha = atand(xDot(midpoint) / yDot(midpoint));
end
%%
Rvertical = [cosd(-alpha) -sind(-alpha); sind(-alpha) cosd(-alpha)]; % rotation matrix in top view
refPoint = [x(midpoint); y(midpoint)];
projected = [(Rvertical * [x' - refPoint(1); y' - refPoint(2)] + refPoint)', z];
plot3(projected(:,1), projected(:,2), projected(:,3), 'b')


%% horizontal kappa should be the same


px = polyfit(q',projected(:,1),2);
py = polyfit(q',projected(:,2),2);
pz = polyfit(q',projected(:,3),2);

% 2. horizontal & vertical kappa
pxDot = polyder(px);
pxDoubleDot = polyder(pxDot);

pyDot = polyder(py);
pyDoubleDot = polyder(pyDot);

pzDot = polyder(pz);
pzDoubleDot = polyder(pzDot);

xDot = polyval(pxDot,q);
xDoubleDot = polyval(pxDoubleDot,q);

yDot = polyval(pyDot,q);
yDoubleDot = polyval(pyDoubleDot,q);

zDot = polyval(pzDot,q);
zDoubleDot = polyval(pzDoubleDot,q);

kappasH = (xDot.*yDoubleDot - yDot.*xDoubleDot) ./ ((xDot.^2 + yDot.^2).^(3/2)) * w3.pxPerMm; % SIGNED CURVATURE, in 1/mm.
kappasV = (zDot.*yDoubleDot - yDot.*zDoubleDot) ./ ((zDot.^2 + yDot.^2).^(3/2)) * w3.pxPerMm; % SIGNED CURVATURE, in 1/mm.

w3.kappaH(frame)
kappasH(midpoint)

w3.kappaV(frame)
kappasV(midpoint)



%%
tn = '217';
% tn = '218';

load([tn,'_WL_2pad.mat'])
load([tn,'_WST.mat'])
load([tn,'_WF_2pad.mat'])
w3 = Whisker.Whisker3D_2pad(ws);
wfnew = Whisker.WhiskerFinal_2pad(wl,w3);

%
poleFrames = union(wf.poleUpFrames, wf.poleMovingFrames);
ylim = [nanmin(union(wf.kappaV, wfnew.kappaV)), nanmax(union(wf.kappaV, wfnew.kappaV))];
figure('units', 'normalized', 'outerposition', [0 0.2 1 0.6]);
subplot(211)
hold on,
plot(wf.kappaV, 'k-')
plot(wfnew.kappaV, 'r-')
patch([poleFrames(1), poleFrames(1), poleFrames(end), poleFrames(end)], [ylim(1), ylim(2), ylim(2), ylim(1)], [0.9 0.9 0.9], 'linestyle', 'none')
for i = 1 : length(wf.protractionTFchunks)
    frames = wf.protractionTFchunks{i};
    patch([frames(1), frames(1), frames(end), frames(end)], [ylim(1), ylim(2), ylim(2), ylim(1)], [0.7 0.7 0.7], 'linestyle', 'none')
end
plot(wf.kappaV, 'k-')
plot(wfnew.kappaV, 'r-')

legend({'old', 'new'})
title('\kappa_V')
ylabel('Curvature (1/mm)')

subplot(212)
hold on
plot(wf.phi, 'k-')
plot(wfnew.phi, 'r.')
ylim = [nanmin(union(wf.phi, wfnew.phi)), nanmax(union(wf.phi, wfnew.phi))];
patch([poleFrames(1), poleFrames(1), poleFrames(end), poleFrames(end)], [ylim(1), ylim(2), ylim(2), ylim(1)], [0.9 0.9 0.9], 'linestyle', 'none')
for i = 1 : length(wf.protractionTFchunks)
    frames = wf.protractionTFchunks{i};
    patch([frames(1), frames(1), frames(end), frames(end)], [ylim(1), ylim(2), ylim(2), ylim(1)], [0.7 0.7 0.7], 'linestyle', 'none')
end
plot(wf.phi, 'k-')
plot(wfnew.phi, 'r.')
title('\phi')
ylabel('Angle (\circ)')
%%
