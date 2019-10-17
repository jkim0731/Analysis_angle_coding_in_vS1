mouse = 39;
session = 1;
plane = 5;
load(sprintf('Y:\\Whiskernas\\JK\\suite2p\\%03d\\F_%03d_%03d_plane%d_proc_final_spikes_noiseCorrected',mouse,mouse,session,plane), 'dat')
load(sprintf('Y:\\Whiskernas\\JK\\suite2p\\%03d\\UberJK%03dS%02d',mouse,mouse,session), 'u')

%% load all the matching frames
frames = dat.ops.frame_to_use{plane};

fn = sprintf('H:\\%03d\\%03d_%03d_000',mouse,mouse,session);
z = jksbxreadframes(fn, frames);

%% settings for frames, # of spikes, and angle indices
imagingLayer = ceil((plane-4)/4)+1;
trialInds = u.planeTrialInds{imagingLayer};

frameNums = cellfun(@(x) size(x.spk,2), u.trials(trialInds));
% sum(frameNums) % total number of frames match with frame_to_use
trialFramesBoundary = [0; cumsum(frameNums)];
trialFrames = cell(length(trialInds),1);
for i = 1 : length(trialInds)
    trialFrames{i} = trialFramesBoundary(i) + 1 : trialFramesBoundary(i+1);
end


cinds = find(u.trials{trialInds(1)}.neuindSession > plane * 1000 & u.trials{trialInds(1)}.neuindSession < (plane+1) * 1000);
spkNums = cellfun(@(x) sum(sum(x.spk(cinds,:))), u.trials(trialInds));

angles = 45:15:135;
angleTrialInds = cell(length(angles),1);
for ai = 1 : length(angles)
    angleTrialInds{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(trialInds)));
end

%% Selecting frames
targetAngle = 45;
ai = find(angles == targetAngle);
spkNumsAngle = spkNums(angleTrialInds{ai});
[~, sorti] = sort(spkNumsAngle, 'descend');
targetTrialsIndSorted = angleTrialInds{ai}(sorti);

%%
targetOrder = 1:10;

targetFrames = cell2mat(trialFrames(targetTrialsIndSorted(targetOrder))');
% uTrialInd = trialInds(targetTrialsIndSorted(targetOrder));
% tpmTime = u.trials{uTrialInd}.tpmTime{mod(plane,5) + ceil((plane-4)/4)};
x = zeros([size(z,1), size(z,2), length(targetFrames)]);
% for i = 1 : size(z,3)
for i = 1 : length(targetFrames)
    D = zeros([size(x,1), size(x,2), 2]);
    fi = targetFrames(i);
    D(:,:,1) = dat.ops.DS(fi,2);
    D(:,:,2) = dat.ops.DS(fi,1);

    x(:,:,i) = imwarp(z(:,:,fi), D);
end

xx = x(50:end, 100:end,:);
implay(mat2gray(xx))

%% test filters and smoothing
medxx = zeros(size(xx));
for i = 1 : size(xx,3)
    medxx(:,:,i) = medfilt2(xx(:,:,i));
end
implay(mat2gray(medxx))
%% spatial filter
% avgxx = zeros(size(xx));
gausxx = zeros(size(xx));
% avgh = fspecial('average', [3 3]);
gaush = fspecial('gauss', [3 3], 0.5);
for i = 1 : size(xx,3)    
%     avgxx(:,:,i) = imfilter(xx(:,:,i), avgh);
    gausxx(:,:,i) = imfilter(xx(:,:,i), gaush);
end
% implay(mat2gray(avgxx))
% implay(mat2gray(gausxx))

%% temporal filter

temporalxx = smoothdata(gausxx, 3, 'movmean',5);
implay(mat2gray(temporalxx))

% %% write to video
% % with scale bar and time points
% videoFn = 'test.mp4';
% v = VideoWriter(videoFn, 'MPEG-4');
% v.FrameRate = u.frameRate;
% for i = 1 : size(xx,3)
%     frameIm = imfilt(xx(:,:,i), ;
%     


%% for 4 planes

mouse = 39;
session = 1;
planes = 5:8;
targetAngle = 45;
numTrials = 10; % top 10 # of spikes trials of the top plane

load(sprintf('Y:\\Whiskernas\\JK\\suite2p\\%03d\\UberJK%03dS%02d',mouse,mouse,session), 'u')
fn = sprintf('D:\\TPM\\JK\\%03d\\%03d_%03d_000',mouse,mouse,session);
% fn = sprintf('H:\\%03d\\%03d_%03d_000',mouse,mouse,session);
h = fspecial('gauss', [3 3], 0.5);
resultVid = cell(length(planes),1);

plane = planes(1); % top plane
imagingLayer = ceil((plane-4)/4)+1;
trialInds = u.planeTrialInds{imagingLayer};

frameNums = cellfun(@(x) size(x.spk,2), u.trials(trialInds));
% sum(frameNums) % total number of frames match with frame_to_use
trialFramesBoundary = [0; cumsum(frameNums)];
trialFrames = cell(length(trialInds),1);
for i = 1 : length(trialInds)
    trialFrames{i} = trialFramesBoundary(i) + 1 : trialFramesBoundary(i+1);
end

cinds = find(u.trials{trialInds(1)}.neuindSession > plane * 1000 & u.trials{trialInds(1)}.neuindSession < (plane+1) * 1000);
spkNums = cellfun(@(x) sum(sum(x.spk(cinds,:))), u.trials(trialInds));

angles = 45:15:135;
angleTrialInds = cell(length(angles),1);
for ai = 1 : length(angles)
    angleTrialInds{ai} = find(cellfun(@(x) x.angle == angles(ai), u.trials(trialInds)));
end

ai = find(angles == targetAngle);
spkNumsAngle = spkNums(angleTrialInds{ai});
[~, sorti] = sort(spkNumsAngle, 'descend');
targetTrialsIndSorted = angleTrialInds{ai}(sorti);

targetOrder = 1:numTrials;
targetFrames = cell2mat(trialFrames(targetTrialsIndSorted(targetOrder))');

poleAvailableFramesTrials = cell2mat(cellfun(@(x) (x.tpmTime{1} > min(x.poleUpTime) & x.tpmTime{1} < max(x.poleUpTime)), u.trials(trialInds(targetTrialsIndSorted(targetOrder))), 'un', 0)');
%%
z = cell(length(planes),1);
for pi = 1 : length(planes)
    
    fprintf('Loading plane #%d/%d\n', pi, length(planes))
    plane = planes(pi);
    load(sprintf('Y:\\Whiskernas\\JK\\suite2p\\%03d\\F_%03d_%03d_plane%d_proc_final_spikes_noiseCorrected',mouse,mouse,session,plane), 'dat')
    frames = dat.ops.frame_to_use{plane};
    readFrames = frames(targetFrames);
    z{pi} = jksbxreadframes(fn, readFrames);
end

for pi = 1 : length(planes)
    fprintf('Processing plane #%d/%d\n', pi, length(planes))
    x = zeros(size(z{pi}));
    plane = planes(pi);
    load(sprintf('Y:\\Whiskernas\\JK\\suite2p\\%03d\\F_%03d_%03d_plane%d_proc_final_spikes_noiseCorrected',mouse,mouse,session,plane), 'dat')
    for i = 1 : length(targetFrames)
        D = zeros([size(x,1), size(x,2), 2]);
        fi = targetFrames(i);
        D(:,:,1) = dat.ops.DS(fi,2);
        D(:,:,2) = dat.ops.DS(fi,1);

        tempx = imwarp(z{pi}(:,:,i), D);
        x(:,:,i) = imfilter(tempx,h);
    end

    resultVid{pi} = mat2gray(smoothdata(x(50:end, 150:end,:), 3, 'movmean', 5));    
end

wholeVidCell = cell(length(targetFrames),1);
for fi = 1 : length(targetFrames)
    I = zeros([size(resultVid{1},1),size(resultVid{1},2),length(planes)]);
    for pi = 1 : length(planes)
        I(:,:,pi) = resultVid{pi}(:,:,fi);
    end
    wholeVidCell{fi} = imtile(I, 'BorderSize', [10 10]);
end
%%
wholeVid = zeros([size(wholeVidCell{1}), length(wholeVidCell)]);
for fi = 1 : length(wholeVidCell)
    wholeVid(:,:,fi) = wholeVidCell{fi};
end
%%
scaleBarPix = 100 / u.pixResolution;
margins = 30; % 20 pixels each away from right bottom corner
sizes = size(wholeVidCell{1});
poleUpMarkX = sizes(2)-50;
poleUpMarkY = 50;
poleUpMarkSize = 30;

figure
vidIm = cell(length(wholeVidCell),1);
for i = 1 : length(wholeVidCell)
% for i = 10    
    imshow(wholeVidCell{i}, [-0.02 0.8]), hold on    
    plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
    if poleAvailableFramesTrials(i)
        plot(poleUpMarkX, poleUpMarkY, 'y.', 'markersize', poleUpMarkSize)
    end
    hold off
    drawnow
    F = getframe(gcf);
	vidIm{i} = frame2im(F);
end

%%

videoFn = 'vid_JK039S01_planes5t8_45degrees_10trials.mp4';
v = VideoWriter(videoFn, 'MPEG-4');
v.FrameRate = 2 * u.frameRate;
v.Quality = 100;
open(v)
for i = 1 : length(vidIm)
    writeVideo(v,vidIm{i}(30:855,89:976,:))
end
close(v)


%% time series traces
utrials = u.trials(trialInds(targetTrialsIndSorted(targetOrder)));
traces = cell(length(planes),1);
numCells = 8; % take the most active 8 cells in each plane
for pi = 1 : length(planes)
% pi = 1;
    plane = planes(pi);
    cind = find(utrials{1}.neuindSession > plane * 1000 & utrials{1}.neuindSession < (plane + 1) * 1000);
    stdCalcium = std(cell2mat(cellfun(@(x) x.dF(cind,:), utrials', 'un', 0)),[],2);
    [~,sorti] = sort(stdCalcium, 'descend');
    cindUse = cind(sorti(1:numCells));
    traces{pi} = cell2mat(cellfun(@(x) x.dF(cindUse,:), utrials', 'un', 0));
    for ci = 1 : size(traces{pi})
        traces{pi}(ci,:) = min_max_normalization(traces{pi}(ci,:));
    end
end
%%
colorListTemp = jet(7);
colorList = [colorListTemp([3,4],:); [1 0 1]; colorListTemp(6,:)]; % to remove yellow
scrollFrameStart = 50;
groupYmax = 50;
eachYmax = 5;
ylimMax = groupYmax * length(planes);
ylimMin = -3;
%%
stimPatchColor = [0.5 0.5 0];
timeScaleBarLengthInSec = 1; % in sec
timeScaleBarLength = timeScaleBarLengthInSec * u.frameRate;
timeScaleBarXmargin = 1 ; % how many frames from the right end
timeScaleBarY = -2;

numFrames = size(traces{1},2);
timeSeriesIm = cell(numFrames,1);

figure('units', 'norm', 'outerpos', [0.2 0 0.2 1]), 

for i = 1 : scrollFrameStart
% for i = 20
    plot(traces{length(planes)}(numCells,1:i)*eachYmax), hold on
    poleFrames = find(poleAvailableFramesTrials(1:scrollFrameStart));
    diffFrames = [0,find(diff(poleFrames) > 1), length(poleFrames)];
    numPatch = length(diffFrames)-1;
    for pi = 1 : numPatch
        patchX = [poleFrames(diffFrames(pi)+1), poleFrames(diffFrames(pi+1))];
        if patchX(1) < i
            patchX(2) = min(patchX(2), i);
            patch([patchX(1), patchX(1), patchX(2), patchX(2)], [ylimMin ylimMax ylimMax ylimMin], stimPatchColor, 'edgecolor', 'none')
        end
    end
    
    
    for pi = 1 : length(planes)
        for ci = 1 : numCells
            plot(traces{pi}(ci,1:i)*eachYmax + (length(planes)-pi) * groupYmax + (numCells - ci) * eachYmax, '-', 'color', colorList(pi,:), 'linewidth', 2)
        end
    end
    
    timeScaleBarX = scrollFrameStart - timeScaleBarXmargin - timeScaleBarLength;
    plot([timeScaleBarX, timeScaleBarX + timeScaleBarLength], [timeScaleBarY, timeScaleBarY], 'w-', 'linewidth', 5) 
    
    xticks([])
    yticks([])
    xlim([0 scrollFrameStart])
    ylim([ylimMin ylimMax])
    set(gca,'visible', 'on', 'Color','k')
    
    drawnow
    F = getframe(gcf);
	timeSeriesIm{i} = frame2im(F);
    hold off
end

for i = scrollFrameStart+1 : numFrames 
    plot(0,0), hold on
    poleFrames = find(poleAvailableFramesTrials);
    diffFrames = [0,find(diff(poleFrames) > 1), length(poleFrames)];
    numPatch = length(diffFrames)-1;
    for pi = 1 : numPatch
        patchX = [poleFrames(diffFrames(pi)+1), poleFrames(diffFrames(pi+1))];
        patch([patchX(1), patchX(1), patchX(2), patchX(2)], [ylimMin ylimMax ylimMax ylimMin], stimPatchColor, 'edgecolor', 'none')
    end

    for pi = 1 : length(planes)
        for ci = 1 : numCells
            plot(traces{pi}(ci,:)*eachYmax + (length(planes)-pi) * groupYmax + (numCells - ci) * eachYmax, '-', 'color', colorList(pi,:), 'linewidth', 2)
        end
    end
    
    timeScaleBarX = i-timeScaleBarXmargin - timeScaleBarLength;
    plot([timeScaleBarX, timeScaleBarX + timeScaleBarLength], [timeScaleBarY, timeScaleBarY], 'w-', 'linewidth', 5) 
    
    xlim([i-scrollFrameStart, i])
    ylim([ylimMin ylimMax])
    xticks([])
    yticks([])
    set(gca,'visible', 'on', 'Color','k')
    
    drawnow
    F = getframe(gcf);
	timeSeriesIm{i} = frame2im(F);
    hold off
end


videoFn = 'vid_timeSeries_JK039S01_planes5t8_45degrees_10trials.mp4';
v = VideoWriter(videoFn, 'MPEG-4');
v.FrameRate = 2 * u.frameRate;
v.Quality = 100;
open(v)
for i = 1 : length(timeSeriesIm)
    writeVideo(v,timeSeriesIm{i}(85:984,50:332,:))
end
close(v)

%% a giant plot of all cells (from one plane?)
plane = 5;
utrials = u.trials(trialInds(targetTrialsIndSorted(targetOrder)));
cind = find(utrials{1}.neuindSession > plane * 1000 & utrials{1}.neuindSession < (plane + 1) * 1000);
cind = cind(1:round(length(cind)/4));
allDF = cell2mat(cellfun(@(x) x.dF(cind,:), utrials', 'un', 0));

allDF = (allDF - min(allDF,[],2))./(max(allDF,[],2) - min(allDF,[],2))*3;
% allDF = allDF * 2;
yShift = [0:(size(allDF,1)-1)] *1;
allDF = allDF + ones(size(allDF)) .*yShift';
figure('units', 'norm', 'outerpos', [0 0 1 1]), hold on
poleFrames = find(poleAvailableFramesTrials);
diffFrames = [0,find(diff(poleFrames) > 1), length(poleFrames)];
numPatch = length(diffFrames)-1;
for pi = 1 : numPatch
    patchX = [poleFrames(diffFrames(pi)+1), poleFrames(diffFrames(pi+1))];
    patch([patchX(1), patchX(1), patchX(2), patchX(2)], [ylimMin ylimMax ylimMax ylimMin], stimPatchColor, 'edgecolor', 'none')
end
for i = 1 : length(cind)
    plot(allDF(i,:), 'color', colorList(1,:), 'linewidth', 1)
%     plot(allDF(i,:), 'k', 'linewidth', 1)
    
end

ylim([0 max(max(allDF))])
xlim([0 size(allDF,2)])
xticks([])
yticks([])
% axis off
set(gca,'color','k')

%%

saveDir = 'C:\Users\shires\Dropbox\Works\Presentation\materials_angle_S1\';
% fn = 'rawdF_oneForth_withStim_plane5.eps';
% fn = 'rawdF_oneForth_noSim_plane5.eps';
% fn = 'rawdF_oneForth_noStim_black_plane5.eps';
% fn = 'normdF_oneForth_noStim_black_plane5.eps';
% fn = 'normdF_oneForth_noSim_plane5.eps';
fn = 'normdF_oneForth_withStim_plane5.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
fix_eps_fonts([saveDir, fn])