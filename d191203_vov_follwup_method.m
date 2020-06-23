% %% For illustrating how to follow-up same volume-of-view across days
% %% Using z-stack images as a sample.
% 
load('D:\TPM\JK\suite2p\039\zstack_039_998.mat', 'knobbyInfo', 'zstack')
load('D:\TPM\JK\suite2p\039\regops_039_001.mat', 'ops1')
% %cosd(35)
% 
% 
%%
% 
implay(zstack(ops1{1}.useY, ops1{1}.useX, :))
% 
%%
% 
figure, imshow(mat2gray(ops1{1}.mimg1))
% 
% 
%%
% 
zComp = zstack(ops1{1}.useY, ops1{1}.useX, :);
% 
template = ops1{1}.mimg1;

%%

repeat = 51;
interval = 1; % 1 is equal to 2 um
startPlane = 200+interval*(repeat-1)/2;
corval = zeros(1,repeat);
figure
% imshpr = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane - interval * i;
    temp = zComp(:,:,pi);
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
    imshowpair(mat2gray(tempRef), mat2gray(tempIm))
    
    drawnow;
    F = getframe(gcf);
    imshpr{i} = frame2im(F);
end

sizes = cell2mat(cellfun(@size, imshpr, 'un', 0));
yLength = min(sizes(:,1));
xLength = min(sizes(:,2));
%%
videoFn = 'vid_imshowpair.mp4';
v = VideoWriter(videoFn, 'MPEG-4');
v.FrameRate = 1;
v.Quality = 100;
open(v)
for i = 1 : length(imshpr)
    writeVideo(v,imshpr{i}(1:yLength, 1:xLength,:))
end
close(v)
%%
figure
implot = cell(repeat,1);
for i = 1 : repeat
    plot(-interval*(repeat-1):interval*2:-interval*(repeat-1)+interval*(i-1)*2, corval(1:i), 'k.')
    xlim([-interval*repeat, interval*repeat])
    ylim([0.78 0.86])
    xlabel('Depth (um)'), ylabel('Spatial correlation')
    set(gca,'box','off')
    drawnow
    
    F = getframe(gcf);
    implot{i} = frame2im(F);
end

videoFn = 'vid_corrPlot.mp4';
v = VideoWriter(videoFn, 'MPEG-4');
v.FrameRate = 1;
v.Quality = 100;
open(v)
for i = 1 : length(implot)
    writeVideo(v,implot{i})
end
close(v)























%%

baseDir = 'D:\TPM\JK\suite2p\';
mouse = 39;
refSession = 1;
zSession = 998;

load(sprintf('%s%03d\\zstack_%03d_%03d', baseDir, mouse, mouse, zSession), 'knobbyInfo', 'zstack')
load(sprintf('%s%03d\\regops_%03d_%03d', baseDir, mouse, mouse, refSession), 'ops1')

zComp = zstack(ops1{1}.useY, ops1{1}.useX, :);
template = ops1{1}.mimg1;
%%
repeat = 51;
interval = 1; % 1 is equal to 2 um
refPlane = 200;
startPlane = refPlane - interval * (repeat-1)/2;
corval = zeros(1,repeat);
zCompIms = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane + interval * (i-1);
    temp = zComp(:,:,pi);
    zCompIms{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
end
%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')

%%
figure, imshow(mat2gray(template))
hold on
pixResolution = 1.4 / str2double(ops1{1}.info.config.magnification_list(ops1{1}.info.config.magnification,:));
scaleBarPix = 100 / pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(template);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)

%%

saveDir = 'C:\Users\shires\Dropbox\Works\Presentation\materials_angle_S1\';
fn = 'templateFOV.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
fix_eps_fonts([saveDir, fn])

%%
saveDir = 'C:\Users\shires\Dropbox\Works\Presentation\materials_angle_S1\';
fnBase = '10umIntervalIm';
figure, hold on
for i = 1 : repeat
    imshow(mat2gray(zCompIms{i})), hold on
    plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
    fn = sprintf('%s%02d.eps',fnBase,i);
    export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
    fix_eps_fonts([saveDir, fn])
end



%%
[~,maxi] = max(corval);
refPlane2 = startPlane + interval * (maxi-1);
repeat2 = 11;
interval2 = 1; 
startPlane2 = startPlaneMax - interval2 * (repeat2-1)/2;

corval2 = zeros(1,repeat2);
zCompIms2 = cell(repeat2,1);
for i = 1 : repeat2
    pi = startPlane2 + interval2 * (i-1);
    temp = zComp(:,:,pi);
    zCompIms2{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval2(i) = temprho(1,2); 
end
%%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20), hold on
plot([knobbyInfo(startPlane2:interval2:startPlane2+interval2*(repeat2-1)).z]/cosd(35), corval2(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')
%%
saveDir = 'C:\Users\shires\Dropbox\Works\Presentation\materials_angle_S1\';
fnBase = '2umIntervalIm';
figure, hold on
for i = 1 : repeat
    imshow(mat2gray(zCompIms{i})), hold on
    plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
    fn = sprintf('%s%02d.eps',fnBase,i);
    export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
    fix_eps_fonts([saveDir, fn])
end














%% Comparison of the correlation graph between tdTomato and GCaMP6s

make_zstack('054_995');
%%
baseDir = 'D:\TPM\JK\2p\';
mouse = 54;
refSession = 1;
zSession = 995;

load(sprintf('%s%03d\\zstack_%03d_%03d', baseDir, mouse, mouse, zSession), 'knobbyInfo', 'zstack')
load(sprintf('%s%03d\\regops_%03d_%03d', baseDir, mouse, mouse, refSession), 'ops1')

%%
zComp = squeeze(zstack(ops1{1}.useY, ops1{1}.useX, :, 1)); % 1 for GREEN, 2 for RED
% %%
% refPlane = 145;
% figure, imshow(mat2gray(zComp(:,:,refPlane)))

template = ops1{1}.mimg1; % mimg1 for GREEN, mimgRED for RED
% template = ops1{1}.mimgRED;
% figure, imshow(mat2gray(template));
% %%
% template = ops1{1}.mimgRED;

repeat = 51;
interval = 1; % 1 is equal to 2 um
refPlane = 180;
startPlane = refPlane - interval * (repeat-1)/2;
corval = zeros(1,repeat);
zCompIms = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane + interval * (i-1);
    temp = zComp(:,:,pi);
    zCompIms{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
end
%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
% ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')




%%
zComp = squeeze(zstack(ops1{1}.useY, ops1{1}.useX, :, 2)); % 1 for GREEN, 2 for RED
% %%
% refPlane = 145;
% figure, imshow(mat2gray(zComp(:,:,refPlane)))

% template = ops1{1}.mimg1; % mimg1 for GREEN, mimgRED for RED
template = ops1{1}.mimgRED;
% figure, imshow(mat2gray(template));
% %%
% template = ops1{1}.mimgRED;

repeat = 51;
interval = 1; % 1 is equal to 2 um
refPlane = 180;
startPlane = refPlane - interval * (repeat-1)/2;
corval = zeros(1,repeat);
zCompIms = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane + interval * (i-1);
    temp = zComp(:,:,pi);
    zCompIms{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
end
%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
% ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')


%%

figure, imshow(mat2gray(zComp(:,:,180)))