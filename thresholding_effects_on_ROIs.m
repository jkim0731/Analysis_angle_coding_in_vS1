baseDir = 'D:\TPM\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1],[3],[3,21],[3],[3],[3]}; 

mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]}; 
borderCells = cell(length(mice),2); %(:,1) naive, (:,2) expert
for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : 2
        session = sessions{mi}(si);
        if mi == 2
            planes = 5:8;
        else
            planes = 1:8;
        end
        borderCells{mi,si} = zeros(length(planes),1);
        for pi = 1 : length(planes)
            plane = planes(pi);
            beforeFn = sprintf('%s%03d\\F_%03d_%03d_plane%d_proc',baseDir, mouse, mouse,session,plane);
            afterFn = sprintf('%s%03d\\F_%03d_%03d_plane%d_proc_final_spikes_noiseCorrected',baseDir, mouse, mouse,session,plane);
            beforeDat = load(beforeFn, 'dat');
            afterDat = load(afterFn, 'dat');
            b = [beforeDat.dat.stat.iscell];
            a = [afterDat.dat.stat.iscell];

            borderCells{mi,si}(pi) = sum(b) - sum(a);
        end
    end
end

%%
baseDir = 'D:\TPM\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1],[3],[3,21],[3],[3],[3]}; 

mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]}; 
changingIscell = cell(length(mice), 2); 
rollingWindowForBaseF = 100; % in s
baseFprctile = 5;

responseThreshold = 0.3; % 
responsePercentile = 95; % percentile over the threshold that assigns a cell

stdwindow = 5;
lowerprct = 5; % 5th percentile

for mi = 1:6
    mouse = mice(mi);
    cd([baseDir, sprintf('%03d',mice(mi))])
    for si = 1 : 2
        session = sessions{mi}(si);
        fnbase = sprintf('F_%03d_%03d_plane',mouse, session);
        fList = dir([fnbase,'*_proc.mat']);
        changingIscell{mi,si} = zeros(length(fList),1);
        for i = 1 : length(fList)
%         for i = 1
            fprintf('mouse %03d session %d\n', mice(mi), i)
            load(fList(i).name)
            dat.rollingWindowForBaseF = rollingWindowForBaseF;
            dat.baseFprctile = baseFprctile;
            dat.responseThreshold = responseThreshold;
            dat.responsePercentile = responsePercentile;

            inds = find([dat.stat.iscell]);
            npcoeffs = min(min(dat.Fcell{1}(inds,:) ./ dat.FcellNeu{1}(inds,:), [], 2), ones(length(inds),1)*0.7);
            npcoeffs = repmat(npcoeffs, 1, size(dat.Fcell{1},2));

            len = size(dat.Fcell{1},2);
            window = round(rollingWindowForBaseF*(dat.ops.imageRate/dat.ops.num_plane));

            tempF = dat.Fcell{1}(inds,:) - dat.FcellNeu{1}(inds,:) .* npcoeffs;
            tempFF = [tempF, tempF(:,end-window+1:end)];
            ttbase = zeros(length(inds),size(tempF,2));
            parfor k = 1 : len
                ttbase(:,k) = prctile(tempFF(:,k:k+window),baseFprctile,2);
            end        
            dF = (tempF - ttbase)./ttbase;

            % final sorting cell based on dF/F0
            response = prctile(dF, responsePercentile, 2);
            isNotCell = find(response < responseThreshold);
            changingIscell{mi,si}(i) = length(isNotCell);
        end
    end
end


%%

baseDir = 'D:\TPM\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1],[3],[3,21],[3],[3],[3]}; 

mice = [25,27,30,36,39,52];
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]}; 
changingIscell = cell(length(mice), 2); 
rollingWindowForBaseF = 100; % in s
baseFprctile = 5;

responseThreshold = 0.3; % 
responsePercentile = 95; % percentile over the threshold that assigns a cell

stdwindow = 5;
lowerprct = 5; % 5th percentile

for mi = 1:6
    mouse = mice(mi);
    cd([baseDir, sprintf('%03d',mice(mi))])
    for si = 1 : 2
        session = sessions{mi}(si);
        fnbase = sprintf('F_%03d_%03d_plane',mouse, session);
        fList = dir([fnbase,'*_proc_final_spikes_noiseCorrected.mat']);
        dFDist{mi,si} = cell(length(fList),1);        
        for i = 1 : length(fList)
%         for i = 1
            fprintf('mouse %03d session %d\n', mice(mi), i)
            load(fList(i).name, 'dat')
            dFDist{mi,si}{i} = sum(dat.dF>=0.3, 2)/size(dat.dF,2);
        end
    end
end
%%
dFDist{2,1} = dFDist{2,1}(5:8);
%%
histRange = 0:0.01:1;
dist = cell(6,2);
for mi = 1 : 6
    for si = 1 : 2
        dist{mi,si} = histcounts(cell2mat(dFDist{mi,si}),histRange, 'norm', 'prob');
    end
end
%%
figure, hold on
plot(histRange(2:end), mean(cell2mat(dist(:,1))), 'r')
plot(histRange(2:end), mean(cell2mat(dist(:,2))), 'b')
boundedline(histRange(2:end), mean(cell2mat(dist(:,1))), sem(cell2mat(dist(:,1))), 'r')
boundedline(histRange(2:end), mean(cell2mat(dist(:,2))), sem(cell2mat(dist(:,2))), 'b')
xlabel('Proportion dF/F0 > threshold (0.3)')
ylabel('Proportion')
legend({'Naive', 'Expert'})

%%




        
