%% proportion of tuned cells
%% in L2/3 C2, L2/3 non-C2, L4 C2, L4 non-C2
clear
baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]}; 

naiveInds = 1:12;
expertInds = [1:4,7,9];
scnnInds = 13:16;

for ni = 1 : length(naiveInds)
    mouse = mice(naiveInds(ni));
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{naiveInds(ni)}(1);
    load(sprintf('JK%03dS%02ddistance_tuning',mouse,session), 'spk')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
%         if ~strcmp(fieldnames{fi}, 'val')
            naive(ni).(fieldnames{fi}) = spk.(fieldnames{fi});
%         end
    end
end

for ei = 1 : length(expertInds)
    mouse = mice(expertInds(ei));
    cd(sprintf('%s%03d',baseDir, mouse))
    session = sessions{expertInds(ei)}(2);
    load(sprintf('JK%03dS%02ddistance_tuning',mouse,session), 'spk')
    fieldnames = fields(spk);
    for fi = 1 : length(fieldnames)
        
        expert(ei).(fieldnames{fi}) = spk.(fieldnames{fi});
        
    end
end

% for si = 1 : length(scnnInds)
%     mouse = mice(scnnInds(si));
%     cd(sprintf('%s%03d',baseDir, mouse))
%     session = sessions{scnnInds(si)}(1);
%     load(sprintf('JK%03dS%02dangle_tuning',mouse,session), 'spk')
%     fieldnames = fields(spk);
%     for fi = 1 : length(fieldnames)
%         
%         L4(si).(fieldnames{fi}) = spk.(fieldnames{fi});
%         
%     end
% end

cd(baseDir)
save('angle_tuning_summary.mat','naive','expert')