clear
mice = [25,27,30,36,37,39,52,53,54,56];
expertmiceind = [1:4,6,7];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
baseDir = 'D:\TPM\JK\suite2p\';
positiveOnly = 1;
negativeOnly = 0;
% range = [1, 136, 243, 350, 700];
range = [1, 136, 350, 700]; % for L2, L3, L4
layerGroup = {[2],[3],[4],[2:3],[2,4],[3,4],[2:4]};
% % layers = 2:4; % layers = 2; layers = 3; layers = 2:3; layers = 4;

scinaive = zeros(length(mice),length(layerGroup));
scipvalnaive = zeros(length(mice),length(layerGroup));
sciexpert = zeros(length(expertmiceind), length(layerGroup));
scipvalexpert = zeros(length(expertmiceind), length(layerGroup));

for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        resultsfn = sprintf('JK%03dS%02dsingleCell_spatial_clustering_index_withdepth.mat', mouse, session);
        targetDir = sprintf('%s%03d',baseDir,mouse);
        cd(targetDir)
        load(resultsfn, 'sci', 'scipval')
        
        if si == 1 % naive
            scinaive(mi,:) = sci';
            scipvalnaive(mi,:) = scipval';
        else % expert
            eind = find(expertmiceind==mi);
            sciexpert(eind,:) = sci';
            scipvalexpert(eind,:) = scipval';
        end
    end
end