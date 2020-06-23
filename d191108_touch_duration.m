%% average touch duration.
%% in ms
%% across angles

baseDir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,39,52];
sessions = {[18,19],[7,10],[20,21],[16,17],[21,23],[20,21]}; 
touchDurations = cell(length(mice),1);
touchAngles = cell(length(mice),1);
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions{mi}(2);
    load(sprintf('%s%03d\\UberJK%03dS%02d', baseDir, mouse, mouse, session),'u')
    touchDurations{mi} = cell2mat(cellfun(@(x) x.protractionTouchDurationByWhisking, u.trials', 'un', 0));
    touchAngles{mi} = cell2mat(cellfun(@(x) x.angle*ones(1,length(x.protractionTouchDurationByWhisking)), u.trials', 'un', 0));
end

%%

mean(cellfun(@mean, touchDurations))
sem(cellfun(@mean, touchDurations))