%% DE comparison between full model with only touch information and with only whisker touch variables
%% Only within touch cells
%% Gather 10 iterations and compare in statistical methods (t-test?)
%% summarize all mice together (16 mice, 12 naive, 6 learned, 4 L4 (among which only 2 are usable for now)

% Inputs:
%     mice, sessions
%
% Files:
%     glmWhiskerTouchVariablesONLY_JKxxxSxx_Rxx
%     glmResponseType_JKxxxSxx_m45_Rxx
%     (use cIDAll from glmWhiskerTouchVariablesONLY files)
%     JKxxxSxxangle_tuning
% 
% Outputs:
%     cIDAll
%     tuned
%     tuned angle
%     complex tuning: multimodal, tuneDirection >= 2, leaveOneOut
%     tuning sharpness
%     ramp tuning
%     touchDE (length(cIDAll), repeat)
%     whiskerDE (length(cIDAll), repeat)
%     touchVSwhiskerRatio: from the average
%     touchDependence: calculated by t-test between whiskerDE and touchDE. if touchDE is statistically larger than whiskerDE.

clear
tic
baseDir = 'D:\TPM\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);

for ni = 1 : length(naiveInd)
    mouse = mice(naiveInd(ni));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{naiveInd(ni)}(1);    
    naive(ni) = glm_results_DE_comparison(mouse, session, baseDir);
end


for ei = 1 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
    expert(ei) = glm_results_DE_comparison(mouse, session, baseDir);
end
% 
% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
% 
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_DE_comparison(mouse, session, baseDir);
% end

save('Y:\Whiskernas\JK\suite2p\glm_results_DEcomparison_lassoVSlasso.mat', 'naive', 'expert')
toc