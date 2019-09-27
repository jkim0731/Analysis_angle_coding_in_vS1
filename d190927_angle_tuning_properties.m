%% 190927 Angle tuning properties

%% Purpose: Detailed description of angle tuning.
%     Modulation distribution
%     - How big is the modulation, if it’s angle-tuned?
%     - Is there a difference across tuned-angle?

%     Sharpness distribution
%     - How sharp is the response difference?
%     - How is it related to modulation, and tuned-angle?

%     % of response 
%     - How sparse is the response? 
%     - In both tuned-angle touches and all touches
 
%     Depth profile
%     - How are all these features distributed across depth?


% First, look at all naive. And then add expert and matching naive. 

%% basic settings and data loading

baseDir = 'D:\TPM\JK\suite2p\';
loadFn = 'angle_tuning_summary_predecision_NC';
load([baseDir, loadFn], 'naive', 'expert')


%% Modulation distribution

mi = 1;
naive(mi).
