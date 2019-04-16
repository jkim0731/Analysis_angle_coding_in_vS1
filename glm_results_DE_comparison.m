function results = glm_results_DE_comparison(mouse, session, baseDir)

% For d190416_glmDEcomparison_results_summary
% Inputs:
%     mouse, session, baseDir
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

% basic setting
repeat = 10;
cd(sprintf('%s%03d',baseDir, mouse))
wtvFnBase = sprintf('glmWhiskerTouchVariablesONLY_JK%03dS%02d_R', mouse, session);
touchFnBase = sprintf('glmResponseType_JK%03dS%02d_m45_R', mouse, session);
tuningFn = sprintf('JK%03dS%02dangle_tuning',mouse,session);
load(tuningFn)

% basic output (from touch and tuning information)
results.cIDAll = spk.touchID;
results.tuned = spk.tuned;
results.tunedAngle = spk.tunedAngle;
complex = zeros(length(spk.touchID),1);
complex(union(union(find(spk.multimodal), find(spk.leaveOneOut)), find(spk.tuneDirection >= 2))) = 1;
results.complex = complex;
results.sharpness = spk.sharpness;
results.ramp = spk.ramp;

cIDAll = results.cIDAll;
touchDE = zeros(length(cIDAll), repeat);
whiskerDE = zeros(length(cIDAll), repeat);
touchWhiskerRatio = zeros(length(cIDAll), 1);
touchDependence = zeros(length(cIDAll),1);
parfor ci = 1 : length(cIDAll)
    tempTouch = zeros(1,repeat);
    tempWhisker = zeros(1,repeat);
    for ri = 1 : repeat
        touchFn = sprintf('%s%02d', touchFnBase, ri);
        wtvFn = sprintf('%s%02d', wtvFnBase, ri);
        touchDat = load(touchFn,'fitDevExplained','cIDAll');
        wtvDat = load(wtvFn,'fitDevExplained','cIDAll'); % for this one, cIDAll is just for safety
        tempTouch(ri) = touchDat.fitDevExplained(touchDat.cIDAll == cIDAll(ci));        
        tempWhisker(ri) = wtvDat.fitDevExplained(wtvDat.cIDAll == cIDAll(ci));        
    end
    tempTouch(~isfinite(tempTouch)) = 0; % very rare (1 case only, maybe)
    tempWhisker(~isfinite(tempWhisker)) = 0; % very rare (1 case only, maybe, in JK053 S03 116th entry)
    touchDE(ci,:) = tempTouch;    
    whiskerDE(ci,:) = tempWhisker;
    touchWhiskerRatio(ci) = mean(tempTouch) / mean(tempWhisker);
    if mean(tempTouch) > mean(tempWhisker)
        if ttest(tempTouch, tempWhisker)
            touchDependence(ci) = 1;
        end
    end
end
results.touchDE = touchDE;
results.whiskerDE = whiskerDE;
results.touchWhiskerRatio = touchWhiskerRatio;
results.touchDependence = touchDependence;
