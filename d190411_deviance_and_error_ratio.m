% clear
% tic
% baseDir = 'C:\JK\';
% 
% mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]}; 
% 
% naiveInd = 1:length(mice);
% expertInd = find(cellfun(@length, sessions)==2);
% 
% 
% for ni = 1 : length(naiveInd)
%     mouse = mice(naiveInd(ni));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{naiveInd(ni)}(1);    
%     naive(ni) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% end
% 
% % expert = struct;
% for ei = 1 : length(expertInd)
%     mouse = mice(expertInd(ei));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{expertInd(ei)}(2);    
%     expert(ei) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% end
% 
% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
% % L4mice = [70];
% % L4sessions = [6];
% 
% % L4 = struct;
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% end
% 
% save('Y:\Whiskernas\JK\suite2p\glm_cell_function_error_ratio.mat', 'naive', 'expert', 'L4')
% toc
% %%
% glm = glm_results_cell_function_shuffling(25,4,'C:\JK\');
% 
% 
% load('cellFunctionRidgeDE010.mat')
% %%
% close all
% curr = naive(1);
% touchInd = find(ismember(curr.cellNums, curr.touchID));
% whiskingInd = find(ismember(curr.cellNums, curr.whiskingID));
% figure,
% for i = 1 : 5
%     subplot(2,3,i)
%     histogram(cellfun(@(x) x(i), glm.exclusionER), 1:0.02:1.5, 'normalization', 'probability')
%     xlabel('Error ratio from exclustion method')
% end
% 
% figure, 
% for i = 1 : 5
%     subplot(2,3,i)
%     histogram(cellfun(@(x) mean(x(i,:)), glm.errorRatio), 1:0.02:1.5, 'normalization', 'probability')
%     xlabel('Error ratio from exclustion method')
% end
% 
% figure,
% for i = 1 : 5
%     subplot(2,3,i)
%     a = cellfun(@(x) x(i), glm.exclusionER);
%     b = cellfun(@(x) mean(x(i,:)), glm.errorRatio);
%     plot(a,b,'k.'), xlabel('Exclusion error ratio'), ylabel('Mean of permutation error ratios')
%     title(sprintf('corr: %.4f', corr(a,b)))
%     
%     if i == 1
%         hold on
%         plot(a(touchInd), b(touchInd), 'r.')
%     elseif i == 4
%         hold on
%         plot(a(whiskingInd), b(whiskingInd), 'r.')
%     end
%     
% end
% 
% figure,
% for i = 1 : 5
%     subplot(2,3,i)
%     a = cellfun(@(x) x(i), glm.DEdiff);
%     b = cellfun(@(x) x(i), glm.exclusionER);    
%     plot(a,b,'k.'), xlabel('Deviance explained difference'), ylabel('Exclusion error ratio')
%     title(sprintf('corr: %.4f', corr(a,b)))
% end
% 
% figure,
% for i = 1 : 5
%     subplot(2,3,i)
%     a = cellfun(@(x) x(i), glm.DEdiff);
%     b = cellfun(@(x) mean(x(i,:)), glm.errorRatio);
%     plot(a,b,'k.'), xlabel('Deviance explained difference'), ylabel('Mean of permutation error ratios')
%     title(sprintf('corr: %.4f', corr(a,b)))
% end

% %%
% clear
% tic
% baseDir = 'C:\JK\';
% 
% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 
% 
% naiveInd = 1:length(mice);
% expertInd = find(cellfun(@length, sessions)==2);
% 
% 
% % for ni = 8 : length(naiveInd)
% %     mouse = mice(naiveInd(ni));
% %     cd(sprintf('%s%03d',baseDir,mouse))
% %     session = sessions{naiveInd(ni)}(1);    
% %     naive2(ni) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% % end
% 
% % expert2 = struct;
% for ei = 6 : length(expertInd)
%     mouse = mice(expertInd(ei));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{expertInd(ei)}(2);    
%     expert2(ei) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% end
% 
% % L4mice = [70,74,75,76];
% % L4sessions = [6,4,4,4];
% % L4mice = [70];
% % L4sessions = [6];
% 
% % L4 = struct;
% % for mi = 1 : length(L4mice)
% %     mouse = L4mice(mi);
% %     cd(sprintf('%s%03d',baseDir,mouse))
% %     session = L4sessions(mi);    
% %     L4(mi) = glm_results_cell_function_shuffling(mouse, session, baseDir);
% % end
% 
% save('Y:\Whiskernas\JK\suite2p\glm_cell_function_error_ratio_withWTV_shuffling_expert1.mat', 'expert')
% % toc

%
clear
tic
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

naiveInd = 1:length(mice);
expertInd = find(cellfun(@length, sessions)==2);

%%
% for ni = 1 : length(naiveInd)
%     mouse = mice(naiveInd(ni));
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = sessions{naiveInd(ni)}(1);    
%     naive(ni) = glm_results_cell_function_shuffling_WTV_ONLY(mouse, session, baseDir);
% end

% expert = struct;
for ei = 5 : length(expertInd)
    mouse = mice(expertInd(ei));
    cd(sprintf('%s%03d',baseDir,mouse))
    session = sessions{expertInd(ei)}(2);    
    expert(ei) = glm_results_cell_function_shuffling_WTV_ONLY(mouse, session, baseDir);
end

% L4mice = [70,74,75,76];
% L4sessions = [6,4,4,4];
% L4mice = [70];
% L4sessions = [6];

% L4 = struct;
% for mi = 1 : length(L4mice)
%     mouse = L4mice(mi);
%     cd(sprintf('%s%03d',baseDir,mouse))
%     session = L4sessions(mi);    
%     L4(mi) = glm_results_cell_function_shuffling_WTV_ONLY(mouse, session, baseDir);
% end

save('Y:\Whiskernas\JK\suite2p\glm_cell_function_error_ratio_WTV_ONLYlasso.mat', 'naive', 'expert')
toc