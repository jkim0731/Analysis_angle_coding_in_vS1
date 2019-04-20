%% distribution of tuning sharpness
%% distribution of tuning modulation level
%% and their change after learning


baseDir = 'D:\TPM\JK\suite2p\';
cd(baseDir)
tuning = load('angle_tuning_summary');
info = load('cellFunctionRidgeDE010.mat');

srange = 0:0.05:1.5;
mrange = 0:0.05:2;

i = 1;
sharpness = tuning.naive(i).sharpness(find(tuning.naive(i).tuned));
modulation = tuning.naive(i).modulation(find(tuning.naive(i).tuned));
tempShist = histcounts(sharpness, srange, 'normalization', 'probability');
tempMhist = histcounts(modulation, mrange, 'normalization', 'probability');
figure,
subplot(121), plot(srange(2:end), tempShist)
subplot(122), plot(mrange(2:end), tempMhist)




