mice = [25,27,30,36,37,38,39,41,52,53,54,56,70,74,75,76];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3],[6],[4],[4],[4]};

mi = 1;

baseDir = 'D:\TPM\JK\suite2p\';

mouse = mice(mi);
session = sessions{mi}(1);

load(sprintf('%s%03d\\UberJK%03dS%02d',baseDir, mouse, mouse, session), 'u');


figure, hold on
imshowpair(u.cellmap{3}, u.cellmap{8})

%%
figure, hold on
imshowpair(u.mimg{2}, u.mimg{7})

%%

mi = 4;
mouse = mice(mi);
session = sessions{mi}(1);
load(sprintf('%s%03d\\UberJK%03dS%02d',baseDir, mouse, mouse, session), 'u');
u1 = u;

mi = 16;
mouse = mice(mi);
session = sessions{mi}(1);
load(sprintf('%s%03d\\UberJK%03dS%02d',baseDir, mouse, mouse, session), 'u');
u2 = u;

%%
figure, hold on
imshowpair(u1.mimg{2}, imresize(u2.mimg{1}, 1/0.8235))

%%
figure, hold on
imshowpair(u1.mimg{2}, imresize(u2.mimg{4}, 1/0.8235))

%%
figure, hold on
imshowpair(u1.mimg{2}, u1.mimg{4})
