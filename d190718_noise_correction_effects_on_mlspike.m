basedir = 'D:\TPM\JK\suite2p\';
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

% for mi = 1 : length(mice)
for mi = 2
    mouse = mice(mi);    
%     for si = 1 : length(sessions{mi})
    for si = 1
        session = sessions{mi}(si);
%         for pi = 1 : 8
        for pi = 4
            newfn = sprintf('F_%03d_%03d_plane%d_proc_final_spikes_noiseCorrected.mat',mouse, session, pi);
            oldfn = sprintf('F_%03d_%03d_plane%d_proc_final_spikes.mat',mouse, session, pi);
            load(sprintf('%s%03d\\%s',basedir,mouse,newfn), 'spk')
            newspk = spk.n;
            load(sprintf('%s%03d\\%s',basedir,mouse,oldfn), 'spk')
            oldspk = spk.n;
            tempCorrVal = zeros(size(newspk,1),1);
            for ci = 1 : length(tempCorrVal)
                tempCorrVal(ci) = corr(oldspk(ci,:)', newspk(ci,:)');
            end
        end
    end
end
