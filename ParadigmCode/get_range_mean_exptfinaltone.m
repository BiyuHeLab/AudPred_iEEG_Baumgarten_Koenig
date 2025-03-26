%% get mean and range for expected final tone bins (sfa_expt3)
load('series_selection_scale2Hz_expt3_v1')
for i_bin = 1:3
    for i_beta = 1:3
        exptone(i_beta,i_bin) = series_pred_Hz{i_beta}(i_bin);
    end
end

mintone = min(exptone);
maxtone = max(exptone);
meantone = mean(exptone);