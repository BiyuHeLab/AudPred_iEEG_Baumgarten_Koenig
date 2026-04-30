function plot_subset_kratio(subset_mask, subset_name, ...
    data_selElec_allSubs, sumSignELec, num_elecs_above1, num_elecs_below1, ...
    mean_KprimeRatio, dv_absmax, dv_absmin, ind_infentries, ...
    pval_threshold, FDR_label, ...
    windows, num_win_common, fsample, ...
    subs, input_DataLabel, ...
    save_poststepFigs, path_fig, param)

% Extract subset data
coords     = data_selElec_allSubs.elecpos(subset_mask, :);
vals       = data_selElec_allSubs.RatioTD_ExpKprime_selElecs_avgWin(subset_mask);
sign_index = subset_mask(subset_mask);

% Compute color limits
dv_abs = max([abs(dv_absmax) abs(dv_absmin)]);
clims  = [0 2];
cmap   = 'jet';

% Channel size - fixed for all electrodes in subset
chanSize = 2.5 * ones(size(vals));

%% Plotting
fig = figure('visible','on');
set(gcf,'units','normalized','outerposition',[0 0 1 1]) % full screen

% Project on left hemisphere
coords_projectedL = coords;
coords_projectedL(:,1) = abs(coords(:,1)) * -1;

NASTD_ECoG_Plot_SubplotSignElecsSurf_SubplotLH_avgTW(...
    coords_projectedL, vals, sign_index, ...
    chanSize, clims, cmap, [1, 1], 1, []);

fig.Children(2).Ticks = [0 0.5 1 1.5 2];
fig.Children(2).FontSize = 16;

% Scatter plot
subplot(1,1,1);
hold on;
scatter(...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}(subset_mask), ...
    data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}(subset_mask), ...
    100, vals, '.')

% Mean ratio marker
scatter(...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{1}(subset_mask)), ...
    nanmean(data_selElec_allSubs.kprime_selelecs_perTD_avgWin{2}(subset_mask)), ...
    200, nanmean(vals), 'd', 'filled', ...
    'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth', 2)

xlim([0 15]); ylim([0 15]);
caxis([0 2]);
c = colorbar; c.Label.String = 'Kprime Ratio across TD';
c.Ticks = [0 0.5 1 1.5 2];
c.FontSize = 16;
colormap('jet');
xlabel('kPrime 200ms TD','FontSize',16)
ylabel('kPrime 400ms TD','FontSize',16)
axis('square')
grid on

% Title
title_str = {['Subset: ' subset_name ' (n = ' num2str(length(subs)) ') - ' input_DataLabel], ...
    ['p < ' num2str(pval_threshold) ', ' FDR_label], ...
    ['Mean ratio: ' num2str(round(mean_KprimeRatio,2))], ...
    ['max/min: ' num2str(round(dv_absmax,2)) '/' num2str(round(dv_absmin,2))], ...
    [num2str(num_elecs_above1) ' >1, ' num2str(num_elecs_below1) ' <1'], ...
    [num2str(length(ind_infentries)) ' removed (0/Inf)']};
sgtitle(title_str, 'Interpreter','none')

% Save
if save_poststepFigs
    filename = ['3DSurf_KprimeRatio_' subset_name '_Allsubn' num2str(length(subs)) ...
        '_' input_DataLabel '_p' num2str(pval_threshold) FDR_label ...
        '_' param.ElecSelect 'elec_AvgWin.png'];
    figfile = [path_fig filename];
    saveas(gcf, figfile, 'png');
    close all;
end
end