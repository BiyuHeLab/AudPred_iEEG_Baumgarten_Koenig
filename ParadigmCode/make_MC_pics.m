%% Make multiple choice pictures

load series_selection_scale2Hz_jenn
path_fig = 'C:\Users\localbiyu\Desktop\Brian\sfa\sfa_expt4\useq_pics\';

sers = series_logHz{3};
set(0, 'DefaultFigureColor', 'remove')

for i_s = 1:size(sers,1)
    h = figd(14,8,3);

    plot([1 33],[log(440) log(440)], 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
    hold on
    plot(sers(i_s,:), '-wo');
    
    hAxes = gca;
    hAxes.XRuler.Axle.LineStyle = 'none';
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    box off
    set(gca,'Color', 'k')
    set(gcf, 'InvertHardCopy', 'off');
    
    saveas(h, [path_fig 'uSeq_' num2str(i_s) '.png']);
end

    