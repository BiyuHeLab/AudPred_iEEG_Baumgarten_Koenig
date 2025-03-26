function h = sfa_expt3_plot_stim(stim, filename)

if ~exist('filename','var')
    filename = [];
end

betaIDs = unique(stim.betaID);
predIDs = unique(stim.predID);

h = figd(15,1);

count = 0;
for i_beta = 1:length(betaIDs)
    for i_pred = 1:length(predIDs)

        count = count + 1;
        subplot(3,3,count); hold on;
        set(gca, 'YTickLabel', []);
        set(gca, 'XTickLabel', []);
        
        f = find(stim.betaID == betaIDs(i_beta) ...
               & stim.predID == predIDs(i_pred) );

        plot(1:33, log(stim.series_f{f(1)}(1:end-1)), 'b-');
        plot(34, stim.logf_pred(f(1)), 'rx');
        
        tf = unique(stim.logf_final);
        for j = 1:6
           plot(34, tf(j), 'k.');
        end
        
        xlim([1 length(stim.series_f{f(1)})]);

        if i_pred == 1
            ylabel(['\beta = ' num2str(stim.beta(f(1))) ]);
        end
        
        ttext = {'low', 'med', 'high'};
        if i_beta == 1
            title(['p*_{34} = ' ttext{i_pred}]);
        end
    end
end


if ~isempty(filename)
    save_fig(h, filename, 1, 1, 'png');
end