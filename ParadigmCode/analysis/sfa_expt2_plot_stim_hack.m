function h = sfa_expt2_plot_stim_hack(stim, filename)

if ~exist('filename','var')
    filename = [];
end

betaIDs = unique(stim.betaID);
sigma_e_IDs = unique(stim.sigma_e_ID);
predIDs = unique(stim.predID);

tf = unique(stim.logf_final);

h = figd(25);
styles = {'g.-', 'b.-'};

for i_beta = 5
    for i_pred = 1

        hold on;
        set(gca, 'YTickLabel', []);
        set(gca, 'XTickLabel', []);
        
        for i_sig = 1:2;
        
            if i_sig == 1
            
                f = find(stim.betaID == 5 ...
                       & stim.predID == 1 ...
                       & stim.sigma_e_ID == sigma_e_IDs(2));
            else
                f = find(stim.betaID == 3 ...
                       & stim.predID == -1 ...
                       & stim.sigma_e_ID == sigma_e_IDs(2));
            end
            
            plot(log(stim.series_f{f(1)}(1:end-1)), styles{i_sig});
        end
    end
end


for i_beta = 5
    for i_pred = 1

        hold on;
        set(gca, 'YTickLabel', []);
        set(gca, 'XTickLabel', []);
        
        for i_sig = 1:2;
        
            if i_sig == 1
            
                f = find(stim.betaID == 5 ...
                       & stim.predID == 1 ...
                       & stim.sigma_e_ID == sigma_e_IDs(2));
            else
                f = find(stim.betaID == 3 ...
                       & stim.predID == -1 ...
                       & stim.sigma_e_ID == sigma_e_IDs(2));
            end
            
            plot(34*ones(size(tf)), tf, 'ko');
           
            x = stim.logf_pred(f);
            if i_sig == 1
                plot(34, x(1), 'g*');
            else
                plot(34, x(1), 'b*');
            end
            
            xlim([1 length(stim.series_f{f(1)})]);
        end
    end
end

legend('\beta = 2', '\beta = 1.01', 'Location', 'NorthWest');

ylabel('tone pitch (log Hz)')
xlabel('time')

if ~isempty(filename)
    save_fig(h, filename, 1, 1, 'png');
end