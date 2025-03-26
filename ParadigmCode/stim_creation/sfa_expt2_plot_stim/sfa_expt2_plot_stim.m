clear

load series_selection_scale2Hz

% logHz plot for all 5 beta levels
figd;
i_p = 1;
betas_title = [0 0.5 1 1.5 2];
for i_b = 1:5
    for i_s = 1:3
        subplot(5,3,i_p); hold on;
        title(['\beta = ' num2str(betas_title(i_b))])
        
        % high sigma_e stim
        plot(1:k, series_logHz{i_b}(i_s,:), 'b.-', 'MarkerSize', 15);
        plot(k+1, series_pred_logHz{i_b}(i_s), 'bo');
        plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'k.');
        
        % low sigma_e stim
        plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g.-', 'MarkerSize', 15);
        plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'go');
        plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
%         plot(1:k, log(220)*ones(1,k), 'r-');
%         plot(1:k, log(880)*ones(1,k), 'r-');
        
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
        if i_p == 7
            ylabel('Pitch (log Hz)')
        elseif i_p == 14
            xlabel('Time')
        end

        i_p = i_p + 1;
    end
end