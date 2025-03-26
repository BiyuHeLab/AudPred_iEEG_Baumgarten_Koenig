H = [.75 .995 .25];

k = 2^7+1;

% logHz plot for all 5 beta levels, same as above
figd;
i_p = 1;
betas_title = {'trend strength = 1', 'trend strength = 2', 'trend strength = 3'};
for i_b = 1:3
    for i_s = 1:5
        subplot(3,5,i_p); hold on;
        title(betas_title{i_b})
        
        [B x] = synthfbmcircul(k, H(i_b));
        
        if i_b == 3
            ss=B;
        else
            ss=x;
        end
        
        plot(1:k, ss, 'k-');
        xlim([1 k])
        ylim([min(ss) max(ss)])
%         plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
%         plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
        
%         plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
%         plot(1:k, log(220)*ones(1,k), 'r-');
%         plot(1:k, log(880)*ones(1,k), 'r-');
        
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
%         if i_p == 7
%             ylabel('Pitch (log Hz)')
%         elseif i_p == 14
%             xlabel('Time')
%         end

        i_p = i_p + 1;
    end
end