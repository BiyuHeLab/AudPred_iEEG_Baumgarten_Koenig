load series_selection_sc.mat

xx = [.15 .3 .5];
logHz_offsets = [-xx(end:-1:1) xx];

% define the list of final tones for each series, based on constant log(Hz) intervals
for i_b = 1:length(betas)
    for i_s = 1:3
        series_rp_tf_Hz{i_b}(i_s,:)    = series_rp{i_b}(i_s,end) + logHz_offsets;
        series_rp_sc_tf_Hz{i_b}(i_s,:) = series_rp_sc{i_b}(i_s,end) + logHz_offsets;
    end
end


figure;
i_p = 1;
for i_b = 1:length(betas)
    for i_s = 1:3
        subplot(5,3,i_p); hold on;
        title(['beta = ' num2str(betas(i_b))])
        
        plot(1:k, series_rp{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_rp{i_b}(i_s), 'ro');
        plot(k+1, series_rp_tf_Hz{i_b}(i_s,:),'r.');
        
        plot(1:k, series_rp_sc{i_b}(i_s,:), 'g-');
        plot(k+1, series_pred_rp_sc{i_b}(i_s), 'ko');
        plot(k+1, series_rp_sc_tf_Hz{i_b}(i_s,:),'k.');
        
        i_p = i_p + 1;
    end
end

save series_selection_tf_Hz.mat