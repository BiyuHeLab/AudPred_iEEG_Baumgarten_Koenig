load series_selection_sc.mat

sd_offsets = [-3 -1.5 -.5 .5 1.5 3];

% calculate sigma2_epsilon
for i_b = 1:length(betas)
    
    beta      = betas(i_b);
    if beta > 1, H = (beta-1)/2; else H = (beta+1)/2; end

    sigma2    = sigma2s(i_b);
%     sigma2_sc = sigma2s_sc(i_b);
    
    sigma_e(i_b)    = sqrt( fgn_sigma2_e(H, sigma2, k) );
%     sigma_e_sc(i_b) = sqrt( fgn_sigma2_e(H, sigma2_sc, k) );

end

% define the list of final tones for each series, based on sigma_e
for i_b = 1:length(betas)
    for i_s = 1:3
        
        sigma_e_sc(i_b,i_s) = sqrt( fgn_sigma2_e(H, sigma_rp_sc(i_b, i_s)^2, k) );
        
        series_rp_tf_SD{i_b}(i_s,:)    = series_pred_rp{i_b}(i_s) + sigma_e(i_b) * sd_offsets;
        series_rp_sc_tf_SD{i_b}(i_s,:) = series_pred_rp_sc{i_b}(i_s) + sigma_e_sc(i_b, i_s) * sd_offsets;
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
        plot(k+1, series_rp_tf_SD{i_b}(i_s,:),'r.');
        
        plot(1:k, series_rp_sc{i_b}(i_s,:), 'g-');
        plot(k+1, series_pred_rp_sc{i_b}(i_s), 'ko');
        plot(k+1, series_rp_sc_tf_SD{i_b}(i_s,:),'k.');
        
        i_p = i_p + 1;
    end
end

save series_selection_tf_SD.mat