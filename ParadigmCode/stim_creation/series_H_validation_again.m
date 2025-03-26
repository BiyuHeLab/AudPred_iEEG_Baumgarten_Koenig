%% Validation of auditory stimulus H values...


%% 
% starup

clear

% cd('C:\Users\maniscalcobs\Desktop\sfa\stim_creation');
% % load('C:\Users\maniscalcobs\Desktop\sfa\stim_creation\series_selection_scale2Hz.mat');
% load('C:\Users\maniscalcobs\Desktop\sfa_expt2\series_selection_scale2Hz.mat');

load('C:\Users\maniscalcobs\Desktop\lab\experiments\sfa\sfa_expt2_noconf\series_selection_scale2Hz.mat')

%%
% 1. Verify that series_rp (the originally created sequence) is almost exactly 
% the same as series_logHz (the same sequence, scaled and translated into log Hz 
% values used in experiment) after scaling (allowing for sound roundoff differences
% due to discretization etc)
%
% Results: they look very similar

% % % % plot the original series_rp and modified series_logHz after scaling
% % % for i_b = 1:5
% % %     for i_s = 1:3
% % %         figure; hold on; 
% % %         plot(Scale(series_rp{i_b}(i_s,:)),'b-'); 
% % %         plot(Scale(series_logHz{i_b}(i_s,:)),'r-');
% % %     end
% % % end
% % % 
% % % % do the same for the smaller SD stimuli
% % % for i_b = 1:5
% % %     for i_s = 1:3
% % %         figure; hold on; 
% % %         plot(Scale(series_rp_sc{i_b}(i_s,:)),'b-'); 
% % %         plot(Scale(series_logHz_sc{i_b}(i_s,:)),'r-');
% % %     end
% % % end


%% 
% 2. Retrieve the H values for each waveform and verify they are close to 
% what they should be
%
% Results
%
% mean(H_est,2)
% ans =
%     0.5045
%     0.7534
%     0.0001
%     0.2479
%     0.4875
%     
% mean(H_est_sc,2)
% ans =
%     0.5045
%     0.7534
%    -0.0274
%     0.2529
%     0.5087
%     
% mean(beta_est,2)
% ans =
%     0.0090
%     0.5067
%     1.0003
%     1.4957
%     1.9750
%     
% mean(beta_est_sc,2)
% ans =
%     0.0090
%     0.5067
%     0.9452
%     1.5058
%     2.0175

Hs(1:2) = (betas(1:2)+1)/2;
Hs(3:5) = (betas(3:5)-1)/2;

for i_b = 1:length(betas)
    for i_s = 1:3
        
        if betas(i_b) < 1
            x = series_rp{i_b}(i_s,:);
%             x = series_logHz{i_b}(i_s,:);

            H_est(i_b, i_s) = wavelet_est(x,1);
            H_err(i_b, i_s) = H_est(i_b, i_s) - Hs(i_b);
            
            beta_est(i_b, i_s) = 2*H_est(i_b,i_s) - 1;
            beta_err(i_b, i_s) = beta_est(i_b, i_s) - betas(i_b);
            
            x = series_rp_sc{i_b}(i_s,:);
%             x = series_logHz_sc{i_b}(i_s,:);

            H_est_sc(i_b, i_s) = wavelet_est(x,1);
            H_err_sc(i_b, i_s) = H_est_sc(i_b, i_s) - Hs(i_b);

            beta_est_sc(i_b, i_s) = 2*H_est_sc(i_b,i_s) - 1;
            beta_err_sc(i_b, i_s) = beta_est_sc(i_b, i_s) - betas(i_b);
            
        
        else
            B = series_rp{i_b}(i_s,:);            
%             B = series_logHz{i_b}(i_s,:);

            x = B(2:end) - B(1:end-1);
            H_est(i_b, i_s) = wavelet_est(x,1);
            H_err(i_b, i_s) = H_est(i_b, i_s) - Hs(i_b);
            
            beta_est(i_b, i_s) = 2*H_est(i_b,i_s) + 1;
            beta_err(i_b, i_s) = beta_est(i_b, i_s) - betas(i_b);

            B = series_rp_sc{i_b}(i_s,:);                        
%             B = series_logHz_sc{i_b}(i_s,:);

            x = B(2:end) - B(1:end-1);
            H_est_sc(i_b, i_s) = wavelet_est(x,1);
            H_err_sc(i_b, i_s) = H_est_sc(i_b, i_s) - Hs(i_b);
            
            beta_est_sc(i_b, i_s) = 2*H_est_sc(i_b,i_s) + 1;
            beta_err_sc(i_b, i_s) = beta_est_sc(i_b, i_s) - betas(i_b);

        end
    end
end

%%
% 3. Visually inspect each waveform and its associated beta
%
% Results: looks ok...

figure;
i_p = 1;
for i_b = 1:length(betas)
    for i_s = 1:3
        subplot(5,3,i_p); hold on;
        title(['t_f by logHz, beta = ' num2str(betas(i_b))])
        
        plot(1:k, series_logHz{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
        plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
        
        plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g-');
        plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'ko');
        plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
        plot(1:k, log(220)*ones(1,k), 'r-');
        plot(1:k, log(880)*ones(1,k), 'r-');

        title(['recovered betas = ' num2str(beta_est(i_b, i_s)) ', ' num2str(beta_est_sc(i_b,i_s))])
        
        i_p = i_p + 1;
    end
end
