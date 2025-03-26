clear
load('C:\Users\maniscalcobs\Desktop\sfa\stim_creation\series_selection_scale2Hz.mat')

Hs(1:2) = (betas(1:2)+1)/2;
Hs(3:5) = (betas(3:5)-1)/2;
H_est_sd = .29;

for i_b = 1:length(betas)
    for i_s = 1:3
        
        if betas(i_b) < 1
            x = series_logHz{i_b}(i_s,:);
            H_est(i_b, i_s) = wavelet_est(x,1);
            H_err(i_b, i_s) = H_est(i_b, i_s) - Hs(i_b);
            
            beta_est(i_b, i_s) = 2*H_est(i_b,i_s) - 1;
            beta_err(i_b, i_s) = beta_est(i_b, i_s) - betas(i_b);
            
            x = series_logHz_sc{i_b}(i_s,:);
            H_est_sc(i_b, i_s) = wavelet_est(x,1);
            H_err_sc(i_b, i_s) = H_est_sc(i_b, i_s) - Hs(i_b);

            beta_est_sc(i_b, i_s) = 2*H_est_sc(i_b,i_s) - 1;
            beta_err_sc(i_b, i_s) = beta_est_sc(i_b, i_s) - betas(i_b);
            
        
        else
            B = series_logHz{i_b}(i_s,:);
            x = B(2:end) - B(1:end-1);
            H_est(i_b, i_s) = wavelet_est(x,1);
            H_err(i_b, i_s) = H_est(i_b, i_s) - Hs(i_b);
            
            beta_est(i_b, i_s) = 2*H_est(i_b,i_s) + 1;
            beta_err(i_b, i_s) = beta_est(i_b, i_s) - betas(i_b);

            B = series_logHz_sc{i_b}(i_s,:);
            x = B(2:end) - B(1:end-1);
            H_est_sc(i_b, i_s) = wavelet_est(x,1);
            H_err_sc(i_b, i_s) = H_est_sc(i_b, i_s) - Hs(i_b);
            
            beta_est_sc(i_b, i_s) = 2*H_est_sc(i_b,i_s) + 1;
            beta_err_sc(i_b, i_s) = beta_est_sc(i_b, i_s) - betas(i_b);

        end
    end
end




%% compare to pilot behvaioral data

subs = {'brian' 'biyu'};
data_dir = 'C:\Users\maniscalcobs\Desktop\sfa\expt2\data\';
prefix = 'sfa_expt2_pilot ';
suffix = '.mat';

for i_sub = 2:2 %length(subs)
    load([data_dir prefix subs{i_sub} suffix]);

    stim.i_s(stim.predID==-1) = 1;
    stim.i_s(stim.predID==+1) = 2;
    stim.i_s(stim.predID== 0) = 3;
    
    
    %%%% filter
    filter = data.resp_beta >= 0 & data.resp_prob > 0;

    clear dataf
    clear stimf
    
    dfn = fieldnames(data);
    for j = 1:length(dfn)
        if eval(['length(data.' dfn{j} ') == 360'])
            eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
        end
    end

    sfn = fieldnames(stim);
    for j = 1:length(sfn)
        if eval(['length(stim.' sfn{j} ') == 360'])
            eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
        end
    end    
    %%%%% end filter


    for i_b = 1:length(betas)
        for i_s = 1:3

            f = stimf.sigma_e_ID == 1;
            g = f & stimf.beta==betas(i_b) & stimf.i_s==i_s;
            beta_est_b(i_b, i_s) = mean(dataf.resp_beta(g));
            beta_est_b_sem(i_b, i_s) = std(dataf.resp_beta(g)) / sqrt(length(dataf.resp_beta(g)));
            beta_err_b(i_b, i_s) = beta_est_b(i_b, i_s) - betas(i_b);
            
            
            f = stimf.sigma_e_ID == 0;
            g = f & stimf.beta==betas(i_b) & stimf.i_s==i_s;
            beta_est_b_sc(i_b, i_s) = mean(dataf.resp_beta(g));
            beta_est_b_sem_sc(i_b, i_s) = std(dataf.resp_beta(g)) / sqrt(length(dataf.resp_beta(g)));
            beta_err_b_sc(i_b, i_s) = beta_est_b_sc(i_b, i_s) - betas(i_b);
           
        end
    end
end

h1=figd;
for i_b = 1:length(betas)
    subplot(1,5,i_b); hold on;
    
    [xx ind] = sort(beta_est(i_b,:));
    errorbar(beta_est(i_b,ind), beta_est_b(i_b,ind), beta_est_b_sem(i_b,ind), 'bo-');
    
    [xx ind] = sort(beta_est_sc(i_b,:));
    errorbar(beta_est_sc(i_b,ind), beta_est_b_sc(i_b,ind), beta_est_b_sem_sc(i_b,ind), 'ro-');

    plot([min([beta_est(i_b,ind) beta_est_sc(i_b,ind)]), max([beta_est(i_b,ind) beta_est_sc(i_b,ind)])], ...
         [min([beta_est(i_b,ind) beta_est_sc(i_b,ind)]), max([beta_est(i_b,ind) beta_est_sc(i_b,ind)])], 'k--');
    
    if i_b == 3
        xlabel('wavelet beta')
    end
    
    title(['true beta = ' num2str(betas(i_b))])
    
    if i_b == 1
        ylabel('behavioral beta')

        legend('high \sigma_\epsilon', 'low \sigma_\epsilon');
    end
end


h2=figd; hold on;
col = linspace(0,1,5);

% beta_err(1,:) = []; beta_err(end,:) = [];
% beta_err_sc(1,:) = []; beta_err_sc(end,:) = [];
% beta_err_b(1,:) = []; beta_err_b(end,:) = [];
% beta_err_b_sc(1,:) = []; beta_err_b_sc(end,:) = [];

wb_err = [beta_err(:); beta_err_sc(:)];
bb_err = [beta_err_b(:); beta_err_b_sc(:)];
[r p]  = corr(wb_err, bb_err, 'type', 'Spearman')
for i_b = 1:length(betas)
    errorbar(beta_err(i_b,ind), beta_err_b(i_b,ind), beta_est_b_sem(i_b,ind), 'bo', 'Color',[col(i_b) 0 0]);
    errorbar(beta_err_sc(i_b,ind), beta_err_b_sc(i_b,ind), beta_est_b_sem_sc(i_b,ind), 'bo', 'Color',[col(i_b) 0 0]);
end

xlabel('wavelet beta estimation error')
ylabel('behavioral beta estimation error')
title(['Spearman''s rho = ' num2str(r,2) ', p = ' num2str(p,2)]);