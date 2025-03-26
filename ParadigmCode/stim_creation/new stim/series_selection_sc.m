clear 

% load series_selection.mat
load series_selection_H.01.mat



%% select series for the scaled stim

for i_b = 1:length(betas)

%     % original code for selecting min/max/median series prediction    
%     [m minInd] = min(series_pred_r{i_b});
%     [m maxInd] = max(series_pred_r{i_b});
%     [m medInd] = min( abs(series_pred_r{i_b} - median(series_pred_r{i_b})) );
%
%     inds = [minInd maxInd medInd];
% 
%     series_rp{i_b}      = series_r{i_b}(inds',:);
%     series_pred_rp{i_b} = series_pred_r{i_b}(inds',:);
%     pred_err_rp{i_b}    = pred_err_r{i_b}(inds',:);    
  

    % new code for selecting "penultimate" min/max/median series prediction
    [xx_sort ind_sort] = sort(series_pred_r{i_b});
    minInd1 = ind_sort(1);
    maxInd1 = ind_sort(end);
        
    xx = abs(series_pred_r{i_b} - median(series_pred_r{i_b}));
    [xx_sort ind_sort] = sort(xx);
    medInd1 = ind_sort(1);

    
%     medInd1 = ind_sort( round(length(ind_sort)/2) );


    
    inds1 = [minInd1 maxInd1 medInd1];

    series_rp{i_b}      = series_r{i_b}(inds1',:);
    series_pred_rp{i_b} = series_pred_r{i_b}(inds1',:);
    pred_err_rp{i_b}    = pred_err_r{i_b}(inds1',:);    
    
end





%% original series_selection_sc code
sigma_w = .7;
for i_b = 1:length(betas)

%     sigma           = sqrt(sigma2s(i_b));
%     sigma_sc        = sigma * sigma_w;
%     sigma2s_sc(i_b) = sigma_sc^2;

    for i_s = 1:3

        beta  = betas(i_b);
        if beta > 1, H = (beta-1)/2; else H = (beta+1)/2; end
        
        s = series_rp{i_b}(i_s,:);

        if beta < 1
            sigma_rp(i_b, i_s) = std(s);
        else
            s_fGn = s(2:end) - s(1:end-1);
            sigma_rp(i_b, i_s) = std(s_fGn);
        end
        
        sigma_rp_sc(i_b, i_s) = sigma_rp(i_b, i_s) * sigma_w;
        

%         [s_adj fGn_adj fGn] = fGn_change_sigma(s, beta, sigma_sc);
        [s_adj fGn_adj fGn] = fGn_change_sigma(s, beta, sigma_rp_sc(i_b, i_s));

%         s_adj = s_adj - (s_adj(end) - v_f); % align final sample of s_adj to v_f
        
        adjf = -s_adj(end) + s(end);
%         s_adj = s_adj - s_adj(end) + s(end); % align final sample of s_adj to that of s
        s_adj = s_adj + adjf; % align final sample of s_adj to that of s
        

        series_rp_sc{i_b}(i_s,:) = s_adj;

%         x_pred = fgn_pred(fGn_adj(end:-1:1), H, sigma_sc^2);
        x_pred = fgn_pred(fGn_adj(end:-1:1), H, sigma_rp_sc(i_b,i_s)^2);

        if beta > 1
            series_pred_rp_sc{i_b}(i_s,1) = x_pred + s_adj(end);
        else
            series_pred_rp_sc{i_b}(i_s,1) = x_pred + adjf;
        end
        
        fGn_std(i_b, i_s) = std(fGn);
        fGn_sc_std(i_b, i_s) = std(fGn_adj);
        fGn_mean(i_b, i_s) = mean(fGn);
        fGn_sc_mean(i_b, i_s) = mean(fGn_adj);
        
    end
end






%% select new series for the non-scaled stim

for i_b = 1:length(betas)

%     % original code for selecting min/max/median series prediction    
%     [m minInd] = min(series_pred_r{i_b});
%     [m maxInd] = max(series_pred_r{i_b});
%     [m medInd] = min( abs(series_pred_r{i_b} - median(series_pred_r{i_b})) );
%
%     inds = [minInd maxInd medInd];
% 
%     series_rp{i_b}      = series_r{i_b}(inds',:);
%     series_pred_rp{i_b} = series_pred_r{i_b}(inds',:);
%     pred_err_rp{i_b}    = pred_err_r{i_b}(inds',:);    
  

    % new code for selecting "penultimate" min/max/median series prediction
    [xx_sort ind_sort] = sort(series_pred_r{i_b});
    minInd2 = ind_sort(2);
    maxInd2 = ind_sort(end-1);
        
    xx = abs(series_pred_r{i_b} - median(series_pred_r{i_b}));
    [xx_sort ind_sort] = sort(xx);
    medInd2 = ind_sort(2);
    
%     medInd2 = ind_sort( 1 + round(length(ind_sort)/2) );

    
    inds2 = [minInd2 maxInd2 medInd2];

    series_rp{i_b}      = series_r{i_b}(inds2',:);
    series_pred_rp{i_b} = series_pred_r{i_b}(inds2',:);
    pred_err_rp{i_b}    = pred_err_r{i_b}(inds2',:);    
    
    
    for i_s = 1:3
        
        % re-calculate sigma_rp for newly selected series_rp
        beta  = betas(i_b);
        if beta > 1, H = (beta-1)/2; else H = (beta+1)/2; end
        
        s = series_rp{i_b}(i_s,:);

        if beta < 1
            sigma_rp(i_b, i_s) = std(s);
        else
            s_fGn = s(2:end) - s(1:end-1);
            sigma_rp(i_b, i_s) = std(s_fGn);
        end
        

        % re-calculate fGn_std and fGn_mean for newly selected series_rp
        [s_adj fGn_adj fGn] = fGn_change_sigma(s, beta, sigma_rp_sc(i_b, i_s));
        
        fGn_std(i_b, i_s) = std(fGn);
        fGn_mean(i_b, i_s) = mean(fGn);
        
        
    end
end




    
    
    
%% plot

figure;
i_p = 1;
for i_b = 1:length(betas)
    for i_s = 1:3
        subplot(5,3,i_p); hold on;
        title(['sc, beta = ' num2str(betas(i_b)) ]);
        
        plot(1:k, series_rp{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_rp{i_b}(i_s), 'ro');
        
        plot(1:k, series_rp_sc{i_b}(i_s,:), 'g-');
        plot(k+1, series_pred_rp_sc{i_b}(i_s), 'ko');
        
        i_p = i_p + 1;
    end
end

% save series_selection_sc.mat