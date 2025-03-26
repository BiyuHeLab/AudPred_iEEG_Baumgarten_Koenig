clear

load 'series_selection_5TrialsPerBeta_jenn.mat'
% load series_selection_H.01_beta.99_1000.mat
% load series_selection_tf_Hz.mat

% define allowable tone frequencies
maxFreq       = 880*2;
nOctaves      = 4;
toneRange     = [maxFreq/2^nOctaves maxFreq]; 
nTonesInRange = 12*nOctaves+1;
logToneFreqs  = linspace(log(toneRange(1)), log(toneRange(2)), nTonesInRange);


% scale series to Hz scale centered on 440 Hz and discretize
for i_b = 1:length(betas)
    for i_s = 1:3
    
        if betas(i_b) == 0
            min_s = min(series_rp{i_b}(i_s,:));
            min_v = log(220);

            series_logHz{i_b}(i_s,:)       = series_rp{i_b}(i_s,:) - min_s + min_v;
%             fGns_logHz{i_b}(i_s,:)         = fGns{i_b}(i_s,:) - min_s + min_v;
            series_logHz_tf_SD{i_b}(i_s,:) = series_rp_tf_SD{i_b}(i_s,:) - min_s + min_v;
%             series_logHz_tf_Hz{i_b}(i_s,:) = series_rp_tf_Hz{i_b}(i_s,:) - min_s + min_v;

            series_pred_logHz{i_b}(i_s,1)  = log(440);
            
            
            
            
            min_s = min(series_rp_sc{i_b}(i_s,:));
            r_sc = range(series_rp_sc{i_b}(i_s,:));
            min_v = log(440) - r_sc/2;
            
            series_logHz_sc{i_b}(i_s,:)       = series_rp_sc{i_b}(i_s,:) - min_s + min_v;
%             fGns_logHz_sc{i_b}(i_s,:)         = fGns_sc{i_b}(i_s,:) - min_s + min_v;
            series_logHz_sc_tf_SD{i_b}(i_s,:) = series_rp_sc_tf_SD{i_b}(i_s,:) - min_s + min_v;
%             series_logHz_sc_tf_Hz{i_b}(i_s,:) = series_rp_sc_tf_Hz{i_b}(i_s,:) - min_s + min_v;
            series_pred_logHz_sc{i_b}(i_s,1)  = log(440);
            
            
            

            
        elseif betas(i_b) == 2
            
            min_s = min(series_rp{i_b}(i_s,:));
            min_v = log(220);

            series_logHz{i_b}(i_s,:)       = series_rp{i_b}(i_s,:) - min_s + min_v;
%             fGns_logHz{i_b}(i_s,:)         = fGns{i_b}(i_s,:);
            series_logHz_tf_SD{i_b}(i_s,:) = series_rp_tf_SD{i_b}(i_s,:) - min_s + min_v;
%             series_logHz_tf_Hz{i_b}(i_s,:) = series_rp_tf_Hz{i_b}(i_s,:) - min_s + min_v;

            series_pred_logHz{i_b}(i_s,1)  = series_logHz{i_b}(i_s,end);
            
            
            
            
            
            
            min_s = min(series_rp_sc{i_b}(i_s,:));
            r_sc = range(series_rp_sc{i_b}(i_s,:));
            min_v = log(440) - r_sc/2;
            
            series_logHz_sc{i_b}(i_s,:)       = series_rp_sc{i_b}(i_s,:) - min_s + min_v;
%             fGns_logHz_sc{i_b}(i_s,:)         = fGns_sc{i_b}(i_s,:);
            series_logHz_sc_tf_SD{i_b}(i_s,:) = series_rp_sc_tf_SD{i_b}(i_s,:) - min_s + min_v;
%             series_logHz_sc_tf_Hz{i_b}(i_s,:) = series_rp_sc_tf_Hz{i_b}(i_s,:) - min_s + min_v;
            
            series_pred_logHz_sc{i_b}(i_s,1)  = series_logHz_sc{i_b}(i_s,end);
            
            
            
        else
            
% % % % %             series(i_s,:) = series(i_s,:) + mean_adj;
% % % % %             series_pred(i_s,:) = series_pred(i_s,:) + mean_adj;
% % % % %             if beta < 1
% % % % %                 fGns(i_s,:) = fGns(i_s,:) + mean_adj;
% % % % %             end

            mean_adj = log(440);

            series_logHz{i_b}(i_s,:)       = series_rp{i_b}(i_s,:) + mean_adj;
%             if beta < 1
%                 fGns_logHz{i_b}(i_s,:)     = fGns{i_b}(i_s,:) + mean_adj;
%             end
%             series_logHz_tf_SD{i_b}(i_s,:) = series_rp_tf_SD{i_b}(i_s,:) + mean_adj;
%             series_logHz_tf_Hz{i_b}(i_s,:) = series_rp_tf_Hz{i_b}(i_s,:) + mean_adj;

            series_pred_logHz{i_b}(i_s,1)  = series_pred_rp{i_b}(i_s,1) + mean_adj;
            



            
%             series_logHz_sc{i_b}(i_s,:)       = series_rp_sc{i_b}(i_s,:) + mean_adj;
%             if beta < 1
%                 fGns_logHz_sc{i_b}(i_s,:)     = fGns_sc{i_b}(i_s,:) + mean_adj;
%             end
%             series_logHz_sc_tf_SD{i_b}(i_s,:) = series_rp_sc_tf_SD{i_b}(i_s,:) + mean_adj;
%             series_logHz_sc_tf_Hz{i_b}(i_s,:) = series_rp_sc_tf_Hz{i_b}(i_s,:) + mean_adj;

%             series_pred_logHz_sc{i_b}(i_s,1)  = series_pred_rp_sc{i_b}(i_s,1) + mean_adj;
            
         
        end
        
        
        discr = 1;
        if discr
            % discretize log frequencies            
            %series_logHz{i_b}(i_s,:)       = discretize(series_logHz{i_b}(i_s,:), logToneFreqs);
%             series_logHz_tf_SD{i_b}(i_s,:) = discretize(series_logHz_tf_SD{i_b}(i_s,:), logToneFreqs);
%             series_logHz_tf_Hz{i_b}(i_s,:) = discretize(series_logHz_tf_Hz{i_b}(i_s,:), logToneFreqs);
            %series_pred_logHz{i_b}(i_s,:)  = discretize(series_pred_logHz{i_b}(i_s,:), logToneFreqs);


%             % discretize log frequencies for _sc series
%             series_logHz_sc{i_b}(i_s,:)       = discretize(series_logHz_sc{i_b}(i_s,:), logToneFreqs);
%             series_logHz_sc_tf_SD{i_b}(i_s,:) = discretize(series_logHz_sc_tf_SD{i_b}(i_s,:), logToneFreqs);
% %             series_logHz_sc_tf_Hz{i_b}(i_s,:) = discretize(series_logHz_sc_tf_Hz{i_b}(i_s,:), logToneFreqs);
%             series_pred_logHz_sc{i_b}(i_s,:)  = discretize(series_pred_logHz_sc{i_b}(i_s,:), logToneFreqs);
        end
      
        
        % re-define tf_Hz using discretized offsets
        xx = [4 8 12];
        nSemitones_offsets = [-xx(end:-1:1) xx];

        ss = series_logHz{i_b}(i_s,end);
        [m ind] = min( abs(ss - logToneFreqs) );
        for i_x = 1:length(nSemitones_offsets)
            series_logHz_tf_Hz{i_b}(i_s,i_x) = logToneFreqs(ind + nSemitones_offsets(i_x));
        end

%         ss = series_logHz_sc{i_b}(i_s,end);
%         [m ind] = min( abs(ss - logToneFreqs) );
%         for i_x = 1:length(nSemitones_offsets)
%             series_logHz_sc_tf_Hz{i_b}(i_s,i_x) = logToneFreqs(ind + nSemitones_offsets(i_x));
%         end
                
        
        
        
        series_Hz{i_b}(i_s,:)        = exp(series_logHz{i_b}(i_s,:));
%             fGns_Hz{i_b}(i_s,:)         = exp(fGns_logHz{i_b}(i_s,:));
%         series_Hz_tf_SD{i_b}(i_s,:)  = exp(series_logHz_tf_SD{i_b}(i_s,:));
        series_Hz_tf_Hz{i_b}(i_s,:)  = exp(series_logHz_tf_Hz{i_b}(i_s,:));
%         series_Hz_tf_Hz2{i_b}(i_s,:) = exp(series_logHz_tf_Hz2{i_b}(i_s,:));
        series_pred_Hz{i_b}(i_s,1)   = exp(series_pred_logHz{i_b}(i_s,1));


%         series_Hz_sc{i_b}(i_s,:)        = exp(series_logHz_sc{i_b}(i_s,:));
% %             fGns_Hz_sc{i_b}(i_s,:)         = exp(fGns_logHz_sc{i_b}(i_s,:));
%         series_Hz_sc_tf_SD{i_b}(i_s,:)  = exp(series_logHz_sc_tf_SD{i_b}(i_s,:));
%         series_Hz_sc_tf_Hz{i_b}(i_s,:)  = exp(series_logHz_sc_tf_Hz{i_b}(i_s,:));
% %         series_Hz_sc_tf_Hz2{i_b}(i_s,:) = exp(series_logHz_sc_tf_Hz2{i_b}(i_s,:));        
%         series_pred_Hz_sc{i_b}(i_s,1)   = exp(series_pred_logHz_sc{i_b}(i_s,1));
        
        
    end
end


% find the unselected series
clear inds
for i_b = 1:length(betas)
    inds{i_b} = [];
    
    for i_s = 1:3
        for i_st = 1:100
            
            if all( series{i_b}(i_st,:) == series_rp{i_b}(i_s,:) )
                inds{i_b} = [inds{i_b} i_st];
%                 break;
            end
            
        end
    end
    
    ind2 = setdiff(1:1000, inds{i_b});
    series_u{i_b} = series{i_b}(ind2',:);
end

% scale the unselected series
for i_b = 1:length(betas)
    for i_s = 1:997
        
        if betas(i_b) == 0 || betas(i_b) == 2
            min_s = min(series_u{i_b}(i_s,:));
            min_v = log(220);
            series_u_logHz{i_b}(i_s,:) = series_u{i_b}(i_s,:) - min_s + min_v;
        
        else
            mean_adj = log(440);
            series_u_logHz{i_b}(i_s,:) = series_u{i_b}(i_s,:) + mean_adj;
                   
        end
        
        %%!!%series_u_logHz{i_b}(i_s,:) = discretize(series_u_logHz{i_b}(i_s,:), logToneFreqs);        
        
        %%!!%series_u_Hz{i_b}(i_s,:) = exp(series_u_logHz{i_b}(i_s,:));
        
    end
end
            

% the subject can only judge sigma_e and the next predicted tone based on
% what they actually hear...
%
% define sigma_e and next predicted tone in a log scale based on actually
% presented stimuli
for i_b = 1:length(betas)
    for i_s = 1:3
        
        s = series_logHz{i_b}(i_s,:); % take log of discretized Hz series
        
        if betas(i_b) < 1
            % for fGn, adjust to zero mean
            adjf = mean(s);
            s_adj = s - adjf;
            fGn = s_adj;
        else
            % no need to adjust for the fGn acquired from the fBm
            fGn = s(2:end) - s(1:end-1);
        end
        
        
        sigma_recov(i_b, i_s) = std(fGn); % observer's estimate of sigma
        s_rec = sigma_recov(i_b, i_s);
        
        %H_recov(i_b, i_s) = wavelet_est(fGn,1); % estimate H
        %h_rec = H_recov(i_b, i_s);
        %if h_rec <= 0, h_rec = 1e-10; end
        
        %sigma_e_recov(i_b, i_s) = sqrt( fgn_sigma2_e( h_rec, s_rec^2, k ) ); % estimate sigma_e
        
        
        %x_pred = fgn_pred( fGn(end:-1:1), h_rec, s_rec^2);  % estimate E[ next tone ]
%         if betas(i_b) < 1
%             % for fGn, add the mean back in to get the "true" predicted
%             % value
%             tf_pred_recov(i_b, i_s) = x_pred + adjf;
%         else
%             % for fBm, increment the last fBm value by the predicted fGn
%             tf_pred_recov(i_b, i_s) = x_pred + s(end);
%         end
%         
        
        
%         % repeat for _sc series
%         s = series_logHz_sc{i_b}(i_s,:); % take log of discretized Hz series
%         
%         if betas(i_b) < 1
%             % for fGn, adjust to zero mean
%             adjf = mean(s);
%             s_adj = s - adjf;
%             fGn = s_adj;
%         else
%             % no need to adjust for the fGn acquired from the fBm
%             fGn = s(2:end) - s(1:end-1);
%         end
%         
%         
%         sigma_recov_sc(i_b, i_s) = std(fGn); % observer's estimate of sigma
%         s_rec = sigma_recov_sc(i_b, i_s);
%         
%         H_recov_sc(i_b, i_s) = wavelet_est(fGn,1); % estimate H
%         h_rec = H_recov_sc(i_b, i_s);
%         if h_rec <= 0, h_rec = 1e-10; end
%         
%         sigma_e_recov_sc(i_b, i_s) = sqrt( fgn_sigma2_e( h_rec, s_rec^2, k ) ); % estimate sigma_e
%         
%         
%         x_pred = fgn_pred( fGn(end:-1:1), h_rec, s_rec^2);  % estimate E[ next tone ]
%         if betas(i_b) < 1
%             % for fGn, add the mean back in to get the "true" predicted
%             % value
%             tf_pred_recov_sc(i_b, i_s) = x_pred + adjf;
%         else
%             % for fBm, increment the last fBm value by the predicted fGn
%             tf_pred_recov_sc(i_b, i_s) = x_pred + s(end);
%         end
          
    end
end

% save('series_tf_recov.mat', ...
%     'sigma_recov', 'H_recov', 'sigma_e_recov', 'tf_pred_recov');

        



% % SD plot
% figure;
% i_p = 1;
% for i_b = 1:length(betas)
%     for i_s = 1:3
%         subplot(5,3,i_p); hold on;
%         title(['t_f by SD, beta = ' num2str(betas(i_b))])
%         
%         plot(1:k, series_Hz{i_b}(i_s,:), 'b-');
%         plot(k+1, series_pred_Hz{i_b}(i_s), 'ro');
%         plot(k+1, series_Hz_tf_SD{i_b}(i_s,:),'r.');
%         
%         plot(1:k, series_Hz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_Hz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_Hz_sc_tf_SD{i_b}(i_s,:),'k.');
%         
%         plot(1:k, 220*ones(1,k), 'r-');
%         plot(1:k, 880*ones(1,k), 'r-');
%         
%         disp(['i_b=' num2str(i_b) ', i_s=' num2str(i_s)])
%         disp(['series end=' num2str(series_Hz{i_b}(i_s,end)) ', sc end=' num2str(series_Hz_sc{i_b}(i_s,end))])
%         disp(' ')
%         
%         i_p = i_p + 1;
%     end
% end


% Hz plot
figure;
i_p = 1;
for i_b = 1:length(betas)
    for i_s = 1:3
        subplot(3,3,i_p); hold on;
        title(['t_f by logHz, beta = ' num2str(betas(i_b))])
        
        plot(1:k, series_Hz{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_Hz{i_b}(i_s), 'ro');
        plot(k+1, series_Hz_tf_Hz{i_b}(i_s,:),'r.');
        
%         plot(1:k, series_Hz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_Hz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_Hz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
        plot(1:k, 220*ones(1,k), 'r-');
        plot(1:k, 880*ones(1,k), 'r-');
        
        i_p = i_p + 1;
    end
end



% logHz plot
figure;
i_p = 1;
for i_b = 1:length(betas)
    for i_s = 1:3
        subplot(3,3,i_p); hold on;
        title(['t_f by logHz, beta = ' num2str(betas(i_b))])
        
        plot(1:k, series_logHz{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
        plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
        
%         plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
        plot(1:k, log(220)*ones(1,k), 'r-');
        plot(1:k, log(880)*ones(1,k), 'r-');
        
        i_p = i_p + 1;
    end
end


% % logHz plot w/ 3 beta levels
% figd;
% i_p = 1;
% betas_title = [0 0 1 1 2];
% for i_b = [1 3 5]
%     for i_s = 1:3
%         subplot(3,3,i_p); hold on;
%         title(['\beta = ' num2str(betas_title(i_b))])
%         
%         plot(1:k, series_logHz{i_b}(i_s,:), 'b-');
%         plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
%         plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
%         
%         plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
%         
% %         plot(1:k, log(220)*ones(1,k), 'r-');
% %         plot(1:k, log(880)*ones(1,k), 'r-');
%         
%         set(gca,'XTick',[]);
%         set(gca,'YTick',[]);
%         
%         if i_p == 4
%             ylabel('Pitch (log Hz)')
%         elseif i_p == 8
%             xlabel('Time')
%         end
% 
%         i_p = i_p + 1;
%     end
% end


% logHz plot for all 5 beta levels, same as above
figd;
i_p = 1;
betas_title = [0.5 0.99 1.5];
for i_b = 1:3
    for i_s = 1:3
        subplot(3,3,i_p); hold on;
        title(['\beta = ' num2str(betas_title(i_b))])
        
        plot(1:k, series_logHz{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
        plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
        
%         plot(1:k, series_logHz_sc{i_b}(i_s,:), 'g-');
%         plot(k+1, series_pred_logHz_sc{i_b}(i_s), 'ko');
%         plot(k+1, series_logHz_sc_tf_Hz{i_b}(i_s,:),'k.');
        
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



% % range plot
% figure;
% i_p = 1;
% for i_b = 1:length(betas)
%     for i_s = 1:3
%         subplot(3,3,i_p); hold on;
%         title(['beta = ' num2str(betas(i_b))])
%         
%         logr    = log2( max(series_Hz{i_b}(i_s,:)) / min(series_Hz{i_b}(i_s,:)) );
%         logr_sc = log2( max(series_Hz_sc{i_b}(i_s,:)) / min(series_Hz_sc{i_b}(i_s,:)) );
%         
%         oct_sc(i_b,i_s) = logr_sc;
%         
%         bar([logr logr_sc])
%         ylabel('# octaves')
%         xlabel({'high \sigma','low \sigma'})
%         
%        
%         i_p = i_p + 1;
%     end
% end
% oct_sc
% mean(oct_sc(:))

close all
save series_selection_scale2Hz
% save series_selection_scale2Hz_k100