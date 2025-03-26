clear

load 'series_selection_5TrialsPerBeta_jenn.mat'

% define allowable tone frequencies
maxFreq       = 880*2;
nOctaves      = 4;
toneRange     = [maxFreq/2^nOctaves maxFreq]; 
nTonesInRange = 12*nOctaves+1; 
logToneFreqs  = linspace(log(toneRange(1)), log(toneRange(2)), nTonesInRange);

betas = [0 0.99 1.5];
k = 33;

% scale series to Hz scale centered on 440 Hz and discretize
for i_b = 1:length(betas)
    for i_s = 1:size(series_rp{i_b},1)

        mean_adj = log(440);
        
        series_logHz{i_b}(i_s,:)       = series_rp{i_b}(i_s,:) + mean_adj;

        series_pred_logHz{i_b}(i_s,1)  = series_pred_rp{i_b}(i_s,1) + mean_adj;
        
        
        discr = 1;
        if discr
            % discretize log frequencies            
            series_logHz{i_b}(i_s,:)       = discretize(series_logHz{i_b}(i_s,:), logToneFreqs);
            series_pred_logHz{i_b}(i_s,:)  = discretize(series_pred_logHz{i_b}(i_s,:), logToneFreqs);

        end
      
        
        % re-define tf_Hz using discretized offsets
        xx = [4 8 12];
        nSemitones_offsets = [-xx(end:-1:1) xx];

        ss = series_logHz{i_b}(i_s,end);
        [m ind] = min( abs(ss - logToneFreqs) );
        for i_x = 1:length(nSemitones_offsets)
            series_logHz_tf_Hz{i_b}(i_s,i_x) = logToneFreqs(ind + nSemitones_offsets(i_x));
        end

        
        series_Hz{i_b}(i_s,:)        = exp(series_logHz{i_b}(i_s,:));

        series_Hz_tf_Hz{i_b}(i_s,:)  = exp(series_logHz_tf_Hz{i_b}(i_s,:));

        series_pred_Hz{i_b}(i_s,1)   = exp(series_pred_logHz{i_b}(i_s,1));

    end
end

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


% logHz plot for all 5 beta levels, same as above
figure;
i_p = 1;
betas_title = [0.5 0.99 1.5];
for i_b = 1:3
    for i_s = 1:3
        subplot(3,3,i_p); hold on;
        title(['\beta = ' num2str(betas_title(i_b))])
        
        plot(1:k, series_logHz{i_b}(i_s,:), 'b-');
        plot(k+1, series_pred_logHz{i_b}(i_s), 'ro');
        plot(k+1, series_logHz_tf_Hz{i_b}(i_s,:),'r.');
        
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


%close all
save series_selection_scale2Hz_jenn
