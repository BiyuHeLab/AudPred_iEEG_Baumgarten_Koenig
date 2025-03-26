clear

load series_selection_scale2Hz

nreps = 10;

for i_rep = 1:nreps

%% simulate performance for each of 4 subject/strategy types for probability rating
for i_b = 1:length(betas)
    for i_s = 1:3
        
        % get the relevant stimulus variables for this condition
        t_pen   = series_logHz{i_b}(i_s,end);
        t_pred  = series_pred_logHz{i_b}(i_s,end);
        t_SD    = series_logHz_tf_SD{i_b}(i_s,:); % defined by SDs from t_pred
        t_logHz = series_logHz_tf_Hz2{i_b}(i_s,:); % defined by Hz from t_f
        
        sig_eps = sigma_e(i_b);
        
        t_pen_sc   = series_logHz_sc{i_b}(i_s,end);
        t_pred_sc  = series_pred_logHz_sc{i_b}(i_s,end);
        t_SD_sc    = series_logHz_sc_tf_SD{i_b}(i_s,:); % defined by SDs from t_pred
        t_logHz_sc = series_logHz_sc_tf_Hz2{i_b}(i_s,:); % defined by Hz from t_f
        
        sig_eps_sc = sigma_e_sc(i_b);
        
        
        
        % define method of setting final tone
        t_f_bySD = 1;
        if t_f_bySD
            t_f    = t_SD;
            t_f_sc = t_SD_sc;
        else
            t_f    = t_logHz;
            t_f_sc = t_logHz_sc;
        end
        
        
  
        % define predictor variables     
        % in the actual experiment, we'll repeat each final tone
        % presentation twice, which is why we multiply each vector by [1 1]
        
        % - pitch diff from penultimate tone
        xx = abs( t_f - t_pen );
        del_logHz_pen{i_b}(i_s,:)    = [xx xx];
        
        xx = abs( t_f - t_pen_sc );
        del_logHz_pen_sc{i_b}(i_s,:) = [xx xx];

        
        % - pitch diff from predicted final tone
        xx = abs( t_f - t_pred );
        del_logHz_pred{i_b}(i_s,:)    = [xx xx];
        
        xx = abs( t_f - t_pred_sc );
        del_logHz_pred_sc{i_b}(i_s,:) = [xx xx];

        
        % - sigma_e diff from penultimate tone
        xx = abs( (t_f - t_pen) / sig_eps );
        del_SD_pen{i_b}(i_s,:)    = [xx xx];
        
        xx = abs( (t_f_sc - t_pen_sc) / sig_eps_sc );
        del_SD_pen_sc{i_b}(i_s,:) = [xx xx];        
        
        
        % - sigma_e diff from predicted final tone
        xx = abs( (t_f - t_pred) / sig_eps );
        del_SD_pred{i_b}(i_s,:)    = [xx xx];
        
        xx = abs( (t_f_sc - t_pred_sc) / sig_eps_sc );
        del_SD_pred_sc{i_b}(i_s,:) = [xx xx];

        
        
        % set x-axis for plot        
        t_f_xaxis_bySD = t_f_bySD;
        if t_f_xaxis_bySD
            del_x{i_b}(i_s,:)    = del_SD_pred{i_b}(i_s,:);
            del_x_sc{i_b}(i_s,:) = del_SD_pred_sc{i_b}(i_s,:);
        else
            del_x{i_b}(i_s,:)    = del_logHz_pen{i_b}(i_s,:);
            del_x_sc{i_b}(i_s,:) = del_logHz_pen_sc{i_b}(i_s,:);
        end
        
        
        
        
        %%%% RATING PROBABILITIES
        
        % define method of rating probability        
        max_p = 5;
        if t_f_bySD
            max_SD_diff = max(sd_offsets);
            max_Hz_diff = max_SD_diff * sigma_e(1);
            max_Hz_diff = max_SD_diff * sigma_e_sc(1);
            
            logHz2pr = max_p * (1/max_Hz_diff);   % multiply difference in logHz by this, round, and take absolute value to get probability rating
            SD2pr    = max_p * (1/max_SD_diff);
        
        else
            max_Hz_diff = .6931;
            max_SD_diff = max_Hz_diff / sigma_e_sc(end);
%             max_SD_diff = max_Hz_diff / sigma_e(end);
            
            logHz2pr = max_p * (1/max_Hz_diff);
            SD2pr    = max_p * (1/max_SD_diff);
        end
        
        
        
        
        % calculate probability rating
        pr_noise = 1.5; % noise in converting stimulus to probability
        
        
        % - pitch diff from penultimate tone
        x  = t_f - t_pen;
        p = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;
        
        p2 = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
               
        pr_logHz_pen{i_b}(i_s,:) = [p p2];
            
        
        x  = t_f_sc - t_pen_sc;
        p = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;        
        
        pr_logHz_pen_sc{i_b}(i_s,:) = [p p2];        
        
        
        % - pitch diff from predicted final tone
        x  = t_f - t_pred;
        p = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_logHz_pred{i_b}(i_s,:) = [p p2];
        
        
        x  = t_f_sc - t_pred_sc;
        p = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( logHz2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_logHz_pred_sc{i_b}(i_s,:) = [p p2];        
        
        
        
        % - sigma_e diff from penultimate tone
        x = (t_f - t_pen) / sig_eps;
        p = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_SD_pen{i_b}(i_s,:) = [p p2];
        
        
        x = (t_f_sc - t_pen_sc) / sig_eps_sc;
        p = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_SD_pen_sc{i_b}(i_s,:) = [p p2];
        
        
        
        % - sigma_e diff from predicted final tone
        x = (t_f - t_pred) / sig_eps;
        p = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;

        p2 = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_SD_pred{i_b}(i_s,:) = [p p2];
        
        
        x = (t_f_sc - t_pred_sc) / sig_eps_sc;
        p = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p(p<1) = 1; p(p>max_p) = max_p;
        
        p2 = round( abs( SD2pr * x  + normrnd(0, pr_noise) ) );
        p2(p2<1) = 1; p2(p2>max_p) = max_p;
        
        pr_SD_pred_sc{i_b}(i_s,:) = [p p2];

    end
end


%% plot

% figure;
% ms = 10;
% style = {'ro' 'rd' 'bo' 'bd'};
% style_sc = {'r*' 'rx' 'b*' 'bx'};
% for i_b = 1:length(betas)
%     subplot(1,5,i_b); hold on;
%     
%     for i_s = 2:2
%         
% %         x    = del_SD{i_b}(i_s,:);
% %         x_sc = del_SD_sc{i_b}(i_s,:);
%         x    = del_x{i_b}(i_s,:);
%         x_sc = del_x_sc{i_b}(i_s,:);
% 
%         
%         % - pitch diff from penultimate tone
%         plot(x, pr_logHz_pen{i_b}(i_s,:), style{1}, 'MarkerSize', ms);
%         plot(x_sc, pr_logHz_pen_sc{i_b}(i_s,:), style_sc{1}, 'MarkerSize', ms);
%         
%         
% %         % - pitch diff from predicted final tone
% %         plot(x, pr_logHz_pred{i_b}(i_s,:), style{2}, 'MarkerSize', ms);
% %         plot(x_sc, pr_logHz_pred_sc{i_b}(i_s,:), style_sc{2}, 'MarkerSize', ms);
%         
%         
% %         % - sigma_e diff from penultimate tone
% %         plot(x, pr_SD_pen{i_b}(i_s,:), style{3}, 'MarkerSize', ms);
% %         plot(x_sc, pr_SD_pen_sc{i_b}(i_s,:), style_sc{3}, 'MarkerSize', ms);
%         
%         
%         % - sigma_e diff from predicted final tone
%         plot(x, pr_SD_pred{i_b}(i_s,:), style{4}, 'MarkerSize', ms);
%         plot(x_sc, pr_SD_pred_sc{i_b}(i_s,:), style_sc{4}, 'MarkerSize', ms);
%         
%         ylabel('probability rating')
%         if t_f_xaxis_bySD
%             xlabel('#SDs from expectation')
%         else
%             xlabel('logHz diff from penultimate tone')
%         end
%         title(['beta = ' num2str(betas(i_b))]);
%     end
% end
        

%% model comparison


%%% subject = logHz_pen
probs    = pr_logHz_pen;
probs_sc = pr_logHz_pen_sc;

iv_list    = {del_logHz_pen,  del_logHz_pred,  del_SD_pen,  del_SD_pred};
iv_sc_list = {del_logHz_pen_sc,  del_logHz_pred_sc,  del_SD_pen_sc,  del_SD_pred_sc};

for i_iv = 1:4

predictor    = iv_list{i_iv};
predictor_sc = iv_sc_list{i_iv};

dv = [];
iv = [];
for i_b = 1:length(betas)
    for i_s = 1:3
        dv = [dv probs{i_b}(i_s,:)];
        dv = [dv probs_sc{i_b}(i_s,:)];
        
        iv = [iv predictor{i_b}(i_s,:)];
        iv = [iv predictor_sc{i_b}(i_s,:)];
    end
end

m = [dv' iv'];
in.DATA = m;
out = MATLAB_Ordered_Probit_Estimate(in);

logL(i_iv) = out.Likelihood.LLV;

end

K = 5;
N = length(dv);

BIC = -2*logL + K*log(N);
AIC = -2*logL + 2*K*N/(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ));

logL_logHz_pen = logL;
AIC_logHz_pen  = AIC;
BIC_logHz_pen  = BIC;
AICw_logHz_pen = AICw;
BICw_logHz_pen = BICw;





%%% subject = logHz_pred
probs    = pr_logHz_pred;
probs_sc = pr_logHz_pred_sc;

iv_list    = {del_logHz_pen,  del_logHz_pred,  del_SD_pen,  del_SD_pred};
iv_sc_list = {del_logHz_pen_sc,  del_logHz_pred_sc,  del_SD_pen_sc,  del_SD_pred_sc};

for i_iv = 1:4

predictor    = iv_list{i_iv};
predictor_sc = iv_sc_list{i_iv};

dv = [];
iv = [];
for i_b = 1:length(betas)
    for i_s = 1:3
        dv = [dv probs{i_b}(i_s,:)];
        dv = [dv probs_sc{i_b}(i_s,:)];
        
        iv = [iv predictor{i_b}(i_s,:)];
        iv = [iv predictor_sc{i_b}(i_s,:)];
    end
end

m = [dv' iv'];
in.DATA = m;
out = MATLAB_Ordered_Probit_Estimate(in);

logL(i_iv) = out.Likelihood.LLV;

end

K = 5;
N = length(dv);

BIC = -2*logL + K*log(N);
AIC = -2*logL + 2*K*N/(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ));

logL_logHz_pred = logL;
AIC_logHz_pred  = AIC;
BIC_logHz_pred  = BIC;
AICw_logHz_pred = AICw;
BICw_logHz_pred = BICw;




%%% subject = SD_pen
probs    = pr_SD_pen;
probs_sc = pr_SD_pen_sc;

iv_list    = {del_logHz_pen,  del_logHz_pred,  del_SD_pen,  del_SD_pred};
iv_sc_list = {del_logHz_pen_sc,  del_logHz_pred_sc,  del_SD_pen_sc,  del_SD_pred_sc};

for i_iv = 1:4

predictor    = iv_list{i_iv};
predictor_sc = iv_sc_list{i_iv};

dv = [];
iv = [];
for i_b = 1:length(betas)
    for i_s = 1:3
        dv = [dv probs{i_b}(i_s,:)];
        dv = [dv probs_sc{i_b}(i_s,:)];
        
        iv = [iv predictor{i_b}(i_s,:)];
        iv = [iv predictor_sc{i_b}(i_s,:)];
    end
end

m = [dv' iv'];
in.DATA = m;
out = MATLAB_Ordered_Probit_Estimate(in);

logL(i_iv) = out.Likelihood.LLV;

end

K = 5;
N = length(dv);

BIC = -2*logL + K*log(N);
AIC = -2*logL + 2*K*N/(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ));

logL_SD_pen = logL;
AIC_SD_pen  = AIC;
BIC_SD_pen  = BIC;
AICw_SD_pen = AICw;
BICw_SD_pen = BICw;




%%% subject = SD_pred
probs    = pr_SD_pred;
probs_sc = pr_SD_pred_sc;

iv_list    = {del_logHz_pen,  del_logHz_pred,  del_SD_pen,  del_SD_pred};
iv_sc_list = {del_logHz_pen_sc,  del_logHz_pred_sc,  del_SD_pen_sc,  del_SD_pred_sc};

for i_iv = 1:4

predictor    = iv_list{i_iv};
predictor_sc = iv_sc_list{i_iv};

dv = [];
iv = [];
for i_b = 1:length(betas)
    for i_s = 1:3
        dv = [dv probs{i_b}(i_s,:)];
        dv = [dv probs_sc{i_b}(i_s,:)];
        
        iv = [iv predictor{i_b}(i_s,:)];
        iv = [iv predictor_sc{i_b}(i_s,:)];
    end
end

m = [dv' iv'];
in.DATA = m;
out = MATLAB_Ordered_Probit_Estimate(in);

logL(i_iv) = out.Likelihood.LLV;

end

K = 5;
N = length(dv);

BIC = -2*logL + K*log(N);
AIC = -2*logL + 2*K*N/(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ));

logL_SD_pred = logL;
AIC_SD_pred  = AIC;
BIC_SD_pred  = BIC;
AICw_SD_pred = AICw;
BICw_SD_pred = BICw;


BICw_logHz_pen;
BICw_logHz_pred;
BICw_SD_pen;
BICw_SD_pred;

if t_f_bySD
    % rows are true generating models
    % columns are recovering models
    BICw_tf_SD(:,:,i_rep) = [ ...
        BICw_logHz_pen; ...
        BICw_logHz_pred; ...
        BICw_SD_pen; ...
        BICw_SD_pred;]

else
    BICw_tf_Hz(:,:,i_rep) = [ ...
        BICw_logHz_pen; ...
        BICw_logHz_pred; ...
        BICw_SD_pen; ...
        BICw_SD_pred;]
    
end

end