clear 

%% define script paths

addpath('/data/gogodisk2/brian/scripts/');
addpath(genpath('/data/gogodisk2/brian/scripts/mopttb/'))


%% load behavioral data

path = '/data/gogodisk1/amy/data/sfa_expt2_noconf/';

saveplot = 0;
path_fig = '/data/gogodisk1/amy/data/sfa_expt2_noconf/analysis/figures/avg/';

mkdir(path_fig, 'MEG');
path_fig_avg = [path_fig 'MEG/'];


% subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16' 'S21'};
% subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16' 'S21'};
subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S16' 'S21'};

addpath('/data/gogodisk2/brian/analysis/sfa_expt2_v2/');

   
for i_sub = 1:length(subs)

sub = subs{i_sub};
subject_info_sfa_expt2_v2

load(si.path_behav);

stim = data.stim;

%% manual corrections to data go here

switch sub
    % S1 flipped the response probability rating
    case 'S1'
        f = data.resp_prob > 0;
        data.resp_prob(f) = 6 - data.resp_prob(f);
end
        
%% define filters

% general filter - ensure all behavioral responses were entered properly
filter_resp = data.resp_beta >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;


% manually defined filters
filter_manual = ones(1,360);

switch sub
%     case 'LS'
%         % manual filter for LS
%         % - omit trials in block 2 where r1_RT < .1 (button got stuck)
%         % - omit trials in block 3 (button was stuck for most of block)
%         manualFilter = manualFilter & data.r1_RT > .1;
%         manualFilter = manualFilter & (data.trialNum < 61 | data.trialNum > 90);
end

filter = filter_resp & filter_rt & filter_manual;



%% apply filters
% filtered data will be applied to newly defined structs "dataf" and "stimf"

% definte dataf
dfn = fieldnames(data);
for j = 1:length(dfn)
    if eval(['length(data.' dfn{j} ') == 360'])
        eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
    end
end


% define stimf
sfn = fieldnames(stim);
for j = 1:length(sfn)
    if eval(['length(stim.' sfn{j} ') == 360'])
        eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
    end
end

%% beta discrimination

betas = unique(stimf.beta);
for i_b = 1:5
    for i_r = 1:5
        
        beta_stim = betas(i_b);
        beta_resp = betas(i_r);
        
        ff = stimf.beta == beta_stim;
        beta_matrix(i_b, i_r, i_sub) = sum( dataf.resp_beta(ff) == beta_resp ) / length(dataf.resp_beta(ff));
        beta_matrix_adj(i_b, i_r) = (sum( dataf.resp_beta(ff) == beta_resp ) + 1/5) / (length(dataf.resp_beta(ff))+1);
 
        beta_cum(i_b, i_r) = sum( dataf.resp_beta(ff) <= beta_resp ) / length(dataf.resp_beta(ff));
        
        beta_nt(i_b, i_r) = sum(stimf.beta==beta_stim & dataf.resp_beta==beta_resp);
    end
end

beta_cum = cumsum(beta_matrix_adj,2);
for i_b = 2:5
    for i_r = 1:4
        beta_d(i_b-1, i_r) = norminv(beta_cum(i_b-1,i_r)) - norminv(beta_cum(i_b,i_r));
    end
end

beta_d(abs(beta_d)==Inf) = NaN;
dd = nanmean(beta_d,2);
dcum(i_sub,:) = cumsum(dd');


[r(i_sub) p(i_sub)] = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman');
% [r1 p1] = corr(stimf.sigma_e', dataf.resp_beta', 'type', 'Spearman') 

r_ts(i_sub) = r(i_sub);
p_ts(i_sub) = p(i_sub);

%% final tone probability rating

tf = unique(stimf.logf_final);

predIDs = [-1 0 1];
sigIDs  = [0 1];

for i_pred = 1:3
    for i_sig = 1:2
        
        % filter for sigma_e
        f_sig = stimf.sigma_e_ID == sigIDs(i_sig);
        
        
        % filter for expected final tone
        if predIDs(i_pred) == 0
            % medium expected tones come from beta=0 and beta=2 stim
            f_exp = (stimf.predID == predIDs(i_pred)) | stimf.betaID == 1 | stimf.betaID == 5;
        else
            % low and high expected tones never come from beta=0 or beta=2 stim
            f_exp = (stimf.predID == predIDs(i_pred)) & stimf.betaID ~= 1 & stimf.betaID ~= 5;
        end

        f = f_sig & f_exp;
        
        
        for i_tf = 1:length(tf)

            f_tone  = stimf.logf_final == tf(i_tf);
            ff      = f & f_tone;
            
            
            rp{i_sig, i_pred}(i_sub, i_tf) = mean(dataf.resp_prob(ff));
%             rp_sem{i_sig, i_pred}(i_tf) = std(dataf.resp_prob(ff)) / sqrt(sum(ff));
            
            rp_tot(i_sub, i_tf) = mean(dataf.resp_prob(f_tone));
%             rp_tot_sem(i_tf) = std(dataf.resp_prob(f_tone)) / sqrt(sum(f_tone));

            % get response probabilities by expected final tone, ignoring sigma_e
            rp_exp{i_pred}(i_sub, i_tf) = mean(dataf.resp_prob(f_exp & f_tone));
%             rp_exp{i_pred}(i_sub, i_tf) = mean(dataf.resp_prob(f_exp & f_tone & f_beta));

            for i_beta = 1:5
                f_beta = stimf.betaID == i_beta;
                rp_exp_beta{i_beta}{i_pred}(i_sub, i_tf) = mean(dataf.resp_prob(f_exp & f_tone & f_beta));
            end


        end
        
        
        mean_tf_pred(i_sub,i_pred) = mean(stimf.logf_pred(f_exp));
        
    end
end

%%%% indeces of sensitivity to expected pitch don't correlate with beta
%%%% discrimination
rp_exp_diff = rp_exp{3}(i_sub,:) - rp_exp{1}(i_sub,:);

pp = polyfit(tf, rp_exp_diff, 1);
r_pearson  = corr(tf', rp_exp_diff', 'type', 'Pearson');
r_spearman = corr(tf', rp_exp_diff', 'type', 'Spearman');

rp_exp_effect_slope(i_sub) = pp(1);
rp_exp_effect_z_pearson(i_sub) = r2z(r_pearson);
rp_exp_effect_z_spearman(i_sub) = r2z(r_spearman);


% figure; hold on;
% plot(tf, rp_exp_diff, 'bo-');
% plot(tf, tf*pp(1)+pp(2), 'r-');
% title([subs{i_sub} ', b = ' num2str(pp(1)) ', z_p = ' num2str(r2z(r_pearson)) ', z_sp = ' num2str(r2z(r_spearman))]);



%% measure of how final tone prob rating correlations with mathematically derived prediction error

nTrials = length(dataf.trialNum);
for i_trial = 1:nTrials
    series = stimf.series_f{i_trial}(1 : 33); % tone sequence up to penultimate tone
    beta = stimf.beta(i_trial);

    series_pred(i_trial) = log( find_predicted_tone(series, beta) );
    series_pitch(i_trial) = log( stimf.series_f{i_trial}(34) );

    series_pred_error(i_trial) = series_pitch(i_trial) - series_pred(i_trial);
    series_pred_error_trivial(i_trial) = series_pitch(i_trial) - log(440);
    
    rating_prob(i_trial) = dataf.resp_prob(i_trial);
    
end

% figure; hold on;
% plot(rating_prob, abs(series_pred_error), 'ro');
% plot(rating_prob, abs(series_pred_error_trivial), 'ko');
% title(subs{i_sub})

r_behav_model(i_sub) = corr(rating_prob', abs(series_pred_error)');
r_behav_model_trivial(i_sub) = corr(rating_prob', abs(series_pred_error_trivial)');

end


%% compute averages

n = length(subs);

beta_matrix_avg = mean(beta_matrix, 3);

dcum_avg = mean(dcum);
dcum_sem = std(dcum) / sqrt(n);

r_avg = z2r( mean( r2z(r) ) );


rp_tot_avg = mean(rp_tot);
rp_tot_sem = std(rp_tot) / sqrt(n);

for i_pred = 1:3
    for i_sig = 1:2
        rp_avg{i_sig, i_pred} = mean( rp{i_sig, i_pred} );
        rp_sem{i_sig, i_pred} = std( rp{i_sig, i_pred} ) / sqrt(n);
    end
    
    rp_exp_avg{i_pred} = mean( rp_exp{i_pred} );
    rp_exp_sem{i_pred} = std( rp_exp{i_pred} ) / sqrt(n);
    
    for i_beta = 1:5
        rp_exp_beta_avg{i_beta}{i_pred} = mean( rp_exp_beta{i_beta}{i_pred} );
        rp_exp_beta_sem{i_beta}{i_pred} = std( rp_exp_beta{i_beta}{i_pred} ) / sqrt(n);
    end
end


%% plot heat map of beta discrim

h  = figure; 
fs = 45;

imagesc(beta_matrix_avg)
axis square
colormap hot
colorbar

xlabel('response \beta', 'FontSize', fs)
ylabel('stimulus \beta', 'FontSize', fs)
title(['avg p( response \beta | stim \beta ), avg rho = ' num2str(r_avg,2) ', n = ' num2str(n)], 'FontSize', fs)
title('mean p( resp \beta | stim \beta )', 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

caxis([0 .5])

if saveplot
    filename = [path_fig_avg 'beta_heat_map.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end
    

%% plot cumulative d'

figd; hold on;

errorbar([0 .5:.5:2], [0 dcum_avg], [0 dcum_sem], 'bo');
pp = polyfit(.5:.5:2, dcum_avg, 1);
xx = [.5 2];
% plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('stimulus \beta')
ylabel('cumulative d''');
title(['MEG subjects, n = ' num2str(length(subs))])
xlim([-.2 2.2])

if saveplot
    filename = [path_fig_avg 'cumulative_d.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot final tone probability rating


% overall
xlab = {'220'  '277'      '349'      '554'      '698'      '880'};

figd(15,2,5);
errorbar(tf, rp_tot_avg, rp_tot_sem, 'bo-');
xlabel('final tone pitch (Hz)')
ylabel('final tone likelihood rating')
ylim([1 5])
xlim([log(200) log(900)])
set(gca,'XTick',tf)
set(gca,'XTickLabel', xlab)
title(['n = ' num2str(n)]);

if saveplot
    filename = [path_fig_avg 'rating_prob_all.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end



% sig x exp version 1        
figd(15,2,5);
ind = 0;

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
stitle = {'\sigma_\epsilon = low', '\sigma_\epsilon = high'};


for i_sig = 1:2
    for i_pred = 1:3
        ind = ind+1;
        subplot(2, 3, ind);
        errorbar(tf, rp_avg{i_sig, i_pred}, rp_sem{i_sig, i_pred}, 'bo-');
        xlabel('final tone pitch (Hz)')
        ylabel('final tone likelihood rating')
        ylim([1 5])
        xlim([log(200) log(900)])
        set(gca,'XTick',tf)
        set(gca,'XTickLabel', xlab)
        title([etitle{i_pred} ', ' stitle{i_sig}]);
    end
end

if saveplot
    filename = [path_fig_avg 'rating_prob_breakdown.png'];
    save_fig(h, filename, 1, 2, 'png');
    delete(h);
end



% sig x exp version 2
figd(15,2,5);
ind = 0;

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
stitle = {'\sigma_\epsilon = low', '\sigma_\epsilon = high'};

colors = {'ro-' 'ko-' 'go-'};

for i_sig = 1:2
    for i_pred = 1:3
        ind = ind+1;
        subplot(1, 2, i_sig); hold on;
        errorbar(tf, rp_avg{i_sig, i_pred}, rp_sem{i_sig, i_pred}, colors{i_pred});
        xlabel('final tone pitch (Hz)')
        ylabel('final tone likelihood rating')
        ylim([1 5])
        xlim([log(200) log(900)])
        set(gca,'XTick',tf)
        set(gca,'XTickLabel', xlab)
        title(stitle{i_sig});
        
        if i_sig == 1
            hl = legend(etitle); 
            set(hl, 'box', 'off')
        end
    end
end

if saveplot
    filename = [path_fig_avg 'rating_prob_breakdown2.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end



% exp collapsing across sig
figd(40,3,10);

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
etitle = {'p*_{34} = low', 'p*_{34} = med', 'p*_{34} = high'};
% etitle = {' ', ' ', ' '};

colors = {'ro-' 'ko-' 'go-'};

pred_colors = {'rs', 'ks', 'gs'};

for i_pred = 1:3
    hold on;
    errorbar(tf, rp_exp_avg{i_pred}, rp_exp_sem{i_pred}, colors{i_pred});
    xlabel('final tone pitch (Hz)')
    ylabel('final tone likelihood rating')
    ylim([1 5])
    xlim([log(200) log(900)])
    set(gca,'XTick',tf)
    set(gca,'XTickLabel', xlab)
%     title(['collapsed across \sigma_\epsilon, n = ' num2str(n)]);

    hl = legend(etitle, 'Location', 'NorthEast'); 
    set(hl, 'box', 'off')
    set(hl, 'FontSize', 25)
    
end

for i_pred = 1:3
    plot(mean(mean_tf_pred(:,i_pred)), 2, pred_colors{i_pred});
end

if saveplot
    filename = [path_fig_avg 'rating_prob_breakdown3.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end






% exp collapsing across sig, separating by beta
figd(15,2,5);

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};

colors = {'ro-' 'ko-' 'go-'};

for i_beta = 1:5
for i_pred = 1:3
    hold on;
    errorbar(tf, rp_exp_beta_avg{i_beta}{i_pred}, rp_exp_beta_sem{i_beta}{i_pred}, colors{i_pred});
    xlabel('final tone pitch (Hz)')
    ylabel('final tone likelihood rating')
    ylim([1 5])
    xlim([log(200) log(900)])
    set(gca,'XTick',tf)
    set(gca,'XTickLabel', xlab)
    title(['collapsed across \sigma_\epsilon, n = ' num2str(n)]);

    hl = legend(etitle); 
    set(hl, 'box', 'off')
end
end


%% paired t-tests for resp prob

for i_tf = 1:6

    rp_Elo_siglo = rp{1,1}(:,i_tf);
    rp_Ehi_siglo = rp{1,3}(:,i_tf);

    [h p_siglo(i_tf)] = ttest(rp_Elo_siglo, rp_Ehi_siglo);

    rp_Elo_sighi = rp{2,1}(:,i_tf);
    rp_Ehi_sighi = rp{2,3}(:,i_tf);

    [h p_sighi(i_tf)] = ttest(rp_Elo_sighi, rp_Ehi_sighi);

    rp_Elo = rp_exp{1}(:,i_tf);
    rp_Ehi = rp_exp{3}(:,i_tf);

    [h p_exp(i_tf)] = ttest(rp_Elo, rp_Ehi);

end
