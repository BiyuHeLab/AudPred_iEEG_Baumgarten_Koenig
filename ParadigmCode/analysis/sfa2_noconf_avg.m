clear 

%% define script paths

addpath('/data/gogodisk1/brian/scripts/');
addpath(genpath('/data/gogodisk1/brian/scripts/mopttb/'))


%% load behavioral data

path = '/data/gogodisk1/amy/data/sfa_expt2_noconf/';

saveplot = 1;
path_fig = '/data/gogodisk1/amy/data/sfa_expt2_noconf/analysis/figures/';

mkdir(path_fig, 'avg');
path_fig_avg = [path_fig 'avg/'];


subs = {'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S20' 'S21' 'S22' 'S23' 'S24' 'S25'};

% % MEG subjects
% % except S1, who did behavioral task on the previous paradigm
% subs = {'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16' 'S21'};

bad = {'S5' ...   % rho = .05
       'S12' ...  % rho = .07
       'S20' ...  % rho = .01
       'S23' ...  % rho = .04
       };

subs = setdiff(subs,bad);
   
for i_sub = 1:length(subs)

sub = subs{i_sub};

switch sub
    case 'S2',  load([path 'S2/behavioral/sfa_expt2 2014-09-10 14-27-52.mat']);
    case 'S3',  load([path 'S3/behavioral/sfa_expt2 2014-09-16 10-50-03.mat']);
    case 'S4',  load([path 'S4/behavioral/sfa_expt2 2014-09-19 11-59-16.mat']);
    case 'S5',  load([path 'S5/behavioral/sfa_expt2 2014-09-19 14-53-12.mat']);
    case 'S6',  load([path 'S6/behavioral/sfa_expt2 2014-09-24 13-46-35.mat']);
    case 'S7',  load([path 'S7/behavioral/sfa_expt2 2014-09-25 10-43-38.mat']);
    case 'S8',  load([path 'S8/behavioral/sfa_expt2 2014-09-25 13-49-54.mat']);
    case 'S9',  load([path 'S9/behavioral/sfa_expt2 2014-09-26 10-22-35.mat']);
    case 'S10', load([path 'S10/behavioral/sfa_expt2 2014-10-03 14-51-59.mat']);
    case 'S11', load([path 'S11/behavioral/sfa_expt2 2014-10-08 14-36-36.mat']);
    case 'S12', load([path 'S12/behavioral/sfa_expt2 2014-10-24 11-22-42.mat']);
    case 'S13', load([path 'S13/behavioral/sfa_expt2 2014-10-28 10-41-21.mat']);
    case 'S14', load([path 'S14/behavioral/sfa_expt2 2014-11-06 11-17-11.mat']);
    case 'S15', load([path 'S15/behavioral/sfa_expt2 2014-11-10 09-58-27.mat']);
    case 'S16', load([path 'S16/behavioral/sfa_expt2 2014-11-12 10-56-04.mat']);
    case 'S17', load([path 'S17/behavioral/sfa_expt2 2014-11-22 15-39-30.mat']);
    case 'S18', load([path 'S18/behavioral/sfa_expt2 2014-11-24 11-24-09.mat']);
    case 'S20', load([path 'S20/behavioral/sfa_expt2 2014-12-08 13-21-33.mat']);
    case 'S21', load([path 'S21/behavioral/sfa_expt2 2014-12-12 10-52-51.mat']);
    case 'S22', load([path 'S22/behavioral/sfa_expt2 2014-12-19 12-58-48.mat']);
    case 'S23', load([path 'S23/behavioral/sfa_expt2 2015-01-23 14-31-36.mat']);
    case 'S24', load([path 'S24/behavioral/sfa_expt2 2015-01-28 12-26-39.mat']);
    case 'S25', load([path 'S25/behavioral/sfa_expt2 2015-01-29 12-45-32.mat']);

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


        end
        
    end
        
end

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
end


%% plot heat map of beta discrim

h  = figd(55); 
% fs = 15;
fs = 55;

imagesc(beta_matrix_avg)
axis square
colormap hot
colorbar

xlabel('response \beta', 'FontSize', fs)
ylabel('stimulus \beta', 'FontSize', fs)
% title(['avg p( response beta | stim beta ), avg rho = ' num2str(r_avg,2) ', n = ' num2str(n)], 'FontSize', fs)
title(['p( response \beta | stim \beta )'], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

% caxis([0 .71])

if saveplot
    filename = [path_fig_avg 'beta_heat_map.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end
    

%% plot cumulative d'

% figd; hold on;
figd(55); hold on;

errorbar([0 .5:.5:2], [0 dcum_avg], [0 dcum_sem], 'bo');
pp = polyfit(.5:.5:2, dcum_avg, 1);
xx = [.5 2];
% plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('stimulus \beta')
ylabel('cumulative d''');
title(['behavioral subjects, n = ' num2str(length(subs))]);
xlim([-.2 2.2])
ylim([0 2.5])

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
% figd(15,2,5);
figd(25,2,5);

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};

colors = {'ro-' 'ko-' 'go-'};


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

    hl = legend(etitle); 
    set(hl, 'box', 'off')
end

if saveplot
    filename = [path_fig_avg 'rating_prob_breakdown3.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end



%% paired t-tests for resp prob

for i_sig = 1:2
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
end