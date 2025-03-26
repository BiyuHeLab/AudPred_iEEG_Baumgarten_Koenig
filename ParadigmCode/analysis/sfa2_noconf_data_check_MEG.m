clear 

%% initialize paths for fieldtrip and personal scripts

addpath('/data/gogodisk1/brian/analysis/');
brian_ft_path


%% initialize paths etc for analysis of this subject

subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16' 'S17'};

sub = 'S17';
addpath('/data/gogodisk1/brian/analysis/sfa_expt2_v2/');
subject_info_sfa_expt2_v2


%% load behavioral data

load(si.path_behav);

stim = data.stim;

path = '/data/gogodisk1/amy/data/sfa_expt2_noconf/';

% subs = {'S1' 'S2' 'S3' 'S4'};
% 
% sub = subs{1};
% 
% switch sub
%     case 'S2',  load([path 'S2/behavioral/sfa_expt2 2014-09-10 14-27-52.mat']);
%     case 'S3',  load([path 'S3/behavioral/sfa_expt2 2014-09-16 10-50-03.mat']);
%     case 'S4',  load([path 'S4/behavioral/sfa_expt2 2014-09-19 11-59-16.mat']);
%     case 'S5',  load([path 'S5/behavioral/sfa_expt2 2014-09-19 14-53-12.mat']);
%     case 'S6',  load([path 'S6/behavioral/sfa_expt2 2014-09-24 13-46-35.mat']);
%     case 'S7',  load([path 'S7/behavioral/sfa_expt2 2014-09-25 10-43-38.mat']);
%     case 'S8',  load([path 'S8/behavioral/sfa_expt2 2014-09-25 13-49-54.mat']);
%     case 'S9',  load([path 'S9/behavioral/sfa_expt2 2014-09-26 10-22-35.mat']);
%     case 'S10', load([path 'S10/behavioral/sfa_expt2 2014-10-03 14-51-59.mat']);
%        
% end

saveplot = 1;
path_fig = '/data/gogodisk1/amy/data/sfa_expt2_noconf/analysis/figures/';

mkdir([path_fig sub '/'], 'MEG');
path_fig_sub = [path_fig sub '/MEG/'];


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
        beta_matrix(i_b, i_r) = sum( dataf.resp_beta(ff) == beta_resp ) / length(dataf.resp_beta(ff));
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

[r p] = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman');
[r1 p1] = corr(stimf.sigma_e', dataf.resp_beta', 'type', 'Spearman');



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
            rp{i_sig, i_pred}(i_tf) = mean(dataf.resp_prob(ff));
            rp_sem{i_sig, i_pred}(i_tf) = std(dataf.resp_prob(ff)) / sqrt(sum(ff));
            
            rp_tot(i_tf) = mean(dataf.resp_prob(f_tone));
            rp_tot_sem(i_tf) = std(dataf.resp_prob(f_tone)) / sqrt(sum(f_tone));

            
            % get response probabilities by expected final tone, ignoring sigma_e
            rp_exp{i_pred}(i_tf) = mean(dataf.resp_prob(f_exp & f_tone));
            rp_exp_sem{i_pred}(i_tf) = std(dataf.resp_prob(f_exp & f_tone)) / sqrt(sum(f_exp & f_tone));
            
        end
        
    end
end



%% plot heat map of beta discrim

h  = figd; 
fs = 15;

imagesc(beta_matrix)
axis square
colormap hot
colorbar

xlabel('response beta', 'FontSize', fs)
ylabel('stimulus beta', 'FontSize', fs)
title(['p( response beta | stim beta ), rho = ' num2str(r,2)], 'FontSize', fs)

set(gca,'XTick',1:5)
set(gca,'YTick',1:5)
set(gca,'XTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'YTickLabel',{'0' '0.5' '1' '1.5' '2'})
set(gca,'FontSize',fs)

caxis([0 .71])

if saveplot
    filename = [path_fig_sub 'beta_heat_map.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end
    

%% plot cumulative d'

figd; hold on;

beta_d(abs(beta_d)==Inf) = NaN;
dd = nanmean(beta_d,2);
dcum = cumsum(dd');

plot([0 .5:.5:2], [0 dcum], 'bo');
pp = polyfit(.5:.5:2, dcum, 1);
xx = [.5 2];
plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('beta')
ylabel('cumulative d''');
xlim([-.2 2.2])

if saveplot
    filename = [path_fig_sub 'cumulative_d.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot responses over time

figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.resp_prob, 'b.');
xlabel('trial')
ylabel('final tone prob rating')
xlim([1 360])
set(gca,'YTick',[1:5])
ylim([.5 5.5])
title([num2str(length(dataf.trialNum)) ' of 360 trials included'])

subplot(2,1,2);
plot(dataf.trialNum, dataf.resp_beta, 'b.');
xlabel('trial')
ylabel('beta response')
xlim([1 360])
set(gca,'YTick',[0:.5:2])
ylim([-.25 2.25])

if saveplot
    filename = [path_fig_sub 'timecourse_responses.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot RTs over time

figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.r1_RT, 'b-');
xlabel('trial')
ylabel('RT of final tone prob rating')
xlim([1 360])
title([num2str(length(dataf.trialNum)) ' of 360 trials included'])
ylim([0 5])

subplot(2,1,2);
plot(dataf.trialNum, dataf.r2_RT, 'b-');
xlabel('trial')
ylabel('RT of beta response')
xlim([1 360])
ylim([0 5])

if saveplot
    filename = [path_fig_sub 'timecourse_RTs.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot performance over time

figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.correct_beta, 'b.');
xlabel('trial')
ylabel('beta response accuracy')
xlim([1 360])
ylim([-.25 1.25])
set(gca,'YTick', [0 1]);
title([num2str(length(dataf.trialNum)) ' of 360 trials included'])

subplot(2,1,2);
plot(dataf.trialNum, dataf.diff_beta, 'b.');
xlabel('trial')
ylabel('beta response error')
xlim([1 360])
ylim([-2.25 2.25])
set(gca,'YTick', [-2:2]);

if saveplot
    filename = [path_fig_sub 'timecourse_performance.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end



%% plot stim

if saveplot
    filename = [path_fig_sub 'stim.png'];
    h = sfa_expt2_plot_stim(stim, filename);
    delete(h);
else
    sfa_expt2_plot_stim(stim);
end


%% plot final tone probability rating

xlab = {'220'  '277'      '349'      '554'      '698'      '880'};

figd(15,2,5);
errorbar(tf, rp_tot, rp_tot_sem, 'bo-');
xlabel('final tone pitch (Hz)')
ylabel('final tone likelihood rating')
ylim([1 5])
xlim([log(200) log(900)])
set(gca,'XTick',tf)
set(gca,'XTickLabel', xlab)
title('for all trials');

if saveplot
    filename = [path_fig_sub 'rating_prob_all.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end

        
figd(15,2,5);
ind = 0;

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
stitle = {'\sigma_\epsilon = low', '\sigma_\epsilon = high'};


for i_sig = 1:2
    for i_pred = 1:3
        ind = ind+1;
        subplot(2, 3, ind);
        errorbar(tf, rp{i_sig, i_pred}, rp_sem{i_sig, i_pred}, 'bo-');
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
    filename = [path_fig_sub 'rating_prob_breakdown.png'];
    save_fig(h, filename, 1, 2, 'png');
    delete(h);
end




figd(15,2,5);
ind = 0;

etitle = {'E[pitch] = low', 'E[pitch] = med', 'E[pitch] = high'};
stitle = {'\sigma_\epsilon = low', '\sigma_\epsilon = high'};

colors = {'ro-' 'ko-' 'go-'};

for i_sig = 1:2
    for i_pred = 1:3
        ind = ind+1;
        subplot(1, 2, i_sig); hold on;
        errorbar(tf, rp{i_sig, i_pred}, rp_sem{i_sig, i_pred}, colors{i_pred});
        xlabel('final tone pitch (Hz)')
        ylabel('final tone likelihood rating')
        ylim([1 5])
        xlim([log(200) log(900)])
        set(gca,'XTick',tf)
        set(gca,'XTickLabel', xlab)
        title(stitle{i_sig});
        
        if i_sig == 1
            hl = legend(etitle, 'Location', 'SouthWest'); 
            set(hl, 'box', 'off')
        end
    end
end

if saveplot
    filename = [path_fig_sub 'rating_prob_breakdown2.png'];
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
    errorbar(tf, rp_exp{i_pred}, rp_exp_sem{i_pred}, colors{i_pred});
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
    filename = [path_fig_sub 'rating_prob_breakdown3.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end