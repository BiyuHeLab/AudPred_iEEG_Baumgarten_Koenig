clear 

%% load behavioral data

toneDurs = {'0.15' '0.3' '0.6' 'all'};
toneDur = toneDurs{3};

path_data = '/data/gogodisk2/brian/data/behavioral/sfa_expt3/data/';

subs = {'BM'};

sub = subs{1};

switch sub
    case 'BM',  load([path_data 'brian pilot data/sfa_expt3 2016-05-18 16-09-43.mat']);
end

saveplot = 1;
path_fig = '/data/gogodisk2/brian/data/behavioral/sfa_expt3/analysis/figures/';

path_fig_sub = [path_fig sub '/'];

switch toneDur
    case '0.15'
        path_fig_sub = [path_fig_sub 'toneDur=0.15/'];

    case '0.3'
        path_fig_sub = [path_fig_sub 'toneDur=0.3/'];

    case '0.6'
        path_fig_sub = [path_fig_sub 'toneDur=0.6/'];
end

mkdir(path_fig_sub);


%% check counterbalancing

trials_by_block = reshape(1:324, [27 12]);
for i_block = 1:12
    
    trials = trials_by_block(:, i_block);
    
    count = zeros(3,3,3);
    
    for i_trial = trials
        count(stim.betaID(i_trial), stim.toneDurID(i_trial), stim.predID(i_trial)+2) = ...  
            count(stim.betaID(i_trial), stim.toneDurID(i_trial), stim.predID(i_trial)+2) + 1;
    end
    
    if all(count(:) == 1)
        block_check(i_block) = 1;
    else
        block_check(i_block) = 0;
        disp('bad counterbalancing')
    end
end
    




%% define filters

% general filter - ensure all behavioral responses were entered properly
filter_resp = data.resp_beta >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;

switch toneDur
    case '0.15'
        filter_toneDur = stim.toneDur == 0.15;

    case '0.3'
        filter_toneDur = stim.toneDur == 0.3;

    case '0.6'
        filter_toneDur = stim.toneDur == 0.6;
        
    otherwise
        filter_toneDur = ones(1,324);

end

% manually defined filters
filter_manual = ones(1,324);

switch sub
%     case 'LS'
%         % manual filter for LS
%         % - omit trials in block 2 where r1_RT < .1 (button got stuck)
%         % - omit trials in block 3 (button was stuck for most of block)
%         manualFilter = manualFilter & data.r1_RT > .1;
%         manualFilter = manualFilter & (data.trialNum < 61 | data.trialNum > 90);
end

filter = filter_resp & filter_rt & filter_toneDur & filter_manual;



%% apply filters
% filtered data will be applied to newly defined structs "dataf" and "stimf"

% definte dataf
dfn = fieldnames(data);
for j = 1:length(dfn)
    if eval(['length(data.' dfn{j} ') == 324'])
        eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
    end
end


% define stimf
sfn = fieldnames(stim);
for j = 1:length(sfn)
    if eval(['length(stim.' sfn{j} ') == 324'])
        eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
    end
end

%% beta discrimination

betas = unique(stimf.beta);
for i_b = 1:3
    for i_r = 1:3
        
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
for i_b = 2:3
    for i_r = 1:2
        beta_d(i_b-1, i_r) = norminv(beta_cum(i_b-1,i_r)) - norminv(beta_cum(i_b,i_r));
    end
end

[r p] = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman') 


%% final tone probability rating

tf = unique(stimf.logf_final);

predIDs = [-1 0 1];

for i_pred = 1:3
        
    f_exp = stimf.predID == predIDs(i_pred);

    for i_tf = 1:length(tf)

        f_tone  = stimf.logf_final == tf(i_tf);
        ff      = f_exp & f_tone;

        rp_tot(i_tf) = mean(dataf.resp_prob(f_tone));
        rp_tot_sem(i_tf) = std(dataf.resp_prob(f_tone)) / sqrt(sum(f_tone));

        % get response probabilities by expected final tone
        rp_exp{i_pred}(i_tf) = mean(dataf.resp_prob(f_exp & f_tone));
        rp_exp_sem{i_pred}(i_tf) = std(dataf.resp_prob(f_exp & f_tone)) / sqrt(sum(f_exp & f_tone));

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
title({'p( response beta | stim beta )', ['rho = ' num2str(r,2) ', p = ' num2str(p,3)]}, 'FontSize', fs)

set(gca,'XTick',1:3)
set(gca,'YTick',1:3)
set(gca,'XTickLabel',{'0.5' '0.99' '1.5'})
set(gca,'YTickLabel',{'0.5' '0.99' '1.5'})
set(gca,'FontSize',fs)

caxis([0 .71])

if saveplot
    filename = [path_fig_sub 'beta_heat_map.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end
    

%% plot cumulative d'

h = figd; hold on;

beta_d(abs(beta_d)==Inf) = NaN;
dd = nanmean(beta_d,2);
dcum = cumsum(dd');

plot([0.5 0.99 1.5], [0 dcum], 'bo');
% pp = polyfit(.5:.5:2, dcum, 1);
% xx = [.5 2];
% plot(xx, xx*pp(1)+pp(2), 'b-'); 
xlabel('beta')
ylabel('cumulative d''');
xlim([0.3 1.7])

if saveplot
    filename = [path_fig_sub 'cumulative_d.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot responses over time

h = figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.resp_prob, 'b.');
xlabel('trial')
ylabel('final tone prob rating')
xlim([1 324])
set(gca,'YTick',[1:5])
ylim([.5 5.5])
title([num2str(length(dataf.trialNum)) ' of 324 trials included'])

subplot(2,1,2);
plot(dataf.trialNum, dataf.resp_beta, 'b.');
xlabel('trial')
ylabel('beta response')
xlim([1 324])
set(gca,'YTick',[0.5 0.99 1.5])
ylim([0.3 1.7])

if saveplot
    filename = [path_fig_sub 'timecourse_responses.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot RTs over time

h = figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.r1_RT, 'b-');
xlabel('trial')
ylabel('RT of final tone prob rating')
xlim([1 324])
title([num2str(length(dataf.trialNum)) ' of 324 trials included'])
ylim([0 5])

subplot(2,1,2);
plot(dataf.trialNum, dataf.r2_RT, 'b-');
xlabel('trial')
ylabel('RT of beta response')
xlim([1 324])
ylim([0 5])

if saveplot
    filename = [path_fig_sub 'timecourse_RTs.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end


%% plot performance over time

h = figd(15,1); hold on;
subplot(2,1,1);
plot(dataf.trialNum, dataf.correct_beta, 'b.');
xlabel('trial')
ylabel('beta response accuracy')
xlim([1 324])
ylim([-.25 1.25])
set(gca,'YTick', [0 1]);
title([num2str(length(dataf.trialNum)) ' of 324 trials included'])

subplot(2,1,2);
plot(dataf.trialNum, dataf.diff_beta, 'b.');
xlabel('trial')
ylabel('beta response error')
xlim([1 324])
ylim([-1.25 1.25])
set(gca,'YTick', [-1:1]);

if saveplot
    filename = [path_fig_sub 'timecourse_performance.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end



%% plot stim

if saveplot
    filename = [path_fig_sub 'stim.png'];
    h = sfa_expt3_plot_stim(stim, filename);
    delete(h);
else
    sfa_expt3_plot_stim(stim);
end


%% plot final tone probability rating

xlab = {'220'  '277'      '349'      '554'      '698'      '880'};

h = figd(15,2,5);
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



h = figd(25,2,5);

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
    filename = [path_fig_sub 'rating_prob_by_EV.png'];
    save_fig(h, filename, 1, 1, 'png');
    delete(h);
end