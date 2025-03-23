function stim_corr_anova_topo_avg(dv_type, filter_type, pthresh, maxF, Fthresh, saveplot)



%% define input and output paths etc

path_data = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/stim_corr_results/';
% subs =  {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13'};
subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S16' 'S21'};

nSubs = length(subs);

mkdir(['/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures_avg/' dv_type '/'], filter_type);
path_fig = ['/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures_avg/' dv_type '/'];

switch dv_type
    case 'F'
        figfile = ['avg_stim_corr_ANOVA_within_' filter_type '_p=' num2str(pthresh) '_maxF=' num2str(maxF) '_Fthresh=' num2str(Fthresh) '.png'];
        path_fig = [path_fig filter_type '/'];
    case 'p'
        figfile = ['avg_stim_corr_ANOVA_within_' filter_type '.png'];        
end

        
figfile_rho = ['stim_corr_' filter_type '_rho_corr_p=' num2str(pthresh) '.png'];
figfile_dprime = ['stim_corr_' filter_type '_dprime_corr_p=' num2str(pthresh) '.png'];


%% load and average data

within_F = zeros(length(subs), 271);
within_p = zeros(length(subs), 271);

for i_sub = 1:length(subs)
    load([path_data subs{i_sub} '/' subs{i_sub} '_stim_corr_' filter_type '.mat']);

    for i_sensor = 1:length(atab)

        t = atab{i_sensor};
        
        within_F(i_sub, i_sensor) = t{4, 6};
        within_p(i_sub, i_sensor) = t{4, 7};

    end    
    
end

within_F_avg = mean(within_F);
within_p_avg = mean(within_p);


%% ttest against 0


cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
cfg.zlim      = [0 maxF];


    cfg.highlight          = 'on';
    cfg.highlightchannel   = find(within_F_avg > Fthresh); % Nx1 cell-array with selection of channels, or vector containing channel indices see FT_CHANNELSELECTION
%     cfg.highlightsymbol    = '.'; %highlight marker symbol (default = 'o')
% %     cfg.highlightcolor     = highlight marker color (default = [0 0 0] (black))
%     cfg.highlightsize      = 40; %highlight marker size (default = 6)
%     cfg.highlightfontsize  = 40; % highlight marker size (default = 8)

cfg.highlightsymbol  = 'o';
cfg.highlightcolor   = [1 0 0];
cfg.highlightsize    = 9;

[h, p, ci, stats] = ttest(within_F, 1);

switch dv_type
    case 'F'
        avg_plot = within_F_avg;
        avg_plot(p > pthresh) = 0;
%         avg_plot(within_F_avg < Fthresh) = 0;
    case 'p'
        avg_plot = log10( within_p_avg );
end
        
        
dat = make_ft_struct(avg_plot', 'timelock');


%     cfg.zlim      = [0  500];
%     cfg.style     = 'fill';
%     cfg.marker    = 'labels';
%     cfg.marker    = 'numbers';

% plot
h = figure;
ft_topoplotER(cfg, dat);


switch dv_type
    case 'F'
        if pthresh == 1
            title({['average F (n = ' num2str(nSubs) ') for within vs across stim correlations'], ...
                   ['using MEG DV = ' corr_within{1,1,1}.DV_label]});

        else
            title({['average F (n = ' num2str(nSubs) ') for within vs across stim correlations'], ...
                   ['using MEG DV = ' corr_within{1,1,1}.DV_label], ...
                   ['(plotting p < ' num2str(pthresh) ' uncorrected)']});
        end
        
    case 'p'
        title({['log10 of average p-values (n = ' num2str(nSubs) ') for within vs across stim correlations'], ...
               ['using MEG DV = ' corr_within{1,1,1}.DV_label]});
end


if saveplot
    % save plot
    coord = get(h, 'Position');
    set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

    saveas(gcf, [path_fig figfile], 'png');
    
    savefig(h, [path_fig figfile '.fig']);

    delete(h);
end


return



%% task performance


%% load behavioral performance data

for i_sub = 1:length(subs)
    
    % load data
    sub = subs{i_sub};
    
    subject_info_sfa_expt2_v2    
    load(si.path_behav);
    
    % filter
    filter_resp = data.resp_beta >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
    filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;

    filter = filter_resp & filter_rt;    

    
    % apply filters
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
    
    
    % Spearman rho for stim beta and trend strength rating
    rho(i_sub) = corr(stimf.beta', dataf.resp_beta', 'type', 'Spearman');
    
    
    
    
    % cumulative d'
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
    
    beta_d(abs(beta_d)==Inf) = NaN;
    dd = nanmean(beta_d,2);
    dcum = cumsum(dd');
    
    dprime(i_sub) = dcum(end);
    
end


%% correlate rho and phase dissim at each sensor

nSensors = 271;
for i_sensor = 1:nSensors
    [F_rho_corr(i_sensor) p_rho(i_sensor)] = corr(rho', within_F(:, i_sensor));
    [F_dprime_corr(i_sensor) p_dprime(i_sensor)] = corr(dprime', within_F(:, i_sensor));
    
%     reg = regstats(rho', dissim(:, i_sensor), 'linear');
%     behav_corr(i_sensor) = reg.tstat.beta(2);
%     p(i_sensor) = reg.tstat.pval(2);
end
    
    
%% plot the topo for Spearman rho

F_plot = F_rho_corr;
F_plot(p_rho > pthresh) = 0;

clear dat
dat = make_ft_struct(F_plot', 'freq');



cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';


% plot
h = figure;
ft_topoplotER(cfg, dat);


fontsize = 20;

if pthresh == 1
    title({['stim tracking (within vs across F value)'], ...
           ['correlation with task performance (Spearman''s rho) (n = ' num2str(nSubs) ')'], ...
           ['using MEG DV = ' corr_within{1,1,1}.DV_label]}, ...
            'FontSize', fontsize);
else
    title({['stim tracking (within vs across F value)'], ...
           ['correlation with task performance (Spearman''s rho) (p < ' num2str(pthresh) ' uncorrected) (n = ' num2str(nSubs) ')'], ...
           ['using MEG DV = ' corr_within{1,1,1}.DV_label]}, ...
            'FontSize', fontsize);    
end    


% keyboard

if saveplot
    % save plot
    coord = get(h, 'Position');
    set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

    saveas(gcf, [path_fig figfile_rho], 'png');

    delete(h);
end






%% plot the topo for d'

F_plot = F_dprime_corr;
F_plot(p_dprime > pthresh) = 0;

clear dat
dat = make_ft_struct(F_plot', 'freq');



cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';


% plot
h = figure;
ft_topoplotER(cfg, dat);


fontsize = 20;

if pthresh == 1
    title({['stim tracking (within vs across F value)'], ...
           ['correlation with task performance (cumulative d'') (n = ' num2str(nSubs) ')'], ...
           ['using MEG DV = ' corr_within{1,1,1}.DV_label]}, ...
            'FontSize', fontsize);
else
    title({['stim tracking (within vs across F value)'], ...
           ['correlation with task performance (cumulative d'') (p < ' num2str(pthresh) ' uncorrected) (n = ' num2str(nSubs) ')'], ...
           ['using MEG DV = ' corr_within{1,1,1}.DV_label]}, ...
            'FontSize', fontsize);    
end    


% keyboard

if saveplot
    % save plot
    coord = get(h, 'Position');
    set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

    saveas(gcf, [path_fig figfile_dprime], 'png');

    delete(h);
end