function stim_corr_anova_topo_avg(dv_type, filter_type, pthresh, saveplot, input_betaID)



%% define input and output paths etc

path_data = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/stim_corr_results/';
% subs =  {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13'};
subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S16' 'S21'};

nSubs = length(subs);

mkdir(['/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/figures_avg/' dv_type '/'], filter_type);
path_fig = ['/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/figures_avg/' dv_type '/'];

switch dv_type
    case 'F'
        figfile = ['avg_stim_corr_ANOVA_within_' filter_type '_betaID=' num2str(input_betaID) '_p=' num2str(pthresh) '.png'];
        path_fig = [path_fig filter_type '/'];
%     case 'p'
%         figfile = ['avg_stim_corr_ANOVA_within_' filter_type '.png'];        
end

        

%% load and average data

within_F = zeros(length(subs), 271);
within_p = zeros(length(subs), 271);

for i_sub = 1:length(subs)
    load([path_data subs{i_sub} '/' subs{i_sub} '_stim_corr_' filter_type '_betaID=' num2str(input_betaID) '.mat']);

    for i_sensor = 1:length(atab)

        t = atab{i_sensor};
        
        within_F(i_sub, i_sensor) = t{3, 6};
        within_p(i_sub, i_sensor) = t{3, 7};

    end    
    
end

within_F_avg = mean(within_F);
within_p_avg = mean(within_p);


%% ttest against 0


cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
% cfg.zlim      = [0 100];



[h, p, ci, stats] = ttest(within_F, 1);

switch dv_type
    case 'F'
        avg_plot = within_F_avg;
        avg_plot(p > pthresh) = 0;
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

betas = [0 .5 1 1.5 2];
b_title = num2str(betas(input_betaID));

switch dv_type
    case 'F'
        if pthresh == 1
            title({['average F (n = ' num2str(nSubs) ') for within vs across stim correlations for stim beta = ' b_title], ...
                   ['using MEG DV = ' corr_within{input_betaID,1,1}.DV_label]});

        else
            title({['average F (n = ' num2str(nSubs) ') for within vs across stim correlations for stim beta = ' b_title], ...
                   ['using MEG DV = ' corr_within{input_betaID,1,1}.DV_label], ...
                   ['(plotting p < ' num2str(pthresh) ' uncorrected)']});
        end
        
    case 'p'
        title({['log10 of average p-values (n = ' num2str(nSubs) ') for within vs across stim correlations'], ...
               ['using MEG DV = ' corr_within{input_betaID,1,1}.DV_label]});
end



if saveplot
    % save plot
    coord = get(h, 'Position');
    set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

    saveas(gcf, [path_fig figfile], 'png');

    delete(h);
end
