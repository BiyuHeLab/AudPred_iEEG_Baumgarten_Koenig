function stim_corr_anova_topo(sub, DV_ind, saveplot, phase_option)

if ~exist('phase_option','var') || isempty(phase_option)
    phase_option = 0;
end

% phase_option values
% 0 --> do circular correlation on phase
% 1 --> unwrap phase, detrend, and do linear correlation
% 2 --> unwrap phase, subtract expected phase, and do linear correlation


%% initialize paths for fieldtrip and personal scripts

addpath('/data/gogodisk1/brian/analysis/');
brian_ft_path


%% define input and output paths

path_base  = '/data/gogodisk1/brian/';

% sub = 'KS';

% select which MEG measure we're using for time series correlation across sensors
DVs = {'raw' ...                                     %  1
       'lowpass_5hz' 'lowpass_3hz' 'lowpass_1hz' ... %  2 -  4
       '3.33Hz_amp' '3.33Hz_phase' ...               %  5 -  6
       'delta_amp' 'delta_phase' ...                 %  7 -  8 
       'theta_amp' 'theta_phase' ...                 %  9 - 10
       'alpha_amp' 'alpha_phase' ...                 % 11 - 12
       'beta_amp' 'beta_phase' ...                   % 13 - 14
       'low_gamma_amp' 'low_gamma_phase' ...         % 15 - 16
       'mid_gamma_amp' 'mid_gamma_phase' ...         % 17 - 18
       'high_gamma_amp' 'high_gamma_phase' ...       % 19 - 20
       ...
       '3.33Hz_log_amp' ...                          % 21
       'delta_log_amp' 'theta_log_amp' ...           % 22 - 23
       'alpha_log_amp' 'beta_log_amp' ...            % 24 - 25
       'low_gamma_log_amp' 'mid_gamma_log_amp' ...   % 26 - 27
       'high_gamma_log_amp'};                        % 28
   
DV  = DVs{DV_ind};


isPhase = 0;
if length(DV) >= 5  &&  strcmp( DV(end-4:end), 'phase' )
    isPhase = 1;
end

if phase_option > 0 && isPhase == 0
    error('phase option selected for non-phase MEG DV')
end


path_in  = [path_base 'analysis/sfa_expt2_v2/stim_corr/stim_corr_results/' sub '/'];

switch phase_option
    case 1,
        filename = [sub '_stim_corr_' DV '_unwrap_detrend.mat'];
        figfile  = [sub '_stim_corr_ANOVA_within_' DV '_unwrap_detrend.png'];

    case 2,
        filename = [sub '_stim_corr_' DV '_unwrap_subtr_exp.mat'];     
        figfile  = [sub '_stim_corr_ANOVA_within_' DV '_unwrap_subtr_exp.png'];     
        
    otherwise
        filename = [sub '_stim_corr_' DV '.mat'];
        figfile  = [sub '_stim_corr_ANOVA_within_' DV '.png'];    
        
end


loadfile = [path_in  filename];

mkdir([path_base 'analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures/'], sub);
path_fig    = [path_base 'analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures/' sub '/'];
% saveplot    = 0;


%% load data

disp('loading...')
load(loadfile);
disp('done loading.')


%% get F-values for the ANOVA

anova_F = zeros(length(atab), length(factor_labels));
for i_sensor = 1:length(atab)

    t = atab{i_sensor};
    
    for i_factor = 1:length(factor_labels)
        anova_F(i_sensor, i_factor) = t{i_factor+1, 6};
    end
end
        
    
%% plot topos of p-values for each factor

% factor_labels = 
% 
%     'Beta'    'Sigma'    'Within'    'Beta*Sigma'    'Beta*Within'    'Sigma*Within'    'Beta*Sigma*Within'

plot_pval = 0;

for i_factor = 3 : 3 %1 : nFactors
    
    pval = anova_p(:, i_factor);
    Fval = anova_F(:, i_factor);

    sig_p = .05 / 271; % Bonferroni corrected p across sensors
%     sig_p = .05;
    
    % find critical value for F
%     sig_ps = pval <= sig_p;
%     sig_Fs = Fval(sig_ps);
%     crit_F = min(sig_Fs);
%     
%     if isempty(crit_F)
%         crit_F = max(Fval);
%     end
    Fval(pval > sig_p) = 0;
    
    dat.dimord = 'chan_time';
    
    if plot_pval
        dat.avg    = log10( pval );
    else
        dat.avg    = Fval;
    end
    
    dat.time   = 0;
    dat.var    = 0;
    
    load('/data/gogodisk1/brian/misc/MEG_sensor_setup_271/label.mat');
    dat.label  = label;
    
    
    cfg = [];
    cfg.layout    = 'CTF275.lay';         
    cfg.comment   = 'no';
    cfg.colorbar  = 'yes';
%     cfg.xlim      = [];
    
    if plot_pval
        cfg.zlim = [min(dat.avg)   log10( sig_p )];
    else
%         cfg.zlim = [crit_F   100]; %max(dat.avg)];
%         cfg.zlim = [0   max(Fval)]; %max(dat.avg)];
        cfg.zlim = 'maxabs';
    end
%     cfg.zlim      = [log10(sig_p / 10000)   log10( sig_p )];

%     cfg.style     = 'fill';
%     cfg.marker    = 'numbers';

    % plot
    figure;
    fontsize = 20;
    
    % title
    if plot_pval
        title({[sub ' expt 2 within vs across stim correlations'], ...
               ['using MEG DV = ' corr_within{1,1,1}.DV_label], ...
               ['log10 of ANOVA p-values for ' factor_labels{i_factor}]}, ...
                'FontSize', fontsize);
    else
        title({[sub ' expt 2 within vs across stim correlations'], ...
               ['using MEG DV = ' corr_within{1,1,1}.DV_label], ...
               ['ANOVA F values for ' factor_labels{i_factor} ' (plotting with Bonferroni correction, p < ' num2str(sig_p) ')']}, ...
                'FontSize', fontsize);
    end
    
    ft_topoplotER(cfg, dat);

    
    if saveplot
        % save plot
        coord = get(gcf, 'Position');
        set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

%         figfile = [sub '_expt2_planar_topo_' beta_list{i_beta} '_freq_' num2str(freq_plot_min,2) '_to_' num2str(freq_plot_max,2) '.png'];
        saveas(gcf, [path_fig figfile], 'png');

        delete(gcf);
    end
    
end
% keyboard
end
