clear

%% initialize paths for fieldtrip and personal scripts

addpath('/data/gogodisk2/brian/analysis/');
brian_ft_path


%% define input and output paths, and cluster analysis settings

% inputs
subs = {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14'};

% subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16'};
subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S16' 'S21'};

nSubs = length(subs);

DVs = {'raw' 'lowpass_1hz' 'lowpass_3hz' 'lowpass_5hz'};
DV_titles = {'raw' 'lowpass 1 Hz' 'lowpass 3 Hz' 'lowpass 5 Hz'};
DV_ind = 1;

DV  = DVs{DV_ind};
DV_title = DV_titles{DV_ind};


path_data = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/stim_corr_results/';
% path_fig  = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures_avg/cluster_within_x_beta_within_only/';
% path_data_cluster = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/stim_corr_results/cluster_within_x_beta_within_only/';
path_fig  = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/figures_avg/cluster_within_x_beta/';
path_data_cluster = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr/stim_corr_results/cluster_within_x_beta/';

mkdir(path_fig);
mkdir(path_data_cluster);



% plotting
saveplot = 1;


% cluster analysis settings
p_thresh = .05;
nReps    = 100;

% cluster analysis output
data_cluster_filename = [DV '_cluster_p_thresh=' num2str(p_thresh) '_nReps=' num2str(nReps) '.mat'];
figfile               = [DV '_cluster_topo_p_thresh=' num2str(p_thresh) '_nReps=' num2str(nReps) '.png'];

load_preexisting_data = 0;

nSensors = 271;

%% load data and process data

if load_preexisting_data
    load([path_data_cluster, data_cluster_filename]);

else 
    
    for i_sub = 1:length(subs)

        % load
        sub = subs{i_sub};
        load([path_data sub '/' sub '_stim_corr_' DV '.mat']);


        % construct the sensor x subject x condition data set
        nSensors = size(zvals, 1);
        for i_sensor = 1:nSensors
            zvals_s  = zvals(i_sensor, :);
            within_s = Within(i_sensor, :);
            beta_s   = Beta(i_sensor, :);
            
            ind = 0;
            for i_within = 0:1
                for i_beta = 1:5
                    ind = ind + 1;
                    dv{i_sensor}(i_sub, ind) = mean(zvals_s(within_s == i_within & beta_s == i_beta));
                    f_sub(i_sub, ind) = i_sub;
                    f_within(i_sub, ind) = i_within;
                    f_beta(i_sub, ind) = i_beta;
                end
            end
        end

    end

    f_names = {'within', 'beta'};


    %% do repeated measures ANOVA on each sensor and time window.


    analysis = [];
    analysis = [analysis, 'stats = rm_anova2( dv{i_sensor}(:), f_sub(:), f1(:), f2(:), f_names );  '];
    analysis = [analysis, 'p_val = stats{4,6};  '];
    analysis = [analysis, 'test_stat = stats{4, 5};  '];
%     analysis = [analysis, 'p_val = stats{2,6};  '];
%     analysis = [analysis, 'test_stat = stats{2, 5};  '];

    t0 = GetSecs;
    [clusters_orig, xx, shuffleMaxStat] = cluster_test_rm_anova2(dv, f_sub, f_within, f_beta, f_names, p_thresh, analysis, nReps);
    tf = GetSecs - t0;


    save([path_data_cluster data_cluster_filename], 'clusters_orig', 'shuffleMaxStat', 'nReps', 'analysis', 'tf', '-v7.3');
end


%% plot topo

topo_plot = zeros(1,nSensors);
nClusterPlot = 0;

for i_cluster = 1:clusters_orig.nClusters

    cluster_p = clusters_orig.cluster_pval(i_cluster);
    if cluster_p < .05
        nClusterPlot = nClusterPlot + 1;

        cluster_ind  = clusters_orig.cluster_sensors{i_cluster};

        topo_plot(cluster_ind) = clusters_orig.inputs.topo_stat(cluster_ind);
        cluster_sizes(nClusterPlot) = clusters_orig.cluster_size(i_cluster);
    end

end

dat = make_ft_struct(topo_plot', 'timelock');

cfg = [];
cfg.layout    = 'CTF275.lay';         
cfg.colorbar  = 'yes';
cfg.comment   = 'no';
if nClusterPlot == 0
    cfg.zlim = [0 1];
end
% cfg.zlim = [0 80];

figure;
fontsize = 15;


if nClusterPlot == 0
    cluster_size_text = 'N/A';
else
    cluster_size_text = [];
    for i_cluster = 1:nClusterPlot
        cluster_size_text = [cluster_size_text num2str(cluster_sizes(i_cluster))];
        if i_cluster < nClusterPlot
            cluster_size_text = [cluster_size_text ', '];
        end
    end
end


title({['Effect of Within x Beta on stimulus correlation (n = ' num2str(nSubs) ')'], ...
       ['for ' DV_title ' timecourse'], ...
       ['plotting ANOVA F values for ' num2str(nClusterPlot) ' clusters with p < .05'], ...
       ['p-threshold for cluster definition = ' num2str(p_thresh) ', cluster sizes = ' cluster_size_text]}, ...
        'FontSize', fontsize);

ft_topoplotER(cfg, dat);

if saveplot
    % save plot
    coord = get(gcf, 'Position');
    set(gcf, 'Position', [1 1 coord(3)*2 coord(4)*2]);

    saveas(gcf, [path_fig figfile], 'png');

    delete(gcf);
end
    