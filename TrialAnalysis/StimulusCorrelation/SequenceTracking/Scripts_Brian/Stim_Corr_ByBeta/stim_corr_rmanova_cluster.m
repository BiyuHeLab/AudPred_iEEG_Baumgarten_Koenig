function stim_corr_rmanova_cluster(input_betaID)

%% initialize paths for fieldtrip and personal scripts

addpath('/data/gogodisk2/brian/analysis/');
brian_ft_path


%% define input and output paths, and cluster analysis settings

% inputs
subs = {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14'};

% subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S13' 'S16'};
subs = {'S1' 'S2' 'S3' 'S4' 'S6' 'S8' 'S10' 'S14' 'S15' 'S16' 'S21'};

nSubs = length(subs);

% input_betaID = 1;

DVs = {'raw' 'lowpass_1hz' 'lowpass_3hz' 'lowpass_5hz'};
DV_titles = {'raw' 'lowpass 1 Hz' 'lowpass 3 Hz' 'lowpass 5 Hz'};
DV_ind = 1;

DV  = DVs{DV_ind};
DV_title = DV_titles{DV_ind};


path_data = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/stim_corr_results/';
path_fig  = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/figures_cluster/';
path_data_cluster = '/data/gogodisk2/brian/analysis/sfa_expt2_v2/stimulus_tracking/stim_corr_by_beta/stim_corr_results/cluster/';




% plotting
saveplot = 1;


% cluster analysis settings
p_thresh = .05;
nReps    = 1000;

% cluster analysis output
data_cluster_filename = [DV '_cluster_p_thresh=' num2str(p_thresh) '_nReps=' num2str(nReps) '_nSubs=' num2str(nSubs) '_betaID=' num2str(input_betaID) '.mat'];
figfile               = [DV '_cluster_topo_p_thresh=' num2str(p_thresh) '_nReps=' num2str(nReps) '_nSubs=' num2str(nSubs) '_betaID=' num2str(input_betaID) '.png'];

load_preexisting_data = 1;

nSensors = 271;

%% load data and process data

if load_preexisting_data
    load([path_data_cluster, data_cluster_filename]);

else 
    for i_sub = 1:length(subs)

        % load
        sub = subs{i_sub};
        load([path_data sub '/' sub '_stim_corr_' DV '_betaID=' num2str(input_betaID) '.mat']);


        % construct the sensor x subject x condition data set
        nSensors = size(zvals, 1);
        for i_sensor = 1:nSensors
            zvals_s  = zvals(i_sensor, :);
            within_s = Within(i_sensor, :);

            for i_within = 0:1
                dv{i_sensor}(i_sub, i_within+1) = mean(zvals_s(within_s == i_within));
            end
        end

    end



    %% do repeated measures ANOVA on each sensor and time window.


    analysis = [];
    analysis = [analysis, '[p_val, table] = anova_rm(dv{i_sensor}, ''off'');  '];
    analysis = [analysis, 'p_val = p_val(1);  '];
    analysis = [analysis, 'test_stat = table{2, 5};  '];

    t0 = GetSecs;
    [clusters_orig, xx, shuffleMaxStat] = cluster_test(dv, p_thresh, analysis, nReps);
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
cfg.zlim = [0 80];

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

betas = [0 .5 1 1.5 2];
b_title = num2str(betas(input_betaID));
title({['Effect of Within vs Across stimulus correlation (n = ' num2str(nSubs) ')'], ...
       ['on ' DV_title ' timecourse for stim beta = ' b_title], ...
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
    