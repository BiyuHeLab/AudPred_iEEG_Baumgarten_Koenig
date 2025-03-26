function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_temporal(dv_orig, p_thresh, analysis, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_temporal(dv, p_thresh, analysis, nReps)
% 
% format of dv
% 1 x nSamples cell array
% each cell holds a matrix of size (nSubs, nDVs) where nDVs correspond to the number of within-subject experimental conditions
%
% analysis is a string specifying the Matlab commands for conducting the
% statistical test on dv{i_sample}. The depdendent variable in these commands 
% must be referred to by the variable name "dv{i_sample}". The resulting
% p-value must be named "p_val" and the resulting test statistic must be
% named "test_stat".
%
% examples of constructing the analysis string
% ---------------------------------------------
% for repeated measures ANOVA:
%
%         analysis = [];
%         analysis = [analysis, '[p_val, table] = anova_rm(dv{i_sample}, ''off'');  '];
%         analysis = [analysis, 'p_val = p_val(1);  '];
%         analysis = [analysis, 'test_stat = table{2, 5};  '];



% original data
% -------------
% for each sample,
% - compute timecourse of the statistical effect across subjects
% across all samples,
% - find clusters in the original data timecourse
%
% shuffled data
% -------------
% for each sample,
% - shuffle experimental condition means randomly and independently for each subject
% across all samples,
% - compute timecourse of the statistical effect across subjects for the shuffled data set
% - find and save the max(ClusterStat) for thus shuffled set
% - over many iterations, this gives the null distribution for ClusterStat


%% analysis of original data set

nSamples = length(dv_orig);

% compute timecourse of the statistical effect across subjects
for i_sample = 1:nSamples

    
%     [pval, table] = anova_rm(dv{i_sample}, 'off');
%     p_val = pval(1);
%     t_teststat = table{2, 5};

    dv = dv_orig;
    eval(analysis);

    p_timecourse(i_sample)    = p_val;
    stat_timecourse(i_sample) = test_stat;

end

% find clusters in the original data timecourse
clusters_orig = find_temporal_clusters(stat_timecourse, p_timecourse, p_thresh);



%% analysis of repeatedly shuffled data set

% shuffle
nSubs = size(dv{1}, 1);
nDVs  = size(dv{1}, 2);

h = waitbar(0, 'Shuffling...');

for i_rep = 1:nReps

    waitbar(i_rep / nReps);

    % construct shuffle permutation for current iteration
    shuffle_ind = [];
    for i_sub = 1:nSubs
        shuffle_ind(i_sub, :) = randperm(nDVs);
    end


    % 
    for i_sample = 1:nSamples

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            ind = shuffle_ind(i_sub, :);
            dv_shuffle{i_sample}(i_sub, :) = dv_orig{i_sample}(i_sub, ind);
        end

        % compute stat topo for shuffled set
    %     [pval, table] = anova_rm(dv{i_sensor}, 'off');
    %     pval = pval(1);
    %     stat = table{2, 5};

        dv = dv_shuffle;
        eval(analysis);

        p_timecourse_shuffle(i_sample)    = p_val;
        stat_timecourse_shuffle(i_sample) = test_stat;

    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_temporal_clusters(stat_timecourse_shuffle, p_timecourse_shuffle, p_thresh);
    shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs;

end

% close(h);

%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set

for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat > abs( clusters_orig.cluster_statSum(i) )) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end