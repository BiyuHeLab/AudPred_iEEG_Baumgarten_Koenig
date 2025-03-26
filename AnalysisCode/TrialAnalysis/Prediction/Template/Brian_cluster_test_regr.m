function [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test_regr(dv_orig, iv_orig, p_thresh, fn_handle, nReps)
% [clusters_orig, clusters_shuffle, shuffleMaxStat] = cluster_test(dv, iv, p_thresh, fn_handle, nReps)
% 
% dv is a cell array of format
% dv{sensor #}(subject #, condition #)
%
% where condition refers to data a within-subject experimental condition
%
% iv is a cell array of format
% iv(IV #, condition #)
%
% where IV# indicates the possibility of multiple IVs for a multiple
% regression.
%
% fn_handle is a function handle to a function for analyzing the across-condition 
% regression b/t dv and iv for each subject.
%
% The function referred to by fn_handle must be of the form
%
% [test_stat, p_val] = fun(dv, iv)
%
% where test_stat is the relevant summary statistic for the dv/iv relationship
% (e.g. regression beta or Pearson correlation z-score) and p_val is its 
% corresponding p-value.
%
% Shuffling takes place across DVs within each individual subject data.
% Then the summary correlation statistic test_stat is subjected to an 
% across-subject one-sample t-test. The t value and p-value from this
% t-test are then used for cluster based analysis. Clusters are defined as
% sets of neighboring sensors such that every sensor in the cluster has a
% p-value less than p_thresh.



% original data
% -------------
% for each sensor,
% - compute topo of the statistical effect across subjects
% across all sensors,
% - find clusters in the original data topo
%
% shuffled data
% -------------
% for each sensor,
% - shuffle experimental condition means randomly and independently for each subject
% across all sensors,
% - compute topo of the statistical effect across subjects for the shuffled data set
% - find and save the max(ClusterStat) for thus shuffled set
% - over many iterations, this gives the null distribution for ClusterStat


%% analysis of original data set

nSensors = length(dv_orig);
nSubs    = size(dv_orig{1}, 1);
nCond    = size(dv_orig{1}, 2);

% compute topo of the statistical effect across subjects
for i_sensor = 1:nSensors

    for i_sub = 1:nSubs
        test_stat(i_sub) = fn_handle(dv_orig{i_sensor}(i_sub,:)', iv_orig');
    end
    
    [h, p_val, ci, stats] = ttest(test_stat);
    
    topo_p(i_sensor)    = p_val;
    topo_stat(i_sensor) = stats.tstat;    

end

% find clusters in the original data topo
clusters_orig = find_clusters(topo_stat, topo_p, p_thresh);



%% analysis of repeatedly shuffled data set


h = waitbar(0, 'Shuffling...');
for i_rep = 1:nReps

    waitbar(i_rep / nReps);
    
    % construct shuffle permutation for current iteration
    % for each iteration, shuffles are
    % - random and independent across subjects
    % - the same for all sensors within a subject
    shuffle_ind = [];
    for i_sub = 1:nSubs
        shuffle_ind(i_sub, :) = randperm(nCond);
    end


    % 
    for i_sensor = 1:nSensors

        % apply the shuffling to the original data
        for i_sub = 1:nSubs
            ind = shuffle_ind(i_sub, :);
            dv_shuffle{i_sensor}(i_sub, :) = dv_orig{i_sensor}(i_sub, ind);
        end
        
        
        for i_sub = 1:nSubs
            test_stat(i_sub) = fn_handle(dv_shuffle{i_sensor}(i_sub,:)', iv_orig');
        end

        [h, p_val, ci, stats] = ttest(test_stat);

        topo_p_shuffle(i_sensor)    = p_val;
        topo_stat_shuffle(i_sensor) = stats.tstat;    
        
        
    end

    % find clusters in the shuffled data topo and get the max cluster stat
    clusters_shuffle{i_rep} = find_clusters(topo_stat_shuffle, topo_p_shuffle, p_thresh);
    shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs;

end

% delete(h);


%% calculate p-values for each cluster in the original data set 
%  on the basis of the estimated null distribution from the shuffled data set

for i = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat > abs( clusters_orig.cluster_statSum(i) )) / nReps;
    clusters_orig.cluster_pval(i) = pval;
end