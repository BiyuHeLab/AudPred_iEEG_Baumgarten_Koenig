function clusters = find_temporal_clusters(stat_timecourse, p_timecourse, p_thresh)
% clusters = find_clusters(stat_timecourse, p_timecourse, p_thresh)
%
% Find the clusters in a timecourse which surpass some statistical threshold.
%
% inputs
% ------
% stat_timecourse - array of length nSamples which holds the statistic of interest (e.g. t, F) at each sample
% p_timecourse    - array of length nSamples which holds the p-value corresponding to stat_timecourse
% p_thresh        - p-value below which sensors will be considered for cluster analysis
%
% output
% ------
% clusters.inputs             - holds the inputs used in the find_temporal_clusters function call
% clusters.cluster_timecourse - array of length nSamples which holds cluster ID at each sample
%                               (0 for sensors not belonging to any cluster)
% clusters.nClusters          - number of clusters found
% clusters.cluster_samples    - cell array of length nClusters where cell i holds the indeces for
%                               samples belonging to cluster i
% clusters.cluster_size       - array of length nClusters holding the number of samples in each cluster
% clusters.cluster_statSum    - array of length nClusters holding the sum of stat_timecourse for each cluster
% clusters.maxSize            - maximum cluster size
% clusters.maxStatSumPos      - maximum statSum among all clusters with a positive statSum
% clusters.maxStatSumNeg      - minimum statSum among all clusters with a negative statSum
% clusters.maxStatSumAbs      - maximum absolute value of statSum


%% find clusters

sig_timecourse = p_timecourse < p_thresh;
sig_timecourse = sig_timecourse .* sign(stat_timecourse);

%TJB 05.05.21 - Tweak to select direction of effects (-1 or 1), since
%cotiguous.m doesn't work with fixed [-1, 1] input if there no -1 direction in sig_timecourse
effectdirections = unique(sig_timecourse);
effectdirections(find(effectdirections == 0)) = [];

if ~isempty(effectdirections)
    runs = contiguous(sig_timecourse, effectdirections);
    nNegClusters = 0;
    nPosClusters = 0;
    for i_effectdirection = 1:size(runs,1)
        if runs{i_effectdirection,1} == -1
            neg_run_ind = runs{i_effectdirection,2};
            nNegClusters = size(neg_run_ind, 1);
        elseif runs{i_effectdirection,1} == 1    
            pos_run_ind = runs{i_effectdirection,2};
            nPosClusters = size(pos_run_ind, 1);
        end
    end
    
    nClusters = nNegClusters + nPosClusters;
    
    
    %% compile cluster stats
    
    clusterID          = 0;
    cluster_samples    = [];
    cluster_timecourse = zeros(size(stat_timecourse));
    cluster_size       = [];
    cluster_statSum    = [];
    
    % compute cluster statistics for negative clusters
    for i = 1:nNegClusters
        clusterID = clusterID + 1;
        
        cluster_samples{clusterID} = [neg_run_ind(i,1) : neg_run_ind(i,2)];
        
        cluster_timecourse(cluster_samples{clusterID}) = clusterID;
        
        cluster_size(clusterID) = length(cluster_samples{clusterID});
        
        cluster_statSum(clusterID) = sum(stat_timecourse(cluster_samples{clusterID}));
    end
    
    
    % compute cluster statistics for positive clusters
    for i = 1:nPosClusters
        clusterID = clusterID + 1;
        
        cluster_samples{clusterID} = [pos_run_ind(i,1) : pos_run_ind(i,2)];
        
        cluster_timecourse(cluster_samples{clusterID}) = clusterID;
        
        cluster_size(clusterID) = length(cluster_samples{clusterID});
        
        cluster_statSum(clusterID) = sum(stat_timecourse(cluster_samples{clusterID}));
    end
    
    
    %% find max cluster info
    
    maxSize       = -Inf;
    maxStatSumPos = -Inf;
    maxStatSumNeg =  Inf;
    
    for i = 1:nClusters
        
        % update max size
        if cluster_size(i) > maxSize
            maxSize = cluster_size(i);
        end
        
        
        % update max stat sum
        if cluster_statSum(i) >= 0
            if cluster_statSum(i) > maxStatSumPos
                maxStatSumPos = cluster_statSum(i);
            end
            
        else
            if cluster_statSum(i) < maxStatSumNeg
                maxStatSumNeg = cluster_statSum(i);
            end
        end
        
    end
    
    % save info on maximum cluster stats
    if maxStatSumPos == -Inf && maxStatSumNeg == Inf
        maxSize       = 0;
        maxStatSumPos = 0;
        maxStatSumNeg = 0;
        maxStatSumAbs = 0;
        
    elseif maxStatSumPos == -Inf
        maxStatSumPos = 0;
        maxStatSumAbs = abs(maxStatSumNeg);
        
    elseif maxStatSumNeg == Inf
        maxStatSumNeg = 0;
        maxStatSumAbs = maxStatSumPos;
        
    else
        maxStatSumAbs = max([ maxStatSumPos, abs(maxStatSumNeg) ]);
    end
    
    
    
    %% package output
    
    clusters.inputs.stat_timecourse = stat_timecourse;
    clusters.inputs.p_timecourse    = p_timecourse;
    clusters.inputs.p_thresh        = p_thresh;
    
    clusters.cluster_timecourse = cluster_timecourse;
    clusters.nClusters          = nClusters;
    clusters.cluster_samples    = cluster_samples;
    clusters.cluster_size       = cluster_size;
    clusters.cluster_statSum    = cluster_statSum;
    
    clusters.maxSize       = maxSize;
    clusters.maxStatSumPos = maxStatSumPos;
    clusters.maxStatSumNeg = maxStatSumNeg;
    clusters.maxStatSumAbs = maxStatSumAbs;
    
else %if no cluster found, manually set entries to 0
    clusters.inputs.stat_timecourse = stat_timecourse;
    clusters.inputs.p_timecourse    = p_timecourse;
    clusters.inputs.p_thresh        = p_thresh;
    
    clusters.cluster_timecourse = zeros(1,length(stat_timecourse));
    clusters.nClusters          = 0;
    clusters.cluster_samples    = 0;
    clusters.cluster_size       = 0;
    clusters.cluster_statSum    = 0;
    
    clusters.maxSize       = 0;
    clusters.maxStatSumPos = 0;
    clusters.maxStatSumNeg = 0;
    clusters.maxStatSumAbs = 0;    
end

end