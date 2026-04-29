function NASTD_ECoG_Predict_CompPred_Sample_Subs_FTPL...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG, FTPLratings)

%Aim: Compute prediction effect (How is neural activity at tone 33
%modulated by the predicted pitch of tone 34 (predp34)?) and
%prediction error effect (How is neural activity at tone 34
%modulated by the difference between predicted (predp34) and presented
%(p34) pitch of tone 34?)

%Method: Linear regression of pitch values p1-32 on sample-wise neural activity during p33.
%Cluster-correction of lin reg. t-values across samples to get
%sign. prediction effect sample-clusters within p33.

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_dataoutput = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

path_fig_cluster = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ClusterStat/' FuncInput_DataType '/']);
if (~exist(path_fig_cluster, 'dir')); mkdir(path_fig_cluster); end
path_fig_ERF = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ERFs/PredEffect/' FuncInput_DataType '/']);
if (~exist(path_fig_ERF, 'dir')); mkdir(path_fig_ERF); end

%% 1. Select current input data in FT-struct
FuncInput_InputData.trial = [];
FuncInput_InputData.trial = FuncInput_InputData.(FuncInput_DataType);

%% 2) Extract p33 and p34 per trial
%2.1 Define parameters depending on TD
nTrials             = length(FuncInput_InputData.trial);
nSensors            = size(FuncInput_InputData.trial{1},1);
SampleFreq          = FuncInput_InputData.fsample;
ToneDur_Sec         = str2num(FuncInput_ToneDur_text);
nSamples_perTone    = ToneDur_Sec * SampleFreq;

%2.2 Determine TP/samples for each tone start+end
TP_Tone_StartStop = NaN(36,2);
TP_Tone_StartStop(1,1) = 0; %p1 set as t = 0 in trial definition
for i_tone = 2:36
    Dist = abs(FuncInput_InputData.time{1} - ...
        ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
    minDist = min(Dist);
    i_minDist = find(Dist == minDist);
    TP_Tone_StartStop(i_tone,1) = FuncInput_InputData.time{1}(i_minDist);
end
for i_tone = 1:35
    i_LastSampleTone = find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone + 1,1));
    TP_Tone_StartStop(i_tone,2) = FuncInput_InputData.time{1}(i_LastSampleTone);
end
TP_Tone_StartStop = TP_Tone_StartStop(1:34,:);
%Check if all tones are of equal length
%(if not then choose min length by deleting the last sample of longer trials)
minSeqLength_sec = min(TP_Tone_StartStop(:,2)-TP_Tone_StartStop(:,1));
for i_tone = 1:34
    if (TP_Tone_StartStop(i_tone,2) - TP_Tone_StartStop(i_tone,1)) > minSeqLength_sec
        TP_Tone_StartStop(i_tone,2) = ...
            FuncInput_InputData.time{1}(...
            find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,2))-1);
    end
end
%Determine samples corresponding to TP
Sample_Tone_StartStop = NaN(34,2);
for i_tone = 1:34
    Sample_Tone_StartStop(i_tone,1) = ...
        find(FuncInput_InputData.time{1} == TP_Tone_StartStop(i_tone,1));
    Sample_Tone_StartStop(i_tone,2) = ...
        find(TP_Tone_StartStop(i_tone,2) == FuncInput_InputData.time{1});
end

%2.3) Initialize and fill 3D data arrays (nSens, nSamples per tone, nTrials)
Data_p33 = zeros(nSensors, ...
    length(Sample_Tone_StartStop(1,1):Sample_Tone_StartStop(1,2)), ...
    nTrials);
Data_p34 = Data_p33;
%copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
for i_trial = 1:nTrials
    Data_p33(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(param.ToneIndex,1) : Sample_Tone_StartStop(param.ToneIndex,2));
    Data_p34(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(param.ToneIndex+1,1) : Sample_Tone_StartStop(param.ToneIndex+1,2));
    Data_p1to33(:, :, i_trial) = ...
        FuncInput_InputData.trial{i_trial}(:, ...
        Sample_Tone_StartStop(1,1) : Sample_Tone_StartStop(param.ToneIndex,2));
end

%% 3) Extract Predp34 and Prediction Error
%For each trial, compute/read out
%1) Predp34 (discretized to next full tone)
%2) p34 (i.e., the really presented last tone)
%3) Simple Prediction error (i.e., the distance between the final (p34) and penultimate (p33( tone
%4) Complex Prediction error (i.e., the distance between presented (p34) and predicted (Predp34) tone

for i_k = 1:length(param.PredSeqrange) %for each tone for which predition is analyzed
    
    series_start_ind = param.ToneIndex - param.PredSeqrange(i_k) + 1; %beginning of tone sequence part used for prediction
    series_end_ind   = param.ToneIndex; %end of tone sequence part used for prediction
    
    for i_trial = 1:nTrials
        series = FuncInput_InputData.behav.stim.series_f{i_trial}...
            (series_start_ind : series_end_ind); %(non-log) selected sequence tones in Hz
        beta = FuncInput_InputData.behav.stim.beta(i_trial); %always 1.5
        
        LogFreq_p33(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex) ); %log(p33)
        LogFreq_p34(i_trial) = log( FuncInput_InputData.behav.stim.series_f{i_trial}(param.ToneIndex+1) ); %log(p34)
        LogFreq_predp34_discretized(i_trial) = FuncInput_InputData.behav.stim.logf_pred(i_trial); %log(Predp34) discretized
        
        LogFreq_ComplexPredError(i_k, i_trial) = ...
            LogFreq_p34(i_trial) - LogFreq_predp34_discretized(i_trial);
        LogFreq_SimplePredError(i_k, i_trial) = ...
            LogFreq_p34(i_trial) - LogFreq_p33(i_trial);
    end    
end

% Extract FTPLratings for each trial, making sure to align the trial IDs
trialID_p34 = FuncInput_InputData.behav.trialNum';
trialID_FTPL = table2array(FTPLratings(:,3));
[~, idx_match] = ismember(trialID_FTPL, trialID_p34);
FTPL_pertrial = FTPLratings.FTPLrating; %Extract the FTPLratings for each trial
FTPL_pertrial = FTPL_pertrial(idx_match(idx_match > 0));


%% 4) Compute linear regression between sample-wise neural signal and prediction / prediction error
NumEffects.Elecswithcluster.PredErrFTPL = 0;

%4.1 Compute effects for experimental data
for i_elec = 1:nSensors
    for i_samples = 1:size(Data_p33,2)
        %Prediction Error (FTPL) Effect modulating p34
        dv_ERF = squeeze(Data_p34(i_elec,i_samples,:));  %p34 ERF at current sensor and sample for all trials
        iv_FTPL = FTPL_pertrial; %FTPL ratings for all trials    
        % linear regression
        stats = regstats(dv_ERF, iv_FTPL, 'linear', 'tstat');
        PredErrFTPL.stats{i_elec}{i_samples} = stats;
        PredErrFTPL.tval{i_elec}(i_samples)  = stats.tstat.t(2);
        PredErrFTPL.pval{i_elec}(i_samples)  = stats.tstat.pval(2);
        PredErrFTPL.beta{i_elec}(i_samples)  = stats.tstat.beta(2);
        stats = [];
    end
    
    %4.1.2 Determine temporally adjacent samples showing a sign. effect (q < 0.05 after FDR) 
    %Prediction Error (FTPL) Effect modulating p34
    stat_timecourse = PredErrFTPL.tval{i_elec}; %Statistics over samples
    p_timecourse = PredErrFTPL.pval{i_elec}; %p-values over samples
    [~,~,~,PredErrFTPL.pval_FDR{i_elec}] = fdr_bh(p_timecourse);
    PredErrFTPL.clusterstat{i_elec} = ...
        find_temporal_clusters(stat_timecourse, p_timecourse, param.clusteralpha);
    if PredErrFTPL.clusterstat{i_elec}.nClusters > 0
        NumEffects.Elecswithcluster.PredErrFTPL = NumEffects.Elecswithcluster.PredErrFTPL + 1;
    end
end

disp(['For ' FuncInput_DataType ' & ' FuncInput_ToneDur_text 'TD, ']);
disp([num2str(NumEffects.Elecswithcluster.PredErrFTPL) '/' num2str(nSensors) ...
    ' electrodes show any Prediction-effect cluster-output at cluster-alpha = ' ...
    num2str(param.clusteralpha)]);

NumEffects.NumCluster.PredErrFTPL = 0;

for i_elec = 1:nSensors
    NumEffects.NumCluster.PredErrFTPL = NumEffects.NumCluster.PredErrFTPL + PredErrFTPL.clusterstat{i_elec}.nClusters;
end

%4.2. Compute effects for shuffled null distribution
h = waitbar(0, 'Creating shuffled null distribution...');

for i_rep = 1:param.numreps
    
    waitbar(i_rep / param.numreps);
    
    %4.2.1 Shuffle neural data across trials to destroy trial-wise relationship
    %between neural data and tone pitch effect
    %Construct shuffle permutation for current iteration
    ind_shuff = [];
    ind_shuff = randperm(nTrials);
    
    %4.2.2 Apply shuffling to the original data and compute lin reg stats for each sample
    for i_elec = 1:nSensors
        for i_samples = 1:size(Data_p33,2)
            %%Prediction Error (FTPL) Effect modulating p34
            dv_ERF_shuff = squeeze(Data_p34(i_elec,i_samples,ind_shuff));  %p34 ERF at current sensor and sample for shuffled trials
            iv_FTPL = FTPL_pertrial; %FTPL ratings for all trials    
            % linear regression
            stats = regstats(dv_ERF_shuff, iv_FTPL, 'linear', 'tstat');
            PredErrFTPL_Shuff.tval{i_elec}(i_samples)  = stats.tstat.t(2);
            PredErrFTPL_Shuff.pval{i_elec}(i_samples)  = stats.tstat.pval(2);
            stats = [];
        end
        
        %4.2.3 Determine temporally adjacent samples showing a sign. effects and
        %group them into clusters, sum up test statistic for each cluster
        %Prediction Error (FTPL) Effect
        stat_timecourse_Shuff = PredErrFTPL_Shuff.tval{i_elec}; %Statistics over samples
        p_timecourse_Shuff = PredErrFTPL_Shuff.pval{i_elec}; %p-values over samples
        PredErrFTPL.cluster_shuff{i_elec}{i_rep} = ...
            find_temporal_clusters(stat_timecourse_Shuff, p_timecourse_Shuff, param.clusteralpha);
        PredErrFTPL.ClusterMaxStat_shuff{i_elec}(i_rep) = ...
            PredErrFTPL.cluster_shuff{i_elec}{i_rep}.maxStatSumAbs;
    end
end
close(h)

%4.3 Summary figure showing if there is a sign. cluster for the experimental data &
%how many repetitions resulted in a cluster-ouput for shuffled data
% if plot_poststepFigs == 1
%     
%     ExpDistr_ClusterInfo_PredEffect = zeros(1,nSensors);
%     ExpDistr_ClusterInfo_SimplePredErrEffect = zeros(1,nSensors);
%     ExpDistr_ClusterInfo_ComplexPredErrEffect = zeros(1,nSensors);
%     
%     ShuffDistr_ClusterInfo_PredEffect = zeros(1,nSensors);
%     ShuffDistr_ClusterInfo_SimplePredErrEffect = zeros(1,nSensors);
%     ShuffDistr_ClusterInfo_ComplexPredErrEffect = zeros(1,nSensors);
%     
%     for i_elec = 1:nSensors
%         ExpDistr_ClusterInfo_PredEffect(i_elec) = ...
%             PredEffect.clusterstat{i_elec}.nClusters;
%         ExpDistr_ClusterInfo_SimplePredErrEffect(i_elec) = ...
%             SimplePredErrEffect.clusterstat{i_elec}.nClusters;
%         ExpDistr_ClusterInfo_ComplexPredErrEffect(i_elec) = ...
%             ComplexPredErrEffect.clusterstat{i_elec}.nClusters;
%         
%         ShuffDistr_ClusterInfo_PredEffect(i_elec) = ...
%             sum(PredEffect.ClusterMaxStat_shuff{i_elec} > 0) / ...
%             length(PredEffect.ClusterMaxStat_shuff{i_elec});
%         ShuffDistr_ClusterInfo_SimplePredErrEffect(i_elec) = ...
%             sum(SimplePredErrEffect.ClusterMaxStat_shuff{i_elec} > 0) / ...
%             length(SimplePredErrEffect.ClusterMaxStat_shuff{i_elec});
%         ShuffDistr_ClusterInfo_ComplexPredErrEffect(i_elec) = ...
%             sum(ComplexPredErrEffect.ClusterMaxStat_shuff{i_elec} > 0) / ...
%             length(ComplexPredErrEffect.ClusterMaxStat_shuff{i_elec});
%     end
%     
%     
%     %To do: add plot for pred error effects
%     h = figure;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     subplot(3,1,1)
%     bar(1:nSensors,ShuffDistr_ClusterInfo_PredEffect,'k')
%     ylim([0 1])
%     xlabel('electrode number');
%     ylabel('reps clusterStatSum > 0/all reps');
%     hold on;
%     bar(find(ExpDistr_ClusterInfo_PredEffect),...
%         ShuffDistr_ClusterInfo_PredEffect(ExpDistr_ClusterInfo_PredEffect > 0),'r')
%     legendtext1 = ['Electrodes without sign. cluster in exp data'];
%     legendtext2 = ['Electrodes with sign. cluster in exp data'];
%     legend(legendtext1, legendtext2, 'Location', 'Best')
%     title('Prediction Effect')
%     subplot(3,1,2)
%     bar(1:nSensors,ShuffDistr_ClusterInfo_SimplePredErrEffect,'k')
%     ylim([0 1])
%     xlabel('electrode number');
%     ylabel('reps clusterStatSum > 0/all reps');
%     hold on;
%     bar(find(ExpDistr_ClusterInfo_SimplePredErrEffect),...
%         ShuffDistr_ClusterInfo_SimplePredErrEffect(ExpDistr_ClusterInfo_SimplePredErrEffect > 0),'r')
%     title('Simple Prediction Error Effect')
%     subplot(3,1,3)
%     bar(1:nSensors,ShuffDistr_ClusterInfo_ComplexPredErrEffect,'k')
%     ylim([0 1])
%     ylabel('reps clusterStatSum > 0/all reps');
%     hold on;
%     bar(find(ExpDistr_ClusterInfo_ComplexPredErrEffect),...
%         ShuffDistr_ClusterInfo_ComplexPredErrEffect(ExpDistr_ClusterInfo_ComplexPredErrEffect > 0),'r')
%     title('Complex Prediction Error Effect')
%     sgtitle({[sub ' - ' FuncInput_ToneDur_text  's TD - ' FuncInput_DataType], ...
%         ['Overview Cluster-results (cluster-alpha = ' num2str(param.clusteralpha) ') for MCC across time, per electrode']})
%     
%     if save_poststepFigs == 1
%         filename     = ['ClusterStatSummary_' sub '_' FuncInput_DataType '_' ...
%             FuncInput_ToneDur_text 'sTD.png'];
%         figfile      = [path_fig_cluster filename];
%         saveas(gcf, [figfile], 'png'); %save png version
%         close
%     end
% end

%4.4 Create (and plot) a histogram of the summed test statistics from shuffled data
%and compare experimental clusterstats against these distributions
for i_elec = 1:nSensors
    maxNumEffects.NumCluster_PredErrFTPL(i_elec) = PredErrFTPL.clusterstat{i_elec}.nClusters;
end
Summary_Clusterpval_PredErrFTPL = nan(max(maxNumEffects.NumCluster_PredErrFTPL),length(maxNumEffects.NumCluster_PredErrFTPL));
Mask_SignCluster_PredErrFTPL = Summary_Clusterpval_PredErrFTPL;

for i_elec = 1:nSensors
    %Prediction error FTPL effect
    if PredErrFTPL.clusterstat{i_elec}.nClusters > 0
        for i_cluster = 1:PredErrFTPL.clusterstat{i_elec}.nClusters
            pval = sum(PredErrFTPL.ClusterMaxStat_shuff{i_elec} > ...
                abs(PredErrFTPL.clusterstat{i_elec}.cluster_statSum(i_cluster)))...
                / param.numreps;
            PredErrFTPL.clusterstat{i_elec}.cluster_pval(i_cluster) = pval;
            Summary_Clusterpval_PredErrFTPL(i_cluster, i_elec) = pval;
        end
    else
        PredErrFTPL.clusterstat{i_elec}.cluster_pval = NaN;       
    end   
end

%Prediction error FTPL effect
Mask_SignCluster_PredErrFTPL(find(Summary_Clusterpval_PredErrFTPL < param.alpha)) = 1;
NumEffects.PredErrFTPL.signcluster = length(find(Summary_Clusterpval_PredErrFTPL < param.alpha));
NumEffects.PredErrFTPL.signclusterlabel = FuncInput_InputData.label(find(nansum(Mask_SignCluster_PredErrFTPL)))';
disp([num2str(NumEffects.PredErrFTPL.signcluster) 'sign. clusters for PredErrFTPL found at alpha = ' num2str(param.alpha) ' in electrodes:']);
disp([FuncInput_InputData.label(find(nansum(Mask_SignCluster_PredErrFTPL)))']);

%Plot summary figure with histograms for all elecs with sign. clusters
if plot_poststepFigs == 1
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    col = 4;
    rows = ceil(NumEfffects.PredEffect.signcluster/4);
    NumSubplot = 0;
    for i_elec = 1:nSensors
        if PredEffect.clusterstat{i_elec}.nClusters > 0
            for i_cluster = 1:PredEffect.clusterstat{i_elec}.nClusters
                if PredEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                    NumSubplot = NumSubplot + 1;
                    subplot(rows,col,NumSubplot)
                    h1 = histogram(PredEffect.ClusterMaxStat_shuff{i_elec},param.numreps/10);
                    hold on; plot([...
                        abs(PredEffect.clusterstat{i_elec}.cluster_statSum(i_cluster)) ...
                        abs(PredEffect.clusterstat{i_elec}.cluster_statSum(i_cluster))], ...
                        [0 max(h1.Values)],'r', 'LineWidth',3)
                    title([FuncInput_InputData.label{i_elec} ...
                        ' (p = ' num2str(PredEffect.clusterstat{i_elec}.cluster_pval(i_cluster)) ')'])
                    xlim([-1 h1.BinLimits(2)+1]);
                    xlabel('SummedStat')
                    ylabel('Frequency')
                    if NumSubplot == 1
                        legend('Shuffled null distribution','Experimental data')
                    end
                    axis tight
                    xlim([-5 h1.BinLimits(2)+5]);
                end
            end
        end
    end
    sgtitle({[sub ' - ' FuncInput_ToneDur_text  's TD - ' FuncInput_DataType], ...
        ['Overview Sign. Cluster for PredEffect (cluster-alpha = ' num2str(param.clusteralpha) ')']})
    if save_poststepFigs == 1
        filename     = ['SignClusterHists_PredEffect_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD.png'];
        figfile      = [path_fig_cluster filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
end

if plot_poststepFigs == 1
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    col = 4;
    rows = ceil(NumEfffects.SimplePredErrEffect.signcluster/4);
    NumSubplot = 0;
    for i_elec = 1:nSensors
        if SimplePredErrEffect.clusterstat{i_elec}.nClusters > 0
            for i_cluster = 1:SimplePredErrEffect.clusterstat{i_elec}.nClusters
                if SimplePredErrEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                    NumSubplot = NumSubplot + 1;
                    subplot(rows,col,NumSubplot)
                    h1 = histogram(SimplePredErrEffect.ClusterMaxStat_shuff{i_elec},param.numreps/10);
                    hold on; plot([...
                        abs(SimplePredErrEffect.clusterstat{i_elec}.cluster_statSum(i_cluster)) ...
                        abs(SimplePredErrEffect.clusterstat{i_elec}.cluster_statSum(i_cluster))], ...
                        [0 max(h1.Values)],'r', 'LineWidth',3)
                    title([FuncInput_InputData.label{i_elec} ...
                        ' (p = ' num2str(SimplePredErrEffect.clusterstat{i_elec}.cluster_pval(i_cluster)) ')'])
                    xlim([-1 h1.BinLimits(2)+1]);
                    xlabel('SummedStat')
                    ylabel('Frequency')
                    if NumSubplot == 1
                        legend('Shuffled null distribution','Experimental data')
                    end
                    axis tight
                    xlim([-5 h1.BinLimits(2)+5]);
                end
            end
        end
    end
    sgtitle({[sub ' - ' FuncInput_ToneDur_text  's TD - ' FuncInput_DataType], ...
        ['Overview Sign. Cluster for SimplePredErrEffect (cluster-alpha = ' num2str(param.clusteralpha) ')']})
    if save_poststepFigs == 1
        filename     = ['SignClusterHists_SimplePredErrEffect_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD.png'];
        figfile      = [path_fig_cluster filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
end

if plot_poststepFigs == 1
    h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    col = 4;
    rows = ceil(NumEfffects.ComplexPredErrEffect.signcluster/4);
    NumSubplot = 0;
    for i_elec = 1:nSensors
        if ComplexPredErrEffect.clusterstat{i_elec}.nClusters > 0
            for i_cluster = 1:ComplexPredErrEffect.clusterstat{i_elec}.nClusters
                if ComplexPredErrEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                    NumSubplot = NumSubplot + 1;
                    subplot(rows,col,NumSubplot)
                    h1 = histogram(ComplexPredErrEffect.ClusterMaxStat_shuff{i_elec},param.numreps/10);
                    hold on; plot([...
                        abs(ComplexPredErrEffect.clusterstat{i_elec}.cluster_statSum(i_cluster)) ...
                        abs(ComplexPredErrEffect.clusterstat{i_elec}.cluster_statSum(i_cluster))], ...
                        [0 max(h1.Values)],'r', 'LineWidth',3)
                    title([FuncInput_InputData.label{i_elec} ...
                        ' (p = ' num2str(ComplexPredErrEffect.clusterstat{i_elec}.cluster_pval(i_cluster)) ')'])
                    xlim([-1 h1.BinLimits(2)+1]);
                    xlabel('SummedStat')
                    ylabel('Frequency')
                    if NumSubplot == 1
                        legend('Shuffled null distribution','Experimental data')
                    end
                    axis tight
                    xlim([-5 h1.BinLimits(2)+5]);
                end
            end
        end
    end
    sgtitle({[sub ' - ' FuncInput_ToneDur_text  's TD - ' FuncInput_DataType], ...
        ['Overview Sign. Cluster for ComplexPredErrEffect (cluster-alpha = ' num2str(param.clusteralpha) ')']})
    if save_poststepFigs == 1
        filename     = ['SignClusterHists_ComplexPredErrEffect_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD.png'];
        figfile      = [path_fig_cluster filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
end

% %% 5. Plotting ERFs per predp34 with sign. clusters
% %Plot timelocked activity during p33 per predp34 and highlight sign. clusters
% 
% %Separate trials by Predp34
% label_Predp34 = [-1 0 1];% low, medium, high;
% for i_Predp34 = 1:length(label_Predp34)
%     TrialFilt_Predp34 = FuncInput_InputData.behav.stim.predID == label_Predp34(i_Predp34);
%     Index_Predp34{i_Predp34} = find(TrialFilt_Predp34 == 1);
%     %Determine trial number in case we want to restrict analysis to equal
%     %amount of trials across predp34conditions
%     NumTrialsperpredp34(i_Predp34) = length(Index_Predp34{i_Predp34});
% end
% MinNumTrialsperpredp34 = min(NumTrialsperpredp34);
% 
% if plot_poststepFigs == 1
%     
%     %Compute AVG + STD across trials per Predp34
%     for i_Predp34 = 1:length(label_Predp34)
%         Avgp33{i_Predp34} = mean(Data_p33(:,:,Index_Predp34{i_Predp34}),3);
%         STDp33{i_Predp34} = std(Data_p33(:,:,Index_Predp34{i_Predp34}),1,3);
%     end
%     
%     %Plot ERF for all selected elecs
%     NumFig = ceil(nSensors/30); %20 elecs per plot
%     plot_param.color = ...
%         {[0, 0.4470, 0.7410, 0.7], ...
%         [0, 0.75, 0.75, 0.7],[0.8500, ...
%         0.3250, 0.0980, 0.7]};
%     
%     %ERF for p33
%     i_elec = 0;
%     for i_fig = 1:NumFig
%         h = figure;
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         for i_subplot = 1:30
%             if i_elec < nSensors
%                 subplot(5,6,i_subplot)
%                 i_elec = i_elec+1;
%                 
%                 %Plot signal traces
%                 plot(1:length(Avgp33{1}(i_elec,:)),Avgp33{1}(i_elec,:),...
%                     'color',plot_param.color{1},'LineWidth', 2)
%                 hold on;
%                 plot(1:length(Avgp33{2}(i_elec,:)),Avgp33{2}(i_elec,:),...
%                     'color',plot_param.color{2},'LineWidth', 2)
%                 hold on;
%                 plot(1:length(Avgp33{3}(i_elec,:)),Avgp33{3}(i_elec,:),...
%                     'color',plot_param.color{3},'LineWidth', 2)
%                 if isnan(PredEffect.clusterstat{i_elec}.cluster_pval)
%                     title(['Elec: ' FuncInput_InputData.label{i_elec}])
%                 elseif isempty(PredEffect.clusterstat{i_elec}.cluster_pval)
%                     title(['Elec: ' FuncInput_InputData.label{i_elec}])
%                 elseif any(PredEffect.clusterstat{i_elec}.cluster_pval < param.alpha)
%                     temp_titletext = [];
%                     for i_cluster = 1:length(PredEffect.clusterstat{i_elec}.cluster_pval)
%                         if PredEffect.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
%                             temp_titletext = ...
%                                 [temp_titletext ...
%                                 num2str(round(...
%                                 PredEffect.clusterstat{i_elec}.cluster_pval(i_cluster),3))...
%                                 ' '];
%                         end
%                     end
%                     title(['Elec: ' FuncInput_InputData.label{i_elec} ...
%                         ' (p = ' temp_titletext ')'])
%                 end
%                 axis tight
%                 xlabel('samples (p33)')
%                 
%                 %Highlight sign. epochs
%                 maxY = max([max(Avgp33{1}(i_elec,:)), max(Avgp33{2}(i_elec,:)), max(Avgp33{3}(i_elec,:))]);
%                 minY = min([min(Avgp33{1}(i_elec,:)), min(Avgp33{2}(i_elec,:)), min(Avgp33{3}(i_elec,:))]);
%                 grey  = [100 100 100]./255;
%                 if any(PredEffect.clusterstat{i_elec}.cluster_pval < param.alpha)
%                     area(1:length(PredEffect.clusterstat{i_elec}.cluster_timecourse), ...
%                         [PredEffect.clusterstat{i_elec}.cluster_timecourse ./ ...
%                         PredEffect.clusterstat{i_elec}.cluster_timecourse * maxY],...
%                         'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
%                     hold on;
%                     area(1:length(PredEffect.clusterstat{i_elec}.cluster_timecourse), ...
%                         [PredEffect.clusterstat{i_elec}.cluster_timecourse ./ ...
%                         PredEffect.clusterstat{i_elec}.cluster_timecourse * minY],...
%                         'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');                    
%                     ylim([minY maxY]);
%                 end
%             end
%         end
%         
%         %Add subtitle
%         sgtitle([sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
%             'msTD - ERF p33 as a function of p*34 (' ...
%             num2str(i_fig) '/' num2str(NumFig) ')'],'Interpreter','none')
%         
%         if save_poststepFigs == 1
%             filename     = ['p33ERF_AllElecs_ClusterStat_' sub '_' FuncInput_DataType '_' ...
%                 FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
%             figfile      = [path_fig_ERF filename];
%             saveas(gcf, [figfile], 'png'); %save png version
%             close
%         end
%     end
% end

%% 6) Save variables
param_LoadedData = param;
labels_loadedData = FuncInput_InputData.label;

%Remove cluster_shuff subfield (too big & we only need maxtstats)
PredErrFTPL = rmfield(PredErrFTPL,'cluster_shuff');

savefile = [path_dataoutput sub '_PredEffectsFTPL_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'];
save(savefile, 'PredErrFTPL', 'NumEffects', ...
    'param_LoadedData', 'labels_loadedData', ...
    'Data_p33', 'Data_p34', ...
    'LogFreq*','FTPL_pertrial');

end