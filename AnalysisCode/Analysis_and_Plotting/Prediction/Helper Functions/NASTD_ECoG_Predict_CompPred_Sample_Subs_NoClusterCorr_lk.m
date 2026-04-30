function NASTD_ECoG_Predict_CompPred_Sample_Subs_NoClusterCorr_lk...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

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


%% 4) Compute linear regression between sample-wise neural signal and prediction / prediction error
% NumEffects.Elecswithcluster.PredEffect = 0;
% NumEffects.Elecswithcluster.SimplePredErrEffect = 0;
% NumEffects.Elecswithcluster.ComplexPredErrEffect = 0;

NumEffects.ElecsWithFDRSig.PredEffect = 0;
NumEffects.ElecsWithFDRSig.SimplePredErrEffect = 0;
NumEffects.ElecsWithFDRSig.ComplexPredErrEffect = 0;

%4.1 Compute effects for experimental data
for i_elec = 1:nSensors
    for i_samples = 1:size(Data_p33,2)
        
        %4.1.1 Compute lin reg stats for experimental data for each sample
        %Prediction Effect
        dv_ERF = squeeze(Data_p33(i_elec,i_samples,:));  %p33 ERF at current sensor and sample for all trials
        iv_predp34 = LogFreq_predp34_discretized(i_k,:)'; %predicted p34 tone pitch for all trials
        % linear regression between ERF (per channel, time window) and predicted tone pitch
        stats = regstats(dv_ERF, iv_predp34, 'linear', 'tstat');
        PredEffect.stats{i_elec}{i_samples}  = stats;
        PredEffect.tval{i_elec}(i_samples)   = stats.tstat.t(2);
        PredEffect.pval{i_elec}(i_samples)   = stats.tstat.pval(2);
        PredEffect.beta{i_elec}(i_samples)   = stats.tstat.beta(2);
        stats = [];
        
        %Simple Prediction Error Effect
        dv_ERF = squeeze(Data_p34(i_elec,i_samples,:));  %p34 ERF at current sensor and sample for all trials
        iv_simplepred_error = abs( LogFreq_SimplePredError(i_k,:)' );
        iv_p34 = LogFreq_p34';
        % linear regression with 2 predictors (simple prediction error + p34)
        stats = regstats(dv_ERF, [iv_simplepred_error, iv_p34], 'linear', 'tstat');
        SimplePredErrEffect.stats{i_elec}{i_samples} = stats;
        SimplePredErrEffect.tval{i_elec}(i_samples)  = stats.tstat.t(2);
        SimplePredErrEffect.pval{i_elec}(i_samples)  = stats.tstat.pval(2);
        SimplePredErrEffect.beta{i_elec}(i_samples)  = stats.tstat.beta(2);
        stats = [];
        
        %Complex Prediction Error Effect
        dv_ERF = squeeze(Data_p34(i_elec,i_samples,:));  %p34 ERF at current sensor and sample for all trials
        iv_pred_error = abs( LogFreq_ComplexPredError(i_k,:)' );
        iv_p34 = LogFreq_p34'; %presented final tone (p34)
        % linear regression with 2 predictors (complex prediction error + p34)
        stats = regstats(dv_ERF, [iv_pred_error, iv_p34], 'linear', 'tstat');
        ComplexPredErrEffect.stats{i_elec}{i_samples}    = stats;
        ComplexPredErrEffect.tval{i_elec}(i_samples)     = stats.tstat.t(2);
        ComplexPredErrEffect.pval{i_elec}(i_samples)     = stats.tstat.pval(2);
        ComplexPredErrEffect.beta{i_elec}(i_samples)     = stats.tstat.beta(2);
        stats = [];
    end
    
    %4.1.2 Determine temporally adjacent samples showing a sign. effect (q < 0.05 after FDR) 
    %Prediction Effect
    stat_timecourse = PredEffect.tval{i_elec}; %Statistics over samples
    p_timecourse = PredEffect.pval{i_elec}; %p-values over samples
    [~,~,~,PredEffect.pval_FDR{i_elec}] = fdr_bh(p_timecourse);
    if sum(PredEffect.pval_FDR{i_elec} < 0.05) > 0
        NumEffects.ElecsWithFDRSig.PredEffect = NumEffects.ElecsWithFDRSig.PredEffect + 1;
    end
    
    %Simple Prediction Error Effect
    stat_timecourse = SimplePredErrEffect.tval{i_elec};
    p_timecourse = SimplePredErrEffect.pval{i_elec};
    [~,~,~,SimplePredEffEffect.pval_FDR{i_elec}] = fdr_bh(p_timecourse);
    if sum(SimplePredEffEffect.pval_FDR{i_elec} < 0.05) > 0
        NumEffects.ElecsWithFDRSig.SimplePredErrEffect = NumEffects.ElecsWithFDRSig.SimplePredErrEffect + 1;
    end
    
    %Complex Prediction Error Effect
    stat_timecourse = ComplexPredErrEffect.tval{i_elec};
    p_timecourse = ComplexPredErrEffect.pval{i_elec};
    [~,~,~,ComplexPredErrEffect.pval_FDR{i_elec}] = fdr_bh(p_timecourse);
    if sum(ComplexPredErrEffect.pval_FDR{i_elec} < 0.05) > 0
        NumEffects.ElecsWithFDRSig.ComplexPredErrEffect = NumEffects.ElecsWithFDRSig.ComplexPredErrEffect + 1;
    end
end

disp(['For ' FuncInput_DataType ' & ' FuncInput_ToneDur_text 'TD, ']);
disp([num2str(NumEffects.ElecsWithFDRSig.PredEffect) '/' num2str(nSensors) ...
    ' electrodes show any Prediction-effect cluster-output at cluster-alpha = ' ...
    num2str(param.clusteralpha)]);
disp([num2str(NumEffects.ElecsWithFDRSig.SimplePredErrEffect) '/' num2str(nSensors) ...
    ' electrodes show any SimplePredError-effect cluster-output at cluster-alpha = ' ...
    num2str(param.clusteralpha)]);
disp([num2str(NumEffects.ElecsWithFDRSig.ComplexPredErrEffect) '/' num2str(nSensors) ...
    ' electrodes show any ComplexPredError-effect cluster-output at cluster-alpha = ' ...
    num2str(param.clusteralpha)]);

%% 6) Save variables
param_LoadedData = param;
labels_loadedData = FuncInput_InputData.label;

%Remove cluster_shuff subfield (too big & we only need maxtstats)
% PredEffect = rmfield(PredEffect,'cluster_shuff');
% SimplePredErrEffect = rmfield(SimplePredErrEffect,'cluster_shuff');
% ComplexPredErrEffect = rmfield(ComplexPredErrEffect,'cluster_shuff');

savefile = [path_dataoutput sub '_PredEffectsCluster_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD_NoClusterCorr.mat'];
save(savefile, 'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect', 'NumEffects', ...
    'param_LoadedData', 'labels_loadedData', ...
    'Data_p33', 'Data_p34', ...
    'LogFreq*');

end