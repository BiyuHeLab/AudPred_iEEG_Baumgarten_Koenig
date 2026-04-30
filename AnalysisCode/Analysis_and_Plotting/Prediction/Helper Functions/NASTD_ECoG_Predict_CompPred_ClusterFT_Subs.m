function NASTD_ECoG_Predict_CompPred_ClusterFT_Subs...
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
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_dataoutput = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];
if (~exist(path_dataoutput, 'dir')); mkdir(path_dataoutput); end

path_fig_cluster = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ClusterStat/' FuncInput_DataType '/']);
if (~exist(path_fig_cluster, 'dir')); mkdir(path_fig_cluster); end
path_fig_ERF = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ERFsperPredp34/' FuncInput_DataType '/']);
if (~exist(path_fig_ERF, 'dir')); mkdir(path_fig_ERF); end

%% 1. Select current input data in FT-struct
FuncInput_InputData.trial = [];
FuncInput_InputData.trial = FuncInput_InputData.(FuncInput_DataType);

%% 2) Extract p33 and p34 per trial
%2.1 Define parameters depending on TD
nTrials     = length(FuncInput_InputData.trial);
nSensors    = size(FuncInput_InputData.trial{1},1);
SampleFreq  = FuncInput_InputData.fsample;
ToneDur_Sec  = str2num(FuncInput_ToneDur_text);
nSamples_perTone = ToneDur_Sec * SampleFreq;

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

%% 4. FT-based approach incl. cluster-corrction across samples (works only
%for prediction effect since herev only 1 predictor)

%Compute cluster-corrected statistic for each electrode
for i_elec = 1:nSensors
    
    %Option 1: Independent samples regression T-statistic
    %with trials treated as independent unit-of-observation. Tests if there
    %is a relation ship between IV/predictor (predp34) and neural data
    %across time.
    
    %Prepare input structure(timelocked ERF)
    cfg                 = [];
    cfg.channel         = FuncInput_InputData.label{i_elec}; %current electrode
    cfg.latency         = ...
        [FuncInput_InputData.time{1}(Sample_Tone_StartStop(33,1)) ...
        FuncInput_InputData.time{1}(Sample_Tone_StartStop(33,2))]; %p33 time points
    cfg.keeptrials      = 'yes';
    StatInput_TLData    = ft_timelockanalysis(cfg, FuncInput_InputData);
    
    %Compute linear regression across samples
    cfg                     = [];
    cfg.avgovertime         = 'no';
    cfg.parameter           = 'trial';
    
    cfg.statistic           = 'ft_statfun_indepsamplesregrT';
    cfg.method              = 'montecarlo';
    cfg.alpha               = 0.05;
    cfg.tail                = 0;   %0 = two-sided test
    cfg.correcttail         = 'prob';   %Correct prob subfield to reflect two-sided test
    cfg.correctm            = 'cluster';
    cfg.clusteralpha        = 0.05;
    cfg.clusterstatistic    = 'maxsum';
    cfg.clustertail         = 0;
    
    cfg.computestat         = 'yes';
    cfg.computecritval      = 'yes';
    cfg.computeprob         = 'yes';
    
    cfg.numrandomization    = param.numreps;
    
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable/predictor (i.e., tone pitch across trials)
    cfg.design(1,1:nTrials) = LogFreq_predp34_discretized; %predp34 tone pitch per trial
    
    PredEffect.ClusterStats{i_elec} = ...
        ft_timelockstatistics(cfg, StatInput_TLData);
    
    %Add subfield listing cluster-stats independent of pos/neg
    if sum(PredEffect.ClusterStats{i_elec}.mask) > 0
        temp_clusterstat = [];
        if ~isempty(PredEffect.ClusterStats{i_elec}.posclusters)
            for i_cluster = 1:length(PredEffect.ClusterStats{i_elec}.posclusters)
                if PredEffect.ClusterStats{i_elec}.posclusters(i_cluster).prob <= 0.05
                    temp_clusterstat = ...
                        [temp_clusterstat PredEffect.ClusterStats{i_elec}.posclusters(i_cluster).prob];
                end
            end
        end
        if ~isempty(PredEffect.ClusterStats{i_elec}.negclusters)
            for i_cluster = 1:length(PredEffect.ClusterStats{i_elec}.negclusters)
                if PredEffect.ClusterStats{i_elec}.negclusters(i_cluster).prob <= 0.05
                    temp_clusterstat = ...
                        [temp_clusterstat PredEffect.ClusterStats{i_elec}.negclusters(i_cluster).prob];
                end
            end
        end
        PredEffect.ClusterStats{i_elec}.SignCluster_pval = temp_clusterstat;
        clear temp_clusterstat
    else
        PredEffect.ClusterStats{i_elec}.SignCluster_pval = [];
    end
    
    %Cleanup
    clear StatInput_TLData
end

%     %Visually Compare linear reg. stats and pvals between regstat and FT-computation
%     figure;
%     plot(PredEffect.ClusterStats{i_elec}.stat,'b.');
%     hold on;
%     plot(PredEffect.tval{i_elec},'r--');
%     hold on;
%     plot(PredEffect.ClusterStats{i_elec}.prob,'b:');
%     hold on;
%     plot(PredEffect.pval{i_elec},'r-.');
%     legend('t-stat FT-computation','t-stat regstats computation',...
%         'p-val FT-computation','p-val regstats computation');

%Count number of electrodes showing sign. effects
signElecs = [];
for i_elec = 1:nSensors
    if sum(PredEffect.ClusterStats{i_elec}.mask) > 0
        signElecs = [signElecs i_elec];
    end
end
disp([num2str(length(signElecs)) '/' num2str(nSensors) ' elecs show sign. prediction effect (cluster-corrected)'])
FuncInput_InputData.label{signElecs}


%% 5. Plotting
%Plot timelocked activity during p33 per predp34 and highlight sign. clusters

%Separate trials by Predp34
label_Predp34 = [-1 0 1];% low, medium, high;
for i_Predp34 = 1:length(label_Predp34)
    TrialFilt_Predp34 = FuncInput_InputData.behav.stim.predID == label_Predp34(i_Predp34);
    Index_Predp34{i_Predp34} = find(TrialFilt_Predp34 == 1);
    %Determine trial number in case we want to restrict analysis to equal
    %amount of trials across predp34conditions
    NumTrialsperpredp34(i_Predp34) = length(Index_Predp34{i_Predp34});
end
MinNumTrialsperpredp34 = min(NumTrialsperpredp34);

if plot_poststepFigs == 1
    
    %Compute AVG + STD across trials per Predp34
    for i_Predp34 = 1:length(label_Predp34)
        Avgp33{i_Predp34} = mean(Data_p33(:,:,Index_Predp34{i_Predp34}),3);
        STDp33{i_Predp34} = std(Data_p33(:,:,Index_Predp34{i_Predp34}),1,3);
    end
    
    %Plot ERF for all selected elecs
    NumFig = ceil(nSensors/30); %20 elecs per plot
    plot_param.color = ...
        {[0, 0.4470, 0.7410, 0.7], ...
        [0, 0.75, 0.75, 0.7],[0.8500, ...
        0.3250, 0.0980, 0.7]};
    
    %ERF for p33
    i_elec = 0;
    for i_fig = 1:NumFig
        h = figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        for i_subplot = 1:30
            if i_elec < nSensors
                subplot(5,6,i_subplot)
                i_elec = i_elec+1;
                
                %Plot signal traces
                plot(1:length(Avgp33{1}(i_elec,:)),Avgp33{1}(i_elec,:),...
                    'color',plot_param.color{1},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp33{2}(i_elec,:)),Avgp33{2}(i_elec,:),...
                    'color',plot_param.color{2},'LineWidth', 2)
                hold on;
                plot(1:length(Avgp33{3}(i_elec,:)),Avgp33{3}(i_elec,:),...
                    'color',plot_param.color{3},'LineWidth', 2)
                if sum(PredEffect.ClusterStats{i_elec}.mask) == 0
                    title(['Elec: ' FuncInput_InputData.label{i_elec}])
                elseif sum(PredEffect.ClusterStats{i_elec}.mask) > 0
                    title(['Elec: ' FuncInput_InputData.label{i_elec} ...
                        ' (clusterpval: ' num2str(round(PredEffect.ClusterStats{i_elec}.SignCluster_pval(1),3)) ')'])
                end
                axis tight
                xlabel('samples (p33)')
                
                %Highlight sign. epochs
                maxY = max([max(Avgp33{1}(i_elec,:)), max(Avgp33{2}(i_elec,:)), max(Avgp33{3}(i_elec,:))]);
                minY = min([min(Avgp33{1}(i_elec,:)), min(Avgp33{2}(i_elec,:)), min(Avgp33{3}(i_elec,:))]);
                grey  = [100 100 100]./255;
                if sum(PredEffect.ClusterStats{i_elec}.mask) > 0
                    area(1:length(PredEffect.ClusterStats{i_elec}.mask), ...
                        PredEffect.ClusterStats{i_elec}.mask * maxY,...
                        'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
                    hold on;
                    %                     area(1:length(PredEffect.ClusterStats{i_elec}.mask), ...
                    %                         PredEffect.ClusterStats{i_elec}.mask * minY,...
                    %                         'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
                    %                     hold on;
                    ylim([minY maxY]);
                end
            end
        end
        
        %Add subtitle
        sgtitle([sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
            'msTD - ERF p33 as a function of p*34 (' ...
            num2str(i_fig) '/' num2str(NumFig) ')'],'Interpreter','none')
        
        if save_poststepFigs == 1
            filename     = ['ERFp33_AllElecs_ClusterStat_' sub '_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
            figfile      = [path_fig_ERF filename];
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
    end
end


%% 6) Save variables
param_LoadedData = param;
labels_loadedData = FuncInput_InputData.label;

savefile = [path_dataoutput sub '_PredEffectsFT_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'];
save(savefile, 'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect', ...
    'param_LoadedData', 'labels_loadedData', ...
    'Data_p33', 'Data_p34', ...
    'Index_Predp34', ...
    'LogFreq*');

%% Dump
%     %Option 2: Dependent-samples regression T-statistic
%     %with trials treated as dependent unit-of-observation.
%     %Tests if there is a difference in the relation ship between
%     %IV/predictor (predp34) and neural data across time between the predp34 conditions.
%
%     %Separate trials by Predp34 and randomly select same number of trials
%     %across conditions
%     for i_Predp34 = 1:length(label_Predp34)
%         cfg                         = [];
%         cfg.channel                 = FuncInput_InputData.label{i_elec}; %current electrode
%         TrialSelection              = Shuffle(Index_Predp34{i_Predp34});
%         TrialSelection              = TrialSelection(1:MinNumTrialsperpredp34);
%         cfg.trials                  = TrialSelection;
%         cfg.latency                 = ...
%             [FuncInput_InputData.time{1}(Sample_Tone_StartStop(33,1)) ...
%             FuncInput_InputData.time{1}(Sample_Tone_StartStop(33,2))]; %p33 time points
%         cfg.keeptrials              = 'yes';
%         StatInput_TLData{i_Predp34} = ft_timelockanalysis(cfg, FuncInput_InputData);
%     end
%
%     %Compute linear regression across samples
%     cfg                     = [];
%     cfg.avgovertime         = 'no';
%     cfg.parameter           = 'trial';
%
%     cfg.statistic           = 'ft_statfun_depsamplesregrT';
%     cfg.method              = 'montecarlo';
%     cfg.alpha               = 0.05;
%     cfg.tail                = 0;        %0 = two-sided test
%     cfg.correcttail         = 'prob';   %Correct prob subfield to reflect two-sided test
%     cfg.correctm            = 'cluster';
%     cfg.clusteralpha        = 0.05;
%     cfg.clusterstatistic    = 'maxsum';
%     cfg.clustertail         = 0;
%
%     cfg.computestat         = 'yes';
%     cfg.computecritval      = 'yes';
%     cfg.computeprob         = 'yes';
%
%     cfg.numrandomization    = param.NumRepsClusterStat;
%
%     cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable/predictor (i.e., predp34 tone pitch)
%     cfg.uvar                = 2; % the 2nd row in cfg.design contains the units of observation (i.e., trials)
%
%     predp34_pitch           = unique(LogFreq_predp34_discretized);
%     NumTrialsStatComp       = ...
%         [size(StatInput_TLData{1}.trial,1) + ...
%         size(StatInput_TLData{2}.trial,1) + ...
%         size(StatInput_TLData{2}.trial,1)];
%     cfg.design(1,1:NumTrialsStatComp) = ones(length(NumTrialsStatComp));
%     cfg.design(1,1:size(StatInput_TLData{1}.trial,1)) = predp34_pitch(1);
%     cfg.design(1,size(StatInput_TLData{1}.trial,1) + 1:...
%         size(StatInput_TLData{1}.trial,1 )+ ...
%         size(StatInput_TLData{2}.trial,1)) = predp34_pitch(2);
%     cfg.design(1,size(StatInput_TLData{1}.trial,1) + ...
%         size(StatInput_TLData{2}.trial,1)+1:...
%         size(StatInput_TLData{1}.trial,1) + ...
%         size(StatInput_TLData{2}.trial,1) + ...
%         size(StatInput_TLData{3}.trial,1)) = predp34_pitch(3);
%     cfg.design(2,1:NumTrialsStatComp) = ...
%         [1:size(StatInput_TLData{1}.trial,1), ...
%         1:size(StatInput_TLData{2}.trial,1), ...
%         1:size(StatInput_TLData{3}.trial,1)];
%
%     PredEffect.ClusterStats{i_elec} = ...
%         ft_timelockstatistics(cfg, ...
%         StatInput_TLData{1}, StatInput_TLData{2}, StatInput_TLData{3});
end