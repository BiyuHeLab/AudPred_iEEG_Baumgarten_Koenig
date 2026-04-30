function NASTD_ECoG_Predict_PlotERF4SignClusterElec_AllSub...
    (subs, FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot a summmary figure showing sign. electrodes on surface brain for
%all subjects. For each of these electrodes, show also the traces per p*34
%conditionduring p1-p33 and during p33. Group ERFs by anatomical location.

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

if param.FDRcorrect == 1
    FDR_label = 'FDRcorr';
else
    FDR_label = 'uncorr';
end

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_ERF/' FuncInput_EffectType ...
    '/ClusterCorr/' param.ElecSelect '/' FDR_label '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) For all subjects, prepare data underlying prediction effect analysis
% Read out preprocessed data for current effect, signal component, TD
tic
for i_sub = 1:length(subs)
    
    %1.1 Load in preprocessed neural and behavioral data
    sub = subs{i_sub};
    disp(['-- Preparing pre-processed data for subject: ' sub ' --'])
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
    tic
    load(loadfile_ECoGpreprocdata);
    
    %1.2 Select relevant electrodes and trials
    temp_DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);
    temp_Index_Valid_Elecs = setdiff(temp_DataClean_CleanTrials.cfg.info_elec.selected.index4EDF, ...
        subs_PreProcSettings.(sub).rejectedChan_index); %via index
    temp_Label_Valid_Elecs = setdiff(temp_DataClean_CleanTrials.cfg.info_elec.selected.Label, ...
        subs_PreProcSettings.(sub).rejectedChan_label); %via label
    if ~isempty(find((strcmp(...
            sort(temp_DataClean_CleanTrials.label(temp_Index_Valid_Elecs)), ...
            sort(temp_Label_Valid_Elecs))) == 0))
        disp(['Mismatch in electrode selection']);
    end
    
    cfg         = [];
    cfg.channel = (temp_Label_Valid_Elecs);
    temp_DataClean_CleanTrialsElecs = ft_selectdata(cfg, temp_DataClean_CleanTrials);
    
    temp_DataClean_CleanTrialsElecs = ...
        NASTD_ECoG_Preproc_SelTrialsperTD...
        (FuncInput_ToneDur_text, temp_DataClean_CleanTrialsElecs);
    
    disp(['Total trial number for TD ' FuncInput_ToneDur_text 's: ' ...
        num2str(length(temp_DataClean_CleanTrialsElecs.trial))])
    
    %Adjust info subfield
    temp_DataClean_CleanTrialsElecs.info.elec = ...
        temp_DataClean_CleanTrialsElecs.cfg.previous.info_elec;
    temp_DataClean_CleanTrialsElecs.info.trigger = ...
        temp_DataClean_CleanTrialsElecs.cfg.previous.info_trigger;
    temp_DataClean_CleanTrialsElecs.info.ref = ...
        temp_DataClean_CleanTrialsElecs.cfg.previous.info_ref;
    
    temp_DataClean_CleanTrialsElecs.cfg.previous = ...
        rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_elec');
    temp_DataClean_CleanTrialsElecs.cfg.previous = ...
        rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_trigger');
    temp_DataClean_CleanTrialsElecs.cfg.previous = ...
        rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_ref');
    
    %1.3 Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus (i.e., tone sequence start) baseline
    temp_BL_win = [-0.5 0];
    for i_trial = 1:length(temp_DataClean_CleanTrialsElecs.trial)
        for i_elec = 1:length(temp_DataClean_CleanTrialsElecs.label)
            temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,:) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,:) - ...
                nanmean(temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,...
                find(temp_DataClean_CleanTrialsElecs.time{1} == temp_BL_win(1)): ...
                find(temp_DataClean_CleanTrialsElecs.time{1} == temp_BL_win(2))));
        end
    end
    
    %1.4 Remove trials where some electrodes have only NaN samples (NY798)
    temp_DataClean_CleanTrialsElecs = ...
        NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
        (sub, FuncInput_ToneDur_text, temp_DataClean_CleanTrialsElecs);
    
    %1.5 Read out current signal component
    if strcmp(FuncInput_DataType, 'LP35Hz')
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.LP35Hz] = ...
            NASTD_ECoG_FiltNaNinterp_LP35Hz...
            (sub, FuncInput_ToneDur_text, ...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    elseif strcmp(FuncInput_DataType, 'HP01toLP30Hz')
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP01toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP01toLP30Hz...
            (sub, FuncInput_ToneDur_text, ...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
    
    elseif strcmp(FuncInput_DataType, 'HP05toLP30Hz')
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, FuncInput_ToneDur_text, ...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
    
    elseif strcmp(FuncInput_DataType, 'HP1toLP30Hz')
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP1toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP1toLP30Hz...
            (sub, FuncInput_ToneDur_text, ...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
    
    elseif strcmp(FuncInput_DataType, 'HP2toLP30Hz')
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP2toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP2toLP30Hz...
            (sub, FuncInput_ToneDur_text, ...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    elseif strcmp(FuncInput_DataType, 'Alpha_LogAmp')
        temp_FrequencyBands = [8 12];
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.Alpha_LogAmp] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel...
            (sub, FuncInput_ToneDur_text, ...
            temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    elseif strcmp(FuncInput_DataType, 'Beta_LogAmp')
        temp_FrequencyBands = [15 30];
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.Beta_LogAmp] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel...
            (sub, FuncInput_ToneDur_text, ...
            temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    elseif strcmp(FuncInput_DataType, 'Gamma_LogAmp')
        temp_FrequencyBands = [30 70];
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.Gamma_LogAmp] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel...
            (sub, FuncInput_ToneDur_text, ...
            temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    elseif strcmp(FuncInput_DataType, 'HighGamma_LogAmp')
        temp_FrequencyBands = [70 150];
        [temp_DataClean_CleanTrialsElecs.trial, ...
            temp_DataClean_CleanTrialsElecs.info.FiltInterp.HighGamma_LogAmp] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel...
            (sub, FuncInput_ToneDur_text, ...
            temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
            temp_DataClean_CleanTrialsElecs, ...
            0, 0, paths_NASTD_ECoG);
        
    end
    
    %1.6 Extract p1-33, p33, and p34 for each trial
    %Define parameters depending on TD
    temp_nTrials             = length(temp_DataClean_CleanTrialsElecs.trial);
    temp_nSensors            = size(temp_DataClean_CleanTrialsElecs.trial{1},1);
    temp_ToneDur_Sec         = str2num(FuncInput_ToneDur_text);
    
    %Determine TP/samples for each tone start+end
    StimTiming.TP_Tone_StartStop{i_sub} = NaN(36,2);
    StimTiming.TP_Tone_StartStop{i_sub}(1,1) = 0; %p1 set as t = 0 in trial definition
    for i_tone = 2:36
        temp_Dist = abs(temp_DataClean_CleanTrialsElecs.time{1} - ...
            ((str2num(FuncInput_ToneDur_text)*i_tone) - str2num(FuncInput_ToneDur_text)));
        temp_minDist = min(temp_Dist);
        i_minDist = find(temp_Dist == temp_minDist);
        StimTiming.TP_Tone_StartStop{i_sub}(i_tone,1) = temp_DataClean_CleanTrialsElecs.time{1}(i_minDist);
    end
    for i_tone = 1:35
        i_LastSampleTone = find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub}(i_tone + 1,1));
        StimTiming.TP_Tone_StartStop{i_sub}(i_tone,2) = temp_DataClean_CleanTrialsElecs.time{1}(i_LastSampleTone);
    end
    StimTiming.TP_Tone_StartStop{i_sub} = StimTiming.TP_Tone_StartStop{i_sub}(1:34,:);
    
    %Check if all tones are of equal length
    %(if not then choose min length by deleting the last sample of longer trials)
    temp_minSeqLength_sec = min(StimTiming.TP_Tone_StartStop{i_sub}(:,2)-StimTiming.TP_Tone_StartStop{i_sub}(:,1));
    for i_tone = 1:34
        if (StimTiming.TP_Tone_StartStop{i_sub}(i_tone,2) - StimTiming.TP_Tone_StartStop{i_sub}(i_tone,1)) > temp_minSeqLength_sec
            StimTiming.TP_Tone_StartStop{i_sub}(i_tone,2) = ...
                temp_DataClean_CleanTrialsElecs.time{1}(...
                find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub}(i_tone,2))-1);
        end
    end
    
    %Determine samples corresponding to TP
    StimTiming.Sample_Tone_StartStop{i_sub} = NaN(34,2);
    for i_tone = 1:34
        StimTiming.Sample_Tone_StartStop{i_sub}(i_tone,1) = ...
            find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub}(i_tone,1));
        StimTiming.Sample_Tone_StartStop{i_sub}(i_tone,2) = ...
            find(StimTiming.TP_Tone_StartStop{i_sub}(i_tone,2) == temp_DataClean_CleanTrialsElecs.time{1});
    end
    
    %Initialize and fill 3D data arrays (dim: nSens, nSamples per tone,
    %nTrials) and copy electrode lebels so that we find electrode indices
    ERFdata.p33{i_sub} = zeros(temp_nSensors, ...
        length(StimTiming.Sample_Tone_StartStop{i_sub}(1,1):StimTiming.Sample_Tone_StartStop{i_sub}(1,2)), ...
        temp_nTrials);
    ERFdata.p34{i_sub} = ERFdata.p33{i_sub};
    ERFdata.elec_labels{i_sub} = temp_DataClean_CleanTrialsElecs.label;
    
    %copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
    for i_trial = 1:temp_nTrials
        ERFdata.p33{i_sub}(:, :, i_trial) = ...
            temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
            StimTiming.Sample_Tone_StartStop{i_sub}(param.ToneIndex,1) : StimTiming.Sample_Tone_StartStop{i_sub}(param.ToneIndex,2));
        ERFdata.p34{i_sub}(:, :, i_trial) = ...
            temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
            StimTiming.Sample_Tone_StartStop{i_sub}(param.ToneIndex+1,1) : StimTiming.Sample_Tone_StartStop{i_sub}(param.ToneIndex+1,2));
        ERFdata.p1p33{i_sub}(:, :, i_trial) = ...
            temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
            StimTiming.Sample_Tone_StartStop{i_sub}(1,1) : StimTiming.Sample_Tone_StartStop{i_sub}(param.ToneIndex,2));
    end
    
    %1.7 Separate trials by Predp34
    label_Predp34 = [-1 0 1];% low, medium, high;
    for i_Predp34 = 1:length(label_Predp34)
        temp_TrialFilt_Predp34 = temp_DataClean_CleanTrialsElecs.behav.stim.predID == label_Predp34(i_Predp34);
        ERFdata.Index_Predp34{i_sub}{i_Predp34} = find(temp_TrialFilt_Predp34 == 1);
        %Determine trial number in case we want to restrict analysis to equal
        %amount of trials across predp34conditions
        ERFdata.NumTrialsperpredp34(i_sub, i_Predp34) = length(ERFdata.Index_Predp34{i_sub}{i_Predp34});
    end
    ERFdata.MinNumTrialsperpredp34(i_sub,1) = ...
        min(ERFdata.NumTrialsperpredp34(i_sub, :));
    
    %1.8 Cleanup
    clear temp* DataClean_AllTrials loadfile_ECoGpreprocdata
    
end
disp(['-- Preparing pre-processed data for all subjects finished after: ' ...
    num2str(round(toc/60,2)) 'min --']) %About 20 min for n = 9


%% 2) Load current effect data and aggregate data across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.filt_StimCorr   = [];
Param2plot.all_subs.label_AnatCat   = [];
Param2plot.all_subs.index_AnatCat   = [];
usedElecs_chanposIndex              = [];

for i_sub = 1:length(subs)
    
    tic
    sub = subs{i_sub};
    disp(['Loading data for sub: ' sub])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    %Load prediction data and select current effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
        sub '/Data/Samplewise/'];
    load([path_inputdata ...
        sub '_PredEffectsCluster_' FuncInput_DataType '_' ...
        FuncInput_ToneDur_text 'sTD.mat'], ...
        FuncInput_EffectType , 'labels_loadedData');
    CurrentEffect = eval(FuncInput_EffectType);
    clear PredEffect SimplePredErrEffect ComplexPredErrEffect
    
    %Also load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
    load(loadfile_ECoGpreprocdata);
    
    %Determine basic parameters
    SampleFreq  = DataClean_AllTrials.fsample;
    nSensors    = size(CurrentEffect.stats,2);
    nSamples    = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
    
    %Load stimulus correlation data and select current effect
    if strcmp(param.ElecSelect, 'StimCorr')
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD.mat'], ...
            'corr_ttest', 'SensorLabels');
        filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elec
    else
        filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
    end
    ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
    nSensors_sel            = sum(filt_signelecs_StimCorr);
    
    %aggregate electrode selection across subjects
    Param2plot.all_subs.filt_StimCorr = ...
        [Param2plot.all_subs.filt_StimCorr; filt_signelecs_StimCorr];
    
    %Read electrode labels, coordinates, and anatomical labels for all
    %analyzed electrodes and aggregate them across subjects
    for i_elec = 1:length(ind_signelecs_StimCorr)
        usedElecs_chanposIndex(i_elec,1) = ...
            find(strcmp(...
            labels_loadedData{ind_signelecs_StimCorr(i_elec)}, ...
            DataClean_AllTrials.elec.label));
    end
    coords_sub{i_sub}   = ...
        DataClean_AllTrials.elec.chanpos...
        (usedElecs_chanposIndex,:);
    coords_allsub       = [coords_allsub; coords_sub{i_sub}];
    anatlabels_allsub   = ...
        [anatlabels_allsub; ...
        DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
    
    for i_elec = 1:length(ind_signelecs_StimCorr)
        labels_allsub{end+1,1} = ...
            [labels_loadedData{ind_signelecs_StimCorr(i_elec)} ' ' sub];
    end
    
    %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
    maxnum_cluster = 10; %Across subject max cluster number estimate
    %     maxnum_cluster = 0; %Individual determination doesn't work since different across subjects
    %     for i_elec = 1:nSensors_sel
    %         temp_signcluster = [];
    %         temp_signcluster = ...
    %             sum(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval ...
    %             < param.pval_plotting);
    %         if temp_signcluster > maxnum_cluster
    %            maxnum_cluster = temp_signcluster;
    %         end
    %     end
    
    %Store cluster information from all
    Param2plot.per_sub.clusterstat{i_sub}           = nan(nSensors_sel, maxnum_cluster);
    Param2plot.per_sub.pval{i_sub}                  = nan(nSensors_sel, maxnum_cluster);
    Param2plot.per_sub.pval_derivative{i_sub}       = nan(nSensors_sel, maxnum_cluster);
    Param2plot.per_sub.effect_onset{i_sub}          = nan(nSensors_sel, maxnum_cluster);
    Param2plot.per_sub.effect_offset{i_sub}         = nan(nSensors_sel, maxnum_cluster);
    Param2plot.per_sub.effect_duration{i_sub}       = nan(nSensors_sel, maxnum_cluster);
    
    for i_elec = 1:nSensors_sel
        %Determine clusterorder based on minimal p-value
        clusterorder_minp = [];
        [~, clusterorder_minp] = ...
            sort(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval);
        %For each cluster, read out corresponding information and order it according to minpval
        index_cluster_placement = 0;
        for i_cluster = clusterorder_minp
            index_cluster_placement = index_cluster_placement + 1;
            Param2plot.per_sub.clusterstat{i_sub}(i_elec, index_cluster_placement) = ... %Clusterstat
                CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
            Param2plot.per_sub.pval{i_sub}(i_elec, index_cluster_placement) = ... %pval
                CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
            %If there is a valid cluster,
            if ~isnan(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster))
                %read out cluster timing info
                Param2plot.per_sub.effect_onset{i_sub}(i_elec,index_cluster_placement ) = ... %First sample of cluster to compute relative onset
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1);
                Param2plot.per_sub.effect_offset{i_sub}(i_elec, index_cluster_placement) = ...%Last sample of cluster to compute relative offset
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end);
                Param2plot.per_sub.effect_duration{i_sub}(i_elec, index_cluster_placement) = ...%Relative cluster duration in samples
                    length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster});
            end
        end
    end
    
    %Restrict cluster data to current electrode selection (StimCorr or All)
    for i_elec = 1:nSensors_sel
        if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
            Param2plot.per_sub.pval{i_sub}(i_elec, :) = NaN;
        end
    end
    
    %Compute derivative of p-value to estimate strength of effect
    Param2plot.per_sub.pval_derivative{i_sub} = ...
        -(log10(Param2plot.per_sub.pval{i_sub}));
    
    %Perform FDR_correction on cluster-p-values
    %NOTE: FDR-correction is performed across electrodes for each cluster
    %seperately, since it is much too strict if we treat different clusters
    %as different electrodes/measurements
    if param.FDRcorrect == 1
        FDR_label = 'FDRcorr';
        cluster_pvalFDR_perelec = [];
        for i_cluster = 1:size(Param2plot.per_sub.pval{i_sub},2)
            [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                fdr_bh(Param2plot.per_sub.pval{i_sub}(:,i_cluster), param.pval_FDR, 'pdep','no');
        end
        %Replace uncorrected p-value matrix with FDR-corrected p-values
        Param2plot.per_sub.pval{i_sub} = cluster_pvalFDR_perelec;
    else
        FDR_label = 'uncorr';
    end
    
    %Read out cluster timecourse for sign. clusters only (to later easily
    %use it for plotting sign. clusters)
    for i_elec = 1:length(Param2plot.per_sub.pval{i_sub})
        %Set empty time course as proxy
        Param2plot.per_sub.cluster_timecourse{i_sub}(i_elec, :) = ...
            zeros(1,length(CurrentEffect.clusterstat...
            {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse));
        
        for i_cluster = 1:length(Param2plot.per_sub.pval{i_sub}(i_elec,:))
            %If a cluster is sign, enter its timecourse in the blank proxy
            %nd denote it with the cluster index
            if Param2plot.per_sub.pval{i_sub}(i_elec,i_cluster) < param.pval_plotting
                Param2plot.per_sub.cluster_timecourse{i_sub}...
                    (i_elec, find(CurrentEffect.clusterstat...
                    {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse ...
                    == i_cluster)) = i_cluster;
            end
        end
    end
    
    %Data aggregation over subjects
    %Array to differentiate subject entries
    Param2plot.all_subs.sub_index = ...
        [Param2plot.all_subs.sub_index; ...
        ones(length(coords_sub{i_sub}),1)*i_sub];
    
    %Clusterinfo aggregation over subjects (Dim: elec * cluster)
    if i_sub == 1
        Param2plot.all_subs.clusterstat     = ...
            Param2plot.per_sub.clusterstat{i_sub};
        Param2plot.all_subs.pval            = ...
            Param2plot.per_sub.pval{i_sub};
        Param2plot.all_subs.pval_derivative = ...
            Param2plot.per_sub.pval_derivative{i_sub};
        Param2plot.all_subs.effect_onset    = ...
            Param2plot.per_sub.effect_onset{i_sub};
        Param2plot.all_subs.effect_offset   = ...
            Param2plot.per_sub.effect_offset{i_sub};
        Param2plot.all_subs.effect_duration = ...
            Param2plot.per_sub.effect_duration{i_sub};
    else
        Param2plot.all_subs.clusterstat     = ...
            [Param2plot.all_subs.clusterstat; ...
            Param2plot.per_sub.clusterstat{i_sub}];
        Param2plot.all_subs.pval            = ...
            [Param2plot.all_subs.pval; ...
            Param2plot.per_sub.pval{i_sub}];
        Param2plot.all_subs.pval_derivative = ...
            [Param2plot.all_subs.pval_derivative; ...
            Param2plot.per_sub.pval_derivative{i_sub}];
        Param2plot.all_subs.effect_onset    = ...
            [Param2plot.all_subs.effect_onset; ...
            Param2plot.per_sub.effect_onset{i_sub}];
        Param2plot.all_subs.effect_offset   = ...
            [Param2plot.all_subs.effect_offset; ...
            Param2plot.per_sub.effect_offset{i_sub}];
        Param2plot.all_subs.effect_duration = ...
            [Param2plot.all_subs.effect_duration; ...
            Param2plot.per_sub.effect_duration{i_sub}];
    end
    
    
    %% 3) Categorize electrodes according to anatomical regions
    AnatReg_allSubs{i_sub} = ...
        NASTD_ECoG_AssignAnatRegions(DataClean_AllTrials, labels_loadedData);
    
    %Restrict anatomical parcellation output to selected electrodes
    %(e.g., sequence tracking selection or all elecs)
    temp_filt_allelecs_currsub = ...
        find(Param2plot.all_subs.sub_index == i_sub);
    temp_filt_selelecs_currsub = ...
        logical(Param2plot.all_subs.filt_StimCorr(temp_filt_allelecs_currsub));
    
    AnatReg_allSubs{i_sub}.CatIndex = ...
        AnatReg_allSubs{i_sub}.CatIndex(temp_filt_selelecs_currsub,:);
    AnatReg_allSubs{i_sub}.Info_perelec = ...
        AnatReg_allSubs{i_sub}.Info_perelec(temp_filt_selelecs_currsub,:);
    
    %Aggregate category labels and indices across subjects
    Param2plot.all_subs.label_AnatCat  = ...
        [Param2plot.all_subs.label_AnatCat; ...
        AnatReg_allSubs{i_sub}.Info_perelec(:,3)];
    Param2plot.all_subs.index_AnatCat  = ...
        [Param2plot.all_subs.index_AnatCat; ...
        AnatReg_allSubs{i_sub}.CatIndex];
    
    %Cleanup
    usedElecs_chanposIndex = [];
    clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
    disp(['done loading in ' num2str(toc) ' sec'])
end


%% 4) Determine sign. electrodes
SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
SignElecs.index         = find(SignElecs.array);%1D index array denoting sign. elecs (independent of number of clusters)
SignElecs.num_elecs     = length(SignElecs.index);
SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

%Determine anatomical parcellations for sign. electrodes
SignElecs.label_AnatCat = Param2plot.all_subs.label_AnatCat(SignElecs.index);
SignElecs.index_AnatCat = max(Param2plot.all_subs.index_AnatCat(SignElecs.index,:),[],2);

%Determine labels of sign. elecs and add number to subtitle
SignElecs.labels        = labels_allsub(SignElecs.index);
SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
SignElecs.fulllabels    = [];
for i_elec = 1:length(SignElecs.labels)
    SignElecs.fulllabels{i_elec,1} = ...
        [SignElecs.labels{i_elec}, ' ', SignElecs.anatlabels{i_elec}];
end

if ~isempty(SignElecs.labels)
    sign_title = ...
        [num2str(SignElecs.num_elecs)...
        ' / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
else
    sign_title = ...
        ['0 / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
end

%Create electrode labels for plotting
for i_elec = 1:length(labels_allsub)%No electrode labels
    labels_plotting_empty{i_elec} = '';
end
counter_elecs           = 0;
labels_plotting_number  = labels_plotting_empty;
for i_elec = SignElecs.index'
    counter_elecs = counter_elecs +1;
    labels_plotting_number{i_elec} = num2str(counter_elecs);
end


if ~isempty(SignElecs.index)
    %% 5) Plot figure showing surf with sign. elecs and p33 ERFs per p*34 condition
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
    end
    
    %Colorlim
    clims               = [0 max(PlotInput)];
    Label_Colorbar      = '- log10 (cluster p-value) ';
    
    dat.dimord          = 'chan_time';
    dat.time            = 0;
    dat.label           = labels_allsub;
    dat.avg             = PlotInput;
    dat.sign_elecs      = SignElecs.array;
    dat.textcolor_rgb   = [1 0 0];
    SizeFactor          = 4;
    
    chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    cmap        = 'parula';
    
    SubplotPosition = [0 0 0 0];
    
    %Project all electrodes on one hemisphere,
    coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
    
    %Plot surface
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
        chanSize, clims, cmap, dat.textcolor_rgb, ...
        DimSubplot, [1 2 DimSubplot(1)+1 DimSubplot(1)+2], SubplotPosition, [], []);
    
    sp_handle_surf = sp_handle_surf_temp.L;
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder'
        subplot_counter = subplot_counter + 1;
        
        %Determine information for current electrode
        ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
        ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
        ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
            find(strcmp(ERFdata.elec_labels{ind_curr_sub}, ...
            strtok(SignElecs.labels{i_signelec},' ')));
        label_curr_elec = SignElecs.labels{i_signelec};
        
        %For current electrode, Compute AVG + STD across trials per Predp34 condition
        for i_Predp34 = 1:length(label_Predp34)
            temp_ERFdata_perp34 = [];
            for i_trial = 1:length(ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34})
                temp_ERFdata_perp34(i_trial,:) = ...
                    ERFdata.p33{ind_curr_sub}...
                    (ind_curr_elec_currsub,:,ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34}(i_trial));
            end
            Avgp33{i_Predp34} = nanmean(temp_ERFdata_perp34,1);
            STDp33{i_Predp34} = nanstd(temp_ERFdata_perp34,1,1);
        end
        
        %Determine trace colors
        plot_param.color = ...
            {[0, 0.4470, 0.7410, 0.7], ...
            [0, 0.75, 0.75, 0.7],[0.8500, ...
            0.3250, 0.0980, 0.7]};
        
        %Plot ERF for sign  elecs
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
        
        plot(1:length(Avgp33{1}(:)),Avgp33{1}(:),...
            'color',plot_param.color{1},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{2}(:)),Avgp33{2}(:),...
            'color',plot_param.color{2},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{3}(:)),Avgp33{3}(:),...
            'color',plot_param.color{3},'LineWidth', 2)
        axis tight
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('p33 [s]')
        
        %Highlight sign. samples/clusters
        grey  = [100 100 100]./255;
        maxY = max([max(Avgp33{1}(:)), max(Avgp33{2}(:)), max(Avgp33{3}(:))]);
        minY = min([min(Avgp33{1}(:)), min(Avgp33{2}(:)), min(Avgp33{3}(:))]);
        
        sign_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub}(ind_curr_elec_currsub,:) > 0;
        area(1:length(sign_samples), ...
            sign_samples * ...
            (max([minY maxY])+abs(max([minY maxY])*0.05)),...
            'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        if min([minY maxY]) < 0
            area(1:length(sign_samples), ...
                sign_samples * ...
                (min([minY maxY])-abs(min([minY maxY])*0.05)),...
                'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        end
        
        %Add cluster pval as text to shading
        for i_cluster = 1:sum(~isnan(...
                Param2plot.per_sub.pval{ind_curr_sub}(ind_curr_elec_currsub,:)))
            if Param2plot.per_sub.pval{ind_curr_sub}(ind_curr_elec_currsub,i_cluster) ...
                    < param.pval_plotting
                shading_startsample = ...
                    find(Param2plot.per_sub.cluster_timecourse{ind_curr_sub}...
                    (ind_curr_elec_currsub,:) == i_cluster);
                text(shading_startsample(1),max([minY maxY] * 0.8), ...
                    [' p = ' num2str(round(...
                    Param2plot.per_sub.pval{ind_curr_sub}...
                    (ind_curr_elec_currsub, i_cluster),3))], ...
                    'FontSize',8);
            end
        end
        
        %Scale y-axis to show sign. text if present
        ylim([min([minY maxY])-+abs(max([minY maxY])*0.05) ...
            max([minY maxY])+abs(max([minY maxY])*0.05)]); %With space for sign. info
        
        %subplot title
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec],'FontSize',8,'Interpreter','none')
        
    end
    
    %Adjust figureheader/title
    switch FuncInput_EffectType
        case 'PredEffect'
            effect_text = ['Effect: Prediction (explain p33 activity by p*34);' ...
                ' p < ' num2str(param.pval_plotting) ' ' FDR_label];
        case 'SimplePredErrEffect'
            effect_text = ['Effect: Simple Prediction error' ...
                ' (explain p34 activity by absolute p33-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
        case 'ComplexPredErrEffect'
            effect_text = ['Effect: Complex Prediction error' ...
                ' (explain p34 activity by absolute p*34-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
    end
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - ' effect_text] ...
        ['Input data: ' FuncInput_DataType ', ' FuncInput_ToneDur_text 'ms TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp33_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' ...
            FuncInput_EffectType '_' FuncInput_DataType ...
            '_' FuncInput_ToneDur_text 'msTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
    
    %% 6) Plot figure showing surf with sign. elecs and p1-p33 ERFs per p*34 condition
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','off'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
    end
    
    %Colorlim
    clims               = [0 max(PlotInput)];
    Label_Colorbar      = '- log10 (cluster p-value) ';
    
    dat.dimord          = 'chan_time';
    dat.time            = 0;
    dat.label           = labels_allsub;
    dat.avg             = PlotInput;
    dat.sign_elecs      = SignElecs.array;
    dat.textcolor_rgb   = [1 0 0];
    SizeFactor          = 4;
    
    chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    cmap        = 'parula';
    
    SubplotPosition = [0 0 0 0];
    
    %Project all electrodes on one hemisphere,
    coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
    
    %Plot surface
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
        chanSize, clims, cmap, dat.textcolor_rgb, ...
        DimSubplot, [1 2 DimSubplot(1)+1 DimSubplot(1)+2], SubplotPosition, [], []);
    
    sp_handle_surf = sp_handle_surf_temp.L;
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder'
        subplot_counter = subplot_counter + 1;
        
        %Determine information for current electrode
        ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
        ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
        ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
            find(strcmp(ERFdata.elec_labels{ind_curr_sub}, ...
            strtok(SignElecs.labels{i_signelec},' ')));
        label_curr_elec = SignElecs.labels{i_signelec};
        
        %For current electrode, Compute AVG + STD across trials per Predp34 condition
        for i_Predp34 = 1:length(label_Predp34)
            temp_ERFdata_perp34 = [];
            for i_trial = 1:length(ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34})
                temp_ERFdata_perp34(i_trial,:) = ...
                    ERFdata.p1p33{ind_curr_sub}...
                    (ind_curr_elec_currsub,:,ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34}(i_trial));
            end
            Avgp1p33{i_Predp34} = nanmean(temp_ERFdata_perp34,1);
            STDp1p33{i_Predp34} = nanstd(temp_ERFdata_perp34,1,1);
        end
        
        %Determine trace colors
        plot_param.color = ...
            {[0, 0.4470, 0.7410, 0.7], ...
            [0, 0.75, 0.75, 0.7],[0.8500, ...
            0.3250, 0.0980, 0.7]};
        
        %Plot ERF for sign  elecs
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
        
        plot(1:length(Avgp1p33{1}(:)),Avgp1p33{1}(:),...
            'color',plot_param.color{1},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp1p33{2}(:)),Avgp1p33{2}(:),...
            'color',plot_param.color{2},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp1p33{3}(:)),Avgp1p33{3}(:),...
            'color',plot_param.color{3},'LineWidth', 2)
        axis tight
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('p1-p33 [s]')
        
        %     %Highlight sign. samples/clusters
        %     grey  = [100 100 100]./255;
        %     maxY = max([max(Avgp1p33{1}(:)), max(Avgp1p33{2}(:)), max(Avgp1p33{3}(:))]);
        %     minY = min([min(Avgp1p33{1}(:)), min(Avgp1p33{2}(:)), min(Avgp1p33{3}(:))]);
        %
        %     sign_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub}(ind_curr_elec_currsub,:) > 0;
        %     %Adjust to entire seq duration (not only p33)
        %     sign_samples = ...
        %         [zeros(1,StimTiming.Sample_Tone_StartStop{ind_curr_sub}(param.ToneIndex,1) - ...
        %         StimTiming.Sample_Tone_StartStop{ind_curr_sub}(1,1)) sign_samples];
        %     area(1:length(sign_samples), ...
        %         sign_samples * ...
        %         (max([minY maxY])+abs(max([minY maxY])*0.05)),...
        %         'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        %     if min([minY maxY]) < 0
        %         area(1:length(sign_samples), ...
        %             sign_samples * ...
        %             (min([minY maxY])-abs(min([minY maxY])*0.05)),...
        %             'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        %     end
        
        %     %Add cluster pval as text to shading
        %     for i_cluster = 1:sum(~isnan(...
        %             Param2plot.per_sub.pval{ind_curr_sub}(ind_curr_elec_currsub,:)))
        %         if Param2plot.per_sub.pval{ind_curr_sub}(ind_curr_elec_currsub,i_cluster) ...
        %                 < param.pval_plotting
        %             shading_startsample = ...
        %                 find(Param2plot.per_sub.cluster_timecourse{ind_curr_sub}...
        %                 (ind_curr_elec_currsub,:) == i_cluster);
        %             %Adjust to entire seq duration (not only p33)
        %             shading_startsample = ...
        %                 shading_startsample + ...
        %                 (StimTiming.Sample_Tone_StartStop{ind_curr_sub}(param.ToneIndex,1) - ...
        %                 StimTiming.Sample_Tone_StartStop{ind_curr_sub}(1,1));
        %             text(shading_startsample(1),max([minY maxY] * 0.8), ...
        %                 [' p = ' num2str(round(...
        %                 Param2plot.per_sub.pval{ind_curr_sub}...
        %                 (ind_curr_elec_currsub, i_cluster),3))], ...
        %                 'FontSize',8);
        %         end
        %     end
        
        %Add demarcation lines for tones
        a = axis;
        %     for i_tone = 1:size(StimTiming.Sample_Tone_StartStop{ind_curr_sub},1)
        %         ind_start = StimTiming.Sample_Tone_StartStop{ind_curr_sub}(i_tone, 1) - StimTiming.Sample_Tone_StartStop{ind_curr_sub}(1, 1);
        %         hold on;
        %         plot([ind_start ind_start],[a(3) a(4)], 'Color',[0.5 0.5 0.5])
        %     end
        for i_tone = [10 20 30 33]
            ind_start = StimTiming.Sample_Tone_StartStop{ind_curr_sub}(i_tone, 1) - StimTiming.Sample_Tone_StartStop{ind_curr_sub}(1, 1);
            hold on;
            plot([ind_start ind_start],[a(3) a(4)], 'Color',[0 0 0], 'LineWidth', 1)
        end
        
        %subplot title
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec],'FontSize',8,'Interpreter','none')
        
    end
    
    %Adjust figureheader/title
    switch FuncInput_EffectType
        case 'PredEffect'
            effect_text = ['Effect: Prediction (explain p33 activity by p*34);' ...
                ' p < ' num2str(param.pval_plotting) ' ' FDR_label];
        case 'SimplePredErrEffect'
            effect_text = ['Effect: Simple Prediction error' ...
                ' (explain p34 activity by absolute p33-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
        case 'ComplexPredErrEffect'
            effect_text = ['Effect: Complex Prediction error' ...
                ' (explain p34 activity by absolute p*34-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
    end
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - ' effect_text] ...
        ['Input data: ' FuncInput_DataType ', ' FuncInput_ToneDur_text 'ms TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp1p33_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' ...
            FuncInput_EffectType '_' FuncInput_DataType ...
            '_' FuncInput_ToneDur_text 'msTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
%     %% 7) Plot figure showing surf with sign. elecs and p34 ERFs per p*34 condition
%     %Determine number of subplots (surf + 1 subplot per sign. elec)
%     num_subplot = SignElecs.num_elecs + 4;
%     DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
%     SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
%     SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
%     
%     %Prepare figure
%     h = figure('visible','off'); %ensures figure doesn't pop up during plotting
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     set(gcf,'renderer','opengl');
%     
%     %Determine parameter to be plotted and set up plotting struct
%     for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
%         PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
%     end
%     
%     %Colorlim
%     clims               = [0 max(PlotInput)];
%     Label_Colorbar      = '- log10 (cluster p-value) ';
%     
%     dat.dimord          = 'chan_time';
%     dat.time            = 0;
%     dat.label           = labels_allsub;
%     dat.avg             = PlotInput;
%     dat.sign_elecs      = SignElecs.array;
%     dat.textcolor_rgb   = [1 0 0];
%     SizeFactor          = 4;
%     
%     chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
%     cmap        = 'parula';
%     
%     SubplotPosition = [0 0 0 0];
%     
%     %Project all electrodes on one hemisphere,
%     coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
%     
%     %Plot surface
%     sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
%         (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
%         chanSize, clims, cmap, dat.textcolor_rgb, ...
%         DimSubplot, [1 2 DimSubplot(1)+1 DimSubplot(1)+2], SubplotPosition, [], []);
%     
%     sp_handle_surf = sp_handle_surf_temp.L;
%     
%     %Colorbar
%     ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
%     h = colorbar;
%     h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
%     h.Label.String      = Label_Colorbar;
%     h.FontSize          = 12;
%     caxis(clims)
%     
%     %Plot ERFs for each sign. electrode
%     %Order sign. elecs based on anatomical parcellation
%     [~, index_anatorder] = sort(SignElecs.index_AnatCat);
%     subplot_counter = 0;
%     for i_signelec = index_anatorder'
%         subplot_counter = subplot_counter + 1;
%         
%         %Determine information for current electrode
%         ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
%         ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
%         ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
%             find(strcmp(ERFdata.elec_labels{ind_curr_sub}, ...
%             strtok(SignElecs.labels{i_signelec},' ')));
%         label_curr_elec = SignElecs.labels{i_signelec};
%         
%         %For current electrode, Compute AVG + STD across trials per Predp34 condition
%         for i_Predp34 = 1:length(label_Predp34)
%             temp_ERFdata_perp34 = [];
%             for i_trial = 1:length(ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34})
%                 temp_ERFdata_perp34(i_trial,:) = ...
%                     ERFdata.p34{ind_curr_sub}...
%                     (ind_curr_elec_currsub,:,ERFdata.Index_Predp34{ind_curr_sub}{i_Predp34}(i_trial));
%             end
%             Avgp34{i_Predp34} = nanmean(temp_ERFdata_perp34,1);
%             STDp34{i_Predp34} = nanstd(temp_ERFdata_perp34,1,1);
%         end
%         
%         %Determine trace colors
%         plot_param.color = ...
%             {[0, 0.4470, 0.7410, 0.7], ...
%             [0, 0.75, 0.75, 0.7],[0.8500, ...
%             0.3250, 0.0980, 0.7]};
%         
%         %Plot ERF for sign  elecs
%         subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
%         
%         plot(1:length(Avgp34{1}(:)),Avgp34{1}(:),...
%             'color',plot_param.color{1},'LineWidth', 2)
%         hold on;
%         plot(1:length(Avgp34{2}(:)),Avgp34{2}(:),...
%             'color',plot_param.color{2},'LineWidth', 2)
%         hold on;
%         plot(1:length(Avgp34{3}(:)),Avgp34{3}(:),...
%             'color',plot_param.color{3},'LineWidth', 2)
%         axis tight
%         
%         %Change x-axis labeling to time [s]
%         temp_xticks = xticks;
%         temp_text = {};
%         for i_xtick = 1:length(temp_xticks)
%             temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
%         end
%         xticklabels(temp_text)
%         xlabel('p34 [s]')
%         
%         %subplot title
%         title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
%             AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
%             label_curr_elec],'FontSize',8,'Interpreter','none')
%     end
%     
%     %Adjust figureheader/title
%     switch FuncInput_EffectType
%         case 'PredEffect'
%             effect_text = ['Effect: Prediction (explain p33 activity by p*34);' ...
%                 ' p < ' num2str(param.pval_plotting) ' ' FDR_label];
%         case 'SimplePredErrEffect'
%             effect_text = ['Effect: Simple Prediction error' ...
%                 ' (explain p34 activity by absolute p33-p34 difference); p < ' ...
%                 num2str(param.pval_plotting) ' ' FDR_label];
%         case 'ComplexPredErrEffect'
%             effect_text = ['Effect: Complex Prediction error' ...
%                 ' (explain p34 activity by absolute p*34-p34 difference); p < ' ...
%                 num2str(param.pval_plotting) ' ' FDR_label];
%     end
%     
%     Fig_title = {['Group level (n = ' num2str(length(subs)) ') - ' effect_text] ...
%         ['Input data: ' FuncInput_DataType ', ' FuncInput_ToneDur_text 'ms TD - Output: ' ...
%         num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
%     sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
%     
%     if save_poststepFigs == 1
%         filename     = ['Surf1HERFp34_SignElec' param.ElecSelect ...
%             '_Allsubn' num2str(length(subs)) '_' ...
%             FuncInput_EffectType '_' FuncInput_DataType ...
%             '_' FuncInput_ToneDur_text 'msTD_Stat' FDR_label '.png'];
%         figfile      = [path_fig filename];
%         saveas(gcf, figfile, 'png'); %save png version
%         close all;
%     end
end