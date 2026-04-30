function NASTD_ECoG_Connectivity_PlotSignElecConnectionsforGC_AllSubTD...
    (subs, ...
    FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, plot_poststepFigs, ...
    paths_NASTD_ECoG)

%Aim: Read out electrodes showing a sign. Prediction (P) or Complex
%Prediciton Error (PE) for each subject, TD, and Input Data Type.
%These electrodes will form the basis for GC analysis.

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

folderlabel_inputsignal = '';
for i = 1:length(FuncInput_DataType)
    folderlabel_inputsignal = strcat(folderlabel_inputsignal, FuncInput_DataType{i});
end
folderlabel_inputsignal = strcat(folderlabel_inputsignal, '-p', num2str(param.pval_plotting), FDR_label);
folderlabel_inputsignal = erase(folderlabel_inputsignal, '_');

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity ...
    'ElecSelect/Allsub_n' num2str(length(subs)) ...
    '/Figs/ElecPairs/' folderlabel_inputsignal '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load prediction effect data and aggregate relevant info across subjects
clear PredictionEffect
usedElecs_chanposIndex                      = [];

tic
for i_inputdatatyp = 1:length(FuncInput_DataType)
    
    PredictionEffect{i_inputdatatyp}.all_subs.sub_index         = [];
    PredictionEffect{i_inputdatatyp}.all_subs.TD_index          = [];
    PredictionEffect{i_inputdatatyp}.all_subs.label_AnatCat     = [];
    PredictionEffect{i_inputdatatyp}.all_subs.index_AnatCat     = [];
    
    PredictionEffect{i_inputdatatyp}.all_subs.label_elec        = [];
    PredictionEffect{i_inputdatatyp}.all_subs.label_anat        = [];
    PredictionEffect{i_inputdatatyp}.all_subs.coords_elec       = [];
    
    for i_sub = 1:length(subs)
        
        sub = subs{i_sub};
        disp(['-- Loading data for sub: ' sub ' --'])
        NASTD_ECoG_subjectinfo %load subject info file (var: si)
        subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
        
        for i_TD = 1:length(FuncInput_ToneDur_text)
            %Load prediction data and select prediction and PE effect
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
                sub '/Data/Samplewise/'];
            
            curr_inputdatatype = FuncInput_DataType{i_inputdatatyp};
            
            load([path_inputdata sub '_PredEffectsCluster_' ...
                curr_inputdatatype '_' ...
                FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'PredEffect' , 'labels_loadedData');
            CurrentEffect = PredEffect;
            clear PredEffect
            
            %Also load ECoG preproc data for channel labels and position
            loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
            load(loadfile_ECoGpreprocdata);
            
            SampleFreq              = DataClean_AllTrials.fsample;
            nSensors_all            = size(CurrentEffect.stats,2);
            nSamples(i_TD)          = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
            
            %Load stimulus correlation data and select current effect
            if strcmp(param.ElecSelect, 'StimCorr')
                path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
                load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
                    FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                    'corr_ttest', 'SensorLabels');
                filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
            else
                filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
            end
            ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
            nSensors_sel            = sum(filt_signelecs_StimCorr);
            
            PredictionEffect{i_inputdatatyp}.per_sub.elec_labels{i_sub} = ...
                labels_loadedData;
            
            %Read electrode labels, coordinates, and anatomical labels for all
            %sign. StimCorr electrodes and aggregate them across subjects
            for i_elec = 1:length(ind_signelecs_StimCorr)
                usedElecs_chanposIndex(i_elec,1) = ...
                    find(strcmp(labels_loadedData{ind_signelecs_StimCorr(i_elec)}, ...
                    DataClean_AllTrials.elec.label));
            end
            coords_sub{i_sub, i_TD}   = ...
                DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);
            PredictionEffect{i_inputdatatyp}.all_subs.coords_elec = ...
                [PredictionEffect{i_inputdatatyp}.all_subs.coords_elec; coords_sub{i_sub, i_TD}];
            PredictionEffect{i_inputdatatyp}.all_subs.label_anat = ...
                [PredictionEffect{i_inputdatatyp}.all_subs.label_anat; DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
            for i_elec = 1:length(ind_signelecs_StimCorr)
                PredictionEffect{i_inputdatatyp}.all_subs.label_elec{end+1,1} = ...
                    [labels_loadedData{ind_signelecs_StimCorr(i_elec)} ' ' sub];
            end
            
            %Categorize electrodes according to anatomical regions
            if i_inputdatatyp == 1
                AnatReg_allSubs{i_sub, i_TD} = ...
                    NASTD_ECoG_AssignAnatRegions...
                    (DataClean_AllTrials, labels_loadedData);
            end
            
            %Aggregate category labels and indices across subjects
            PredictionEffect{i_inputdatatyp}.all_subs.label_AnatCat  = ...
                [PredictionEffect{i_inputdatatyp}.all_subs.label_AnatCat; ...
                AnatReg_allSubs{i_sub, i_TD}.Info_perelec(:,3)];
            PredictionEffect{i_inputdatatyp}.all_subs.index_AnatCat  = ...
                [PredictionEffect{i_inputdatatyp}.all_subs.index_AnatCat; ...
                AnatReg_allSubs{i_sub, i_TD}.CatIndex];
            
            %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
            maxnum_cluster = 10; %Across subject max cluster number estimate
            %Store cluster information from all
            PredictionEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}           = nan(nSensors_sel, maxnum_cluster);
            PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}                  = nan(nSensors_sel, maxnum_cluster);
            PredictionEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
            PredictionEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}          = nan(nSensors_sel, maxnum_cluster);
            PredictionEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}         = nan(nSensors_sel, maxnum_cluster);
            PredictionEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
            
            for i_elec = 1:nSensors_sel
                %Determine clusterorder based on minimal p-value
                clusterorder_minp = [];
                [~, clusterorder_minp] = ...
                    sort(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval);
                %For each cluster, read out corresponding information and order it according to minpval
                index_cluster_placement = 0;
                for i_cluster = clusterorder_minp
                    index_cluster_placement = index_cluster_placement + 1;
                    PredictionEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %Clusterstat
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
                    PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %pval
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
                    %If there is a valid cluster,
                    if ~isnan(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster))
                        %read out cluster timing info
                        PredictionEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}(i_elec,index_cluster_placement ) = ... %First sample of cluster to compute relative onset
                            CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1) / nSamples(i_TD);
                        PredictionEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Last sample of cluster to compute relative offset
                            CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end) / nSamples(i_TD);
                        PredictionEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Relative cluster duration in samples
                            length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}) / nSamples(i_TD);
                    end
                end
            end
            
            %Restrict cluster data to current electrode selection (StimCorr or All)
            for i_elec = 1:nSensors_sel
                if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
                    PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec, :) = NaN;
                end
            end
            
            %Compute derivative of p-value to estimate strength of effect
            PredictionEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD} = ...
                -(log10(PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}));
            
            %Perform FDR_correction on cluster-p-values
            %NOTE: FDR-correction is performed across electrodes for each cluster
            %seperately, since it is much too strict if we treat different clusters
            %as different electrodes/measurements
            if param.FDRcorrect == 1
                FDR_label = 'FDRcorr';
                cluster_pvalFDR_perelec = [];
                for i_cluster = 1:size(PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD},2)
                    [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                        fdr_bh(PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(:,i_cluster), param.pval_FDR, 'pdep','no');
                end
                
                disp(['p-val FDR: ' num2str(cluster_critpFDR)])
                disp(['p-val Bonferroni: ' num2str(param.pval_plotting/nSensors_all)])
                
                %Replace uncorrected p-value matrix with FDR-corrected p-values
                PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD} = cluster_pvalFDR_perelec;
            else
                FDR_label = 'uncorr';
            end
            
            %Read out cluster timecourse for sign. clusters only
            %(to later easily use it for plotting sign. clusters)
            for i_elec = 1:length(PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD})
                %Set empty time course as proxy
                PredictionEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec, :) = ...
                    zeros(1,length(CurrentEffect.clusterstat...
                    {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse));
                
                for i_cluster = 1:length(PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec,:))
                    %If a cluster is sign, enter its timecourse in the blank proxy
                    %nd denote it with the cluster index
                    if PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec,i_cluster) < param.pval_plotting
                        PredictionEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}...
                            (i_elec, find(CurrentEffect.clusterstat...
                            {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse ...
                            == i_cluster)) = i_cluster;
                    end
                end
            end
            
            %Data aggregation over subjects
            %Array to differentiate subject entries
            PredictionEffect{i_inputdatatyp}.all_subs.sub_index = ...
                [PredictionEffect{i_inputdatatyp}.all_subs.sub_index; ...
                ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];
            %Array to differentiate TD entries for each electrode
            for i_elec = 1:nSensors_all
                PredictionEffect{i_inputdatatyp}.all_subs.TD_index    = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.TD_index; ...
                    FuncInput_ToneDur_text{i_TD} 's'];
            end
            
            %Clusterinfo aggregation over subjects (Dim: elec * cluster)
            if i_sub == 1 && i_TD == 1
                PredictionEffect{i_inputdatatyp}.all_subs.clusterstat     = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD};
                PredictionEffect{i_inputdatatyp}.all_subs.pval            = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD};
                PredictionEffect{i_inputdatatyp}.all_subs.pval_derivative = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD};
                PredictionEffect{i_inputdatatyp}.all_subs.effect_onset    = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD};
                PredictionEffect{i_inputdatatyp}.all_subs.effect_offset   = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD};
                PredictionEffect{i_inputdatatyp}.all_subs.effect_duration = ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD};
            else
                PredictionEffect{i_inputdatatyp}.all_subs.clusterstat     = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.clusterstat; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}];
                PredictionEffect{i_inputdatatyp}.all_subs.pval            = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.pval; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}];
                PredictionEffect{i_inputdatatyp}.all_subs.pval_derivative = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.pval_derivative; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD}];
                PredictionEffect{i_inputdatatyp}.all_subs.effect_onset    = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.effect_onset; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}];
                PredictionEffect{i_inputdatatyp}.all_subs.effect_offset   = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.effect_offset; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}];
                PredictionEffect{i_inputdatatyp}.all_subs.effect_duration = ...
                    [PredictionEffect{i_inputdatatyp}.all_subs.effect_duration; ...
                    PredictionEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}];
            end
            
            %Cleanup
            usedElecs_chanposIndex = [];
            clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
        end
        disp([' -- Prediction electrodes processed for ' sub ' (both TW, ' FuncInput_DataType{i_inputdatatyp} ') --'])
    end
end
disp(['-- Prediction electrodes read out for all input data types, subs, and TW after ' num2str(round(toc/60),2) ' minutes --'])


%% 2) Load cmplex PE effect data and aggregate relevant info across subjects
clear PEEffect
usedElecs_chanposIndex              = [];
tic
for i_inputdatatyp = 1:length(FuncInput_DataType)
    
    PEEffect{i_inputdatatyp}.all_subs.sub_index       = [];
    PEEffect{i_inputdatatyp}.all_subs.TD_index        = [];
    PEEffect{i_inputdatatyp}.all_subs.label_AnatCat   = [];
    PEEffect{i_inputdatatyp}.all_subs.index_AnatCat   = [];
    
    for i_sub = 1:length(subs)
        
        sub = subs{i_sub};
        disp(['-- Loading data for sub: ' sub ' --'])
        NASTD_ECoG_subjectinfo %load subject info file (var: si)
        subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
        
        for i_TD = 1:length(FuncInput_ToneDur_text)
            %Load prediction data and select prediction and PE effect
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
                sub '/Data/Samplewise/'];
            
            curr_inputdatatype = FuncInput_DataType{i_inputdatatyp};
            
            load([path_inputdata sub '_PredEffectsCluster_' ...
                curr_inputdatatype '_' ...
                FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'ComplexPredErrEffect' , 'labels_loadedData');
            CurrentEffect = ComplexPredErrEffect;
            clear ComplexPredErrEffect
            
            %Also load ECoG preproc data for channel labels and position
            loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
            load(loadfile_ECoGpreprocdata);
            
            SampleFreq              = DataClean_AllTrials.fsample;
            nSensors_all            = size(CurrentEffect.stats,2);
            nSamples(i_TD)          = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
            
            %Load stimulus correlation data and select current effect
            if strcmp(param.ElecSelect, 'StimCorr')
                path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
                load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
                    FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                    'corr_ttest', 'SensorLabels');
                filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elecs
            else
                filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
            end
            ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
            nSensors_sel            = sum(filt_signelecs_StimCorr);
            
            PEEffect{i_inputdatatyp}.per_sub.elec_labels{i_sub} = ...
                labels_loadedData;
            
            %Aggregate category labels and indices across subjects
            PEEffect{i_inputdatatyp}.all_subs.label_AnatCat  = ...
                [PEEffect{i_inputdatatyp}.all_subs.label_AnatCat; ...
                AnatReg_allSubs{i_sub, i_TD}.Info_perelec(:,3)];
            PEEffect{i_inputdatatyp}.all_subs.index_AnatCat  = ...
                [PEEffect{i_inputdatatyp}.all_subs.index_AnatCat; ...
                AnatReg_allSubs{i_sub, i_TD}.CatIndex];
            
            %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
            maxnum_cluster = 10; %Across subject max cluster number estimate
            
            %Store cluster information from all
            PEEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}           = nan(nSensors_sel, maxnum_cluster);
            PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}                  = nan(nSensors_sel, maxnum_cluster);
            PEEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
            PEEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}          = nan(nSensors_sel, maxnum_cluster);
            PEEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}         = nan(nSensors_sel, maxnum_cluster);
            PEEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
            
            for i_elec = 1:nSensors_sel
                %Determine clusterorder based on minimal p-value
                clusterorder_minp = [];
                [~, clusterorder_minp] = ...
                    sort(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval);
                %For each cluster, read out corresponding information and order it according to minpval
                index_cluster_placement = 0;
                for i_cluster = clusterorder_minp
                    index_cluster_placement = index_cluster_placement + 1;
                    PEEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %Clusterstat
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
                    PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %pval
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
                    %If there is a valid cluster,
                    if ~isnan(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster))
                        %read out cluster timing info
                        PEEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}(i_elec,index_cluster_placement ) = ... %First sample of cluster to compute relative onset
                            CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1) / nSamples(i_TD);
                        PEEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Last sample of cluster to compute relative offset
                            CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end) / nSamples(i_TD);
                        PEEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Relative cluster duration in samples
                            length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}) / nSamples(i_TD);
                    end
                end
            end
            
            %Restrict cluster data to current electrode selection (StimCorr or All)
            for i_elec = 1:nSensors_sel
                if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
                    PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec, :) = NaN;
                end
            end
            
            %Compute derivative of p-value to estimate strength of effect
            PEEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD} = ...
                -(log10(PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}));
            
            %Perform FDR_correction on cluster-p-values
            %NOTE: FDR-correction is performed across electrodes for each cluster
            %seperately, since it is much too strict if we treat different clusters
            %as different electrodes/measurements
            if param.FDRcorrect == 1
                FDR_label = 'FDRcorr';
                cluster_pvalFDR_perelec = [];
                for i_cluster = 1:size(PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD},2)
                    [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                        fdr_bh(PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(:,i_cluster), param.pval_FDR, 'pdep','no');
                end
                %Replace uncorrected p-value matrix with FDR-corrected p-values
                PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD} = cluster_pvalFDR_perelec;
            else
                FDR_label = 'uncorr';
            end
            
            %Read out cluster timecourse for sign. clusters only (to later easily
            %use it for plotting sign. clusters)
            for i_elec = 1:length(PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD})
                %Set empty time course as proxy
                PEEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec, :) = ...
                    zeros(1,length(CurrentEffect.clusterstat...
                    {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse));
                
                for i_cluster = 1:length(PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec,:))
                    %If a cluster is sign, enter its timecourse in the blank proxy
                    %nd denote it with the cluster index
                    if PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}(i_elec,i_cluster) < param.pval_plotting
                        PEEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}...
                            (i_elec, find(CurrentEffect.clusterstat...
                            {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse ...
                            == i_cluster)) = i_cluster;
                    end
                end
            end
            
            %Data aggregation over subjects
            %Array to differentiate subject entries
            PEEffect{i_inputdatatyp}.all_subs.sub_index = ...
                [PEEffect{i_inputdatatyp}.all_subs.sub_index; ...
                ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];
            %Array to differentiate TD entries for each electrode
            for i_elec = 1:nSensors_all
                PEEffect{i_inputdatatyp}.all_subs.TD_index    = ...
                    [PEEffect{i_inputdatatyp}.all_subs.TD_index; ...
                    FuncInput_ToneDur_text{i_TD} 's'];
            end
            
            %Clusterinfo aggregation over subjects (Dim: elec * cluster)
            if i_sub == 1 && i_TD == 1
                PEEffect{i_inputdatatyp}.all_subs.clusterstat     = ...
                    PEEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD};
                PEEffect{i_inputdatatyp}.all_subs.pval            = ...
                    PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD};
                PEEffect{i_inputdatatyp}.all_subs.pval_derivative = ...
                    PEEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD};
                PEEffect{i_inputdatatyp}.all_subs.effect_onset    = ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD};
                PEEffect{i_inputdatatyp}.all_subs.effect_offset   = ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD};
                PEEffect{i_inputdatatyp}.all_subs.effect_duration = ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD};
            else
                PEEffect{i_inputdatatyp}.all_subs.clusterstat     = ...
                    [PEEffect{i_inputdatatyp}.all_subs.clusterstat; ...
                    PEEffect{i_inputdatatyp}.per_sub.clusterstat{i_sub, i_TD}];
                PEEffect{i_inputdatatyp}.all_subs.pval            = ...
                    [PEEffect{i_inputdatatyp}.all_subs.pval; ...
                    PEEffect{i_inputdatatyp}.per_sub.pval{i_sub, i_TD}];
                PEEffect{i_inputdatatyp}.all_subs.pval_derivative = ...
                    [PEEffect{i_inputdatatyp}.all_subs.pval_derivative; ...
                    PEEffect{i_inputdatatyp}.per_sub.pval_derivative{i_sub, i_TD}];
                PEEffect{i_inputdatatyp}.all_subs.effect_onset    = ...
                    [PEEffect{i_inputdatatyp}.all_subs.effect_onset; ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_onset{i_sub, i_TD}];
                PEEffect{i_inputdatatyp}.all_subs.effect_offset   = ...
                    [PEEffect{i_inputdatatyp}.all_subs.effect_offset; ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_offset{i_sub, i_TD}];
                PEEffect{i_inputdatatyp}.all_subs.effect_duration = ...
                    [PEEffect{i_inputdatatyp}.all_subs.effect_duration; ...
                    PEEffect{i_inputdatatyp}.per_sub.effect_duration{i_sub, i_TD}];
            end
            
            %Cleanup
            usedElecs_chanposIndex = [];
            clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
        end
        disp([' -- Prediction Error electrodes processed for ' sub ' (both TW, ' FuncInput_DataType{i_inputdatatyp} ') --'])
    end
end
disp(['-- Prediction Error electrodes read out for all input data types, subs, and TW after ' num2str(round(toc/60),2) ' minutes --'])


%% 3) Determine sign. electrodes for each effect
elec_signindex   = struct;

%3.1 Prediction effects
temp        = struct;
temp.label          = [];
temp.inputsignal    = [];
temp.TD_signeffect  = [];
temp.anatcat        = [];
temp.detailedlabel  = [];
temp.coords_elec    = [];
temp.pval           = [];
temp.index_allsubs  = [];
temp.index_ssub     = [];

for i_inputdatatyp = 1:length(FuncInput_DataType)
    %Determine significance filters and cluster information
    elec_signindex.PredEffect{i_inputdatatyp}.array         = ...
        any(PredictionEffect{i_inputdatatyp}.all_subs.pval < ...
        param.pval_plotting,2);
    %1D filter denoting sign. elecs (independent of number of clusters)
    elec_signindex.PredEffect{i_inputdatatyp}.array = ...
        logical(elec_signindex.PredEffect{i_inputdatatyp}.array);
    elec_signindex.PredEffect{i_inputdatatyp}.index         = ...
        find(elec_signindex.PredEffect{i_inputdatatyp}.array);
    elec_signindex.PredEffect{i_inputdatatyp}.num_elecs     = ...
        length(elec_signindex.PredEffect{i_inputdatatyp}.index);
    elec_signindex.PredEffect{i_inputdatatyp}.num_cluster   = ...
        sum(sum(PredictionEffect{i_inputdatatyp}.all_subs.pval < ...
        param.pval_plotting));
    
    %Determine cluster timing information
    for i_elec = 1:length(elec_signindex.PredEffect{i_inputdatatyp}.index)
        elec_signindex.PredEffect{i_inputdatatyp}.onset_minpcluster(i_elec,1) = ...
            PredictionEffect{i_inputdatatyp}.all_subs.effect_onset(...
            elec_signindex.PredEffect{i_inputdatatyp}.index(i_elec),1);
        elec_signindex.PredEffect{i_inputdatatyp}.offset_minpcluster(i_elec,1) = ...
            PredictionEffect{i_inputdatatyp}.all_subs.effect_offset(...
            elec_signindex.PredEffect{i_inputdatatyp}.index(i_elec),1);
        elec_signindex.PredEffect{i_inputdatatyp}.duration_minpcluster(i_elec,1) = ...
            PredictionEffect{i_inputdatatyp}.all_subs.effect_duration(...
            elec_signindex.PredEffect{i_inputdatatyp}.index(i_elec),1);
    end
    
    %Determine labels of sign. elecs and add number to subtitle
    elec_signindex.PredEffect{i_inputdatatyp}.labels        = ...
        PredictionEffect{i_inputdatatyp}.all_subs.label_elec(...
        elec_signindex.PredEffect{i_inputdatatyp}.index);
    elec_signindex.PredEffect{i_inputdatatyp}.anatlabels    = ...
        PredictionEffect{i_inputdatatyp}.all_subs.label_anat(...
        elec_signindex.PredEffect{i_inputdatatyp}.index);
    elec_signindex.PredEffect{i_inputdatatyp}.fulllabels    = ...
        [];
    for i_elec = 1:length(elec_signindex.PredEffect{i_inputdatatyp}.labels)
        elec_signindex.PredEffect{i_inputdatatyp}.fulllabels{i_elec,1} = ...
            [elec_signindex.PredEffect{i_inputdatatyp}.labels{i_elec}, ' ', ...
            PredictionEffect{i_inputdatatyp}.all_subs.TD_index(...
            elec_signindex.PredEffect{i_inputdatatyp}.index(i_elec),:), ' ', ...
            elec_signindex.PredEffect{i_inputdatatyp}.anatlabels{i_elec}];
    end
    
    %Aggregate data across input signals
    temp.label = ...
        [temp.label; split(elec_signindex.PredEffect{i_inputdatatyp}.labels, ' ')];
    temp2_inputsignal = {};
    for i_elecs = 1:length(elec_signindex.PredEffect{i_inputdatatyp}.labels)
        temp2_inputsignal{i_elecs,1} = ...
            FuncInput_DataType{i_inputdatatyp};
    end
    temp.inputsignal = ...
        [temp.inputsignal; temp2_inputsignal];
    temp.TD_signeffect = ...
        [temp.TD_signeffect; cellstr(PredictionEffect{i_inputdatatyp}.all_subs.TD_index(elec_signindex.PredEffect{i_inputdatatyp}.array,:))];
    temp.anatcat = ...
        [temp.anatcat; PredictionEffect{i_inputdatatyp}.all_subs.label_AnatCat(elec_signindex.PredEffect{i_inputdatatyp}.array)];
    temp.detailedlabel = ...
        [temp.detailedlabel; elec_signindex.PredEffect{i_inputdatatyp}.fulllabels];
    temp.coords_elec    = ...
        [temp.coords_elec; PredictionEffect{i_inputdatatyp}.all_subs.coords_elec(elec_signindex.PredEffect{i_inputdatatyp}.array,:)];
    temp.pval = ...
        [temp.pval; PredictionEffect{i_inputdatatyp}.all_subs.pval(elec_signindex.PredEffect{i_inputdatatyp}.array,1:3)];
    temp.index_allsubs = ...
        [temp.index_allsubs; elec_signindex.PredEffect{i_inputdatatyp}.index];
    
end

temp.cluster_timecourse = nan(length(temp.label), 205);
for i_signelec = 1:length(temp.label)
    i_sub = ...
        find(strcmp(subs, temp.label{i_signelec,2}));
    i_TD = ...
        find(strncmp(FuncInput_ToneDur_text, temp.TD_signeffect{i_signelec},3));
    i_inputdatatyp = ...
        find(strncmp(FuncInput_DataType, temp.inputsignal{i_signelec},3));
    i_elec_ssub = ...
        find(strcmp(PredictionEffect{1}.per_sub.elec_labels{i_sub}, temp.label{i_signelec,1}));
    temp.cluster_timecourse(i_signelec, 1:length(PredictionEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec_ssub,:))) = ...
        PredictionEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec_ssub,:);
    
    temp.index_ssub = ...
        [temp.index_ssub; i_elec_ssub];
end

%Summarize sign. electrodes in table
varnames = ...
    {'Electrode Label', 'Subject Label', 'Input Signal', 'TD of effect', ...
    'AnatCat Label', 'Detailed Label', 'Elec Coords', ...
    'Effect p-value', 'Elec Index (all elecs)', 'Elec Index (Ssub elecs)', ...
    'Cluser Timecourse'};
SelElecs.PredEffect = table(...
    temp.label(:,1), ...
    temp.label(:,2), ...
    temp.inputsignal, ...
    temp.TD_signeffect, ...
    temp.anatcat, ...
    temp.detailedlabel, ...
    temp.coords_elec, ...
    temp.pval, ...
    temp.index_allsubs, ...
    temp.index_ssub, ...
    temp.cluster_timecourse, ...
    'VariableNames', varnames);
SelElecs.PredEffect.Properties.Description = ...
    ['Information for electrodes showing a sign. prediction effects selected at p<' num2str(param.pval_plotting) ' (' FDR_label ')'];

%Delete double entries due to multiple effects (TD, input signal)
temp_filter_double = false(height(SelElecs.PredEffect), height(SelElecs.PredEffect));
for i_elec = 1:height(SelElecs.PredEffect)
    temp_filter_double(:,i_elec) = ... %find double entries with same elec and sub
        strcmp(SelElecs.PredEffect{i_elec,1}, SelElecs.PredEffect{:,1}) & ...
        strcmp(SelElecs.PredEffect{i_elec,2}, SelElecs.PredEffect{:,2});
    temp_filter_double(i_elec, i_elec) = false; %ensure that its not current elec
    %For each double entry, compare input signal & TD. If they are different, append label to show double info
    temp_index_double = find(temp_filter_double(:,i_elec));
    for i_doubleelec = 1:length(temp_index_double)
        if strcmp(SelElecs.PredEffect{i_elec,3}, SelElecs.PredEffect{temp_index_double(i_doubleelec),3}) == false
            if contains(char(SelElecs.PredEffect{i_elec,3}), char(SelElecs.PredEffect{temp_index_double(i_doubleelec),3})) == false
                SelElecs.PredEffect{i_elec,3} = ...
                    strcat(SelElecs.PredEffect{i_elec,3}, '/', ...
                    SelElecs.PredEffect{temp_index_double(i_doubleelec),3});
            end
        end
        if strcmp(SelElecs.PredEffect{i_elec,4}, SelElecs.PredEffect{temp_index_double(i_doubleelec),4}) == false
            if contains(char(SelElecs.PredEffect{i_elec,4}), char(SelElecs.PredEffect{temp_index_double(i_doubleelec),4})) == false
                SelElecs.PredEffect{i_elec,4} = ...
                    strcat(SelElecs.PredEffect{i_elec,4}, '/', ...
                    SelElecs.PredEffect{temp_index_double(i_doubleelec),4});
            end
        end
    end
end
temp_filter_double = any(temp_filter_double,2);
%find 1st entry of double elec to keep
temp_index_double = find(temp_filter_double);
temp_list_doubleeleclabel = [];
for i_doubleentries = 1:sum(temp_filter_double)
    if strcmp(temp_list_doubleeleclabel, SelElecs.PredEffect{temp_index_double(i_doubleentries),1}) == false
        temp_list_doubleeleclabel = ...
            [temp_list_doubleeleclabel, SelElecs.PredEffect{temp_index_double(i_doubleentries),1}];
        temp_filter_double(temp_index_double(i_doubleentries)) = false;
    end
end
SelElecs.PredEffect(temp_filter_double,:) = [];

clear temp*

%3.2 Prediction Error effects
temp        = struct;
temp.label          = [];
temp.inputsignal    = [];
temp.TD_signeffect  = [];
temp.anatcat        = [];
temp.detailedlabel  = [];
temp.coords_elec    = [];
temp.pval           = [];
temp.index_allsubs  = [];
temp.index_ssub     = [];

for i_inputdatatyp = 1:length(FuncInput_DataType)
    
    elec_signindex.PEEffect{i_inputdatatyp}.array         = ...
        any(PEEffect{i_inputdatatyp}.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
    elec_signindex.PEEffect{i_inputdatatyp}.index         = ...
        find(elec_signindex.PEEffect{i_inputdatatyp}.array);
    elec_signindex.PEEffect{i_inputdatatyp}.num_elecs     = ...
        length(elec_signindex.PEEffect{i_inputdatatyp}.index);
    elec_signindex.PEEffect{i_inputdatatyp}.num_cluster   = ...
        sum(sum(PEEffect{i_inputdatatyp}.all_subs.pval < param.pval_plotting));
    
    for i_elec = 1:length(elec_signindex.PEEffect{i_inputdatatyp}.index)
        elec_signindex.PEEffect{i_inputdatatyp}.onset_minpcluster(i_elec,1) = ...
            PEEffect{i_inputdatatyp}.all_subs.effect_onset(...
            elec_signindex.PEEffect{i_inputdatatyp}.index(i_elec),1);
        elec_signindex.PEEffect{i_inputdatatyp}.offset_minpcluster(i_elec,1) = ...
            PEEffect{i_inputdatatyp}.all_subs.effect_offset(...
            elec_signindex.PEEffect{i_inputdatatyp}.index(i_elec),1);
        elec_signindex.PEEffect{i_inputdatatyp}.duration_minpcluster(i_elec,1) = ...
            PEEffect{i_inputdatatyp}.all_subs.effect_duration(...
            elec_signindex.PEEffect{i_inputdatatyp}.index(i_elec),1);
    end
    
    %Determine labels of sign. elecs and add number to subtitle
    elec_signindex.PEEffect{i_inputdatatyp}.labels        = ...
        PredictionEffect{i_inputdatatyp}.all_subs.label_elec(...
        elec_signindex.PEEffect{i_inputdatatyp}.index);
    elec_signindex.PEEffect{i_inputdatatyp}.anatlabels    = ...
        PredictionEffect{i_inputdatatyp}.all_subs.label_anat(...
        elec_signindex.PEEffect{i_inputdatatyp}.index);
    elec_signindex.PEEffect{i_inputdatatyp}.fulllabels    = [];
    for i_elec = 1:length(elec_signindex.PEEffect{i_inputdatatyp}.labels)
        elec_signindex.PEEffect{i_inputdatatyp}.fulllabels{i_elec,1} = ...
            [elec_signindex.PEEffect{i_inputdatatyp}.labels{i_elec}, ' ', ...
            PEEffect{i_inputdatatyp}.all_subs.TD_index(...
            elec_signindex.PEEffect{i_inputdatatyp}.index(i_elec),:), ' ', ...
            elec_signindex.PEEffect{i_inputdatatyp}.anatlabels{i_elec}];
    end
    
    %Aggregate data across input signals
    temp.label = ...
        [temp.label; split(elec_signindex.PEEffect{i_inputdatatyp}.labels, ' ')];
    temp2_inputsignal = {};
    for i_elecs = 1:length(elec_signindex.PEEffect{i_inputdatatyp}.labels)
        temp2_inputsignal{i_elecs,1} = ...
            FuncInput_DataType{i_inputdatatyp};
    end
    temp.inputsignal = ...
        [temp.inputsignal; temp2_inputsignal];
    temp.TD_signeffect = ...
        [temp.TD_signeffect; cellstr(PEEffect{i_inputdatatyp}.all_subs.TD_index(elec_signindex.PEEffect{i_inputdatatyp}.array,:))];
    temp.anatcat = ...
        [temp.anatcat; PEEffect{i_inputdatatyp}.all_subs.label_AnatCat(elec_signindex.PEEffect{i_inputdatatyp}.array)];
    temp.detailedlabel = ...
        [temp.detailedlabel; elec_signindex.PEEffect{i_inputdatatyp}.fulllabels];
    temp.coords_elec    = ...
        [temp.coords_elec; PredictionEffect{i_inputdatatyp}.all_subs.coords_elec(elec_signindex.PEEffect{i_inputdatatyp}.array,:)];
    temp.pval = ...
        [temp.pval; PEEffect{i_inputdatatyp}.all_subs.pval(elec_signindex.PEEffect{i_inputdatatyp}.array,1:3)];
    temp.index_allsubs = ...
        [temp.index_allsubs; elec_signindex.PEEffect{i_inputdatatyp}.index];
    
end

temp.cluster_timecourse = nan(length(temp.label), 205);
for i_signelec = 1:length(temp.label)
    i_sub = ...
        find(strcmp(subs, temp.label{i_signelec,2}));
    i_TD = ...
        find(strncmp(FuncInput_ToneDur_text, temp.TD_signeffect{i_signelec},3));
    i_inputdatatyp = ...
        find(strncmp(FuncInput_DataType, temp.inputsignal{i_signelec},3));
    i_elec_ssub = ...
        find(strcmp(PEEffect{1}.per_sub.elec_labels{i_sub}, temp.label{i_signelec,1}));
    temp.cluster_timecourse(i_signelec, 1:length(PEEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec_ssub,:))) = ...
        PEEffect{i_inputdatatyp}.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec_ssub,:);
    
    temp.index_ssub = ...
        [temp.index_ssub; i_elec_ssub];
end

%Summarize sign. electrodes in table
varnames = ...
    {'Electrode Label', 'Subject Label', 'Input Signal', 'TD of effect', ...
    'AnatCat Label', 'Detailed Label', 'Elec Coords', ...
    'Effect p-value', 'Elec Index (all elecs)', 'Elec Index (Ssub elecs)', ...
    'Cluser Timecourse'};
SelElecs.PEEffect = table(...
    temp.label(:,1), ...
    temp.label(:,2), ...
    temp.inputsignal, ...
    temp.TD_signeffect, ...
    temp.anatcat, ...
    temp.detailedlabel, ...
    temp.coords_elec, ...
    temp.pval, ...
    temp.index_allsubs, ...
    temp.index_ssub, ...
    temp.cluster_timecourse, ...
    'VariableNames', varnames);
SelElecs.PEEffect.Properties.Description = ...
    ['Information for electrodes showing a sign. complex prediction error effects selected at p<' num2str(param.pval_plotting) ' (' FDR_label ')'];

%Delete double entries due to multiple effects (TD, input signal)
temp_filter_double = false(height(SelElecs.PEEffect), height(SelElecs.PEEffect));
for i_elec = 1:height(SelElecs.PEEffect)
    temp_filter_double(:,i_elec) = ... %find double entries with same elc and sub
        strcmp(SelElecs.PEEffect{i_elec,1}, SelElecs.PEEffect{:,1}) & ...
        strcmp(SelElecs.PEEffect{i_elec,2}, SelElecs.PEEffect{:,2});
    temp_filter_double(i_elec, i_elec) = false; %ensure that its not current elec
    %For each double entry, compare input signal & TD. If they are different, append label to show double info
    temp_index_double = find(temp_filter_double(:,i_elec));
    for i_doubleelec = 1:length(temp_index_double)
        if strcmp(SelElecs.PEEffect{i_elec,3}, SelElecs.PEEffect{temp_index_double(i_doubleelec),3}) == false
            if contains(char(SelElecs.PEEffect{i_elec,3}), char(SelElecs.PEEffect{temp_index_double(i_doubleelec),3})) == false
                SelElecs.PEEffect{i_elec,3} = ...
                    strcat(SelElecs.PEEffect{i_elec,3}, '/', ...
                    SelElecs.PEEffect{temp_index_double(i_doubleelec),3});
            end
        end
        if strcmp(SelElecs.PEEffect{i_elec,4}, SelElecs.PEEffect{temp_index_double(i_doubleelec),4}) == false
            if contains(char(SelElecs.PEEffect{i_elec,4}), char(SelElecs.PEEffect{temp_index_double(i_doubleelec),4})) == false
                SelElecs.PEEffect{i_elec,4} = ...
                    strcat(SelElecs.PEEffect{i_elec,4}, '/', ...
                    SelElecs.PEEffect{temp_index_double(i_doubleelec),4});
            end
        end
    end
end
temp_filter_double = any(temp_filter_double,2);
%find 1st entry of double elec to keep
temp_index_double = find(temp_filter_double);
temp_list_doubleeleclabel = [];
for i_doubleentries = 1:sum(temp_filter_double)
    if strcmp(temp_list_doubleeleclabel, SelElecs.PEEffect{temp_index_double(i_doubleentries),1}) == false
        temp_list_doubleeleclabel = ...
            [temp_list_doubleeleclabel, SelElecs.PEEffect{temp_index_double(i_doubleentries),1}];
        temp_filter_double(temp_index_double(i_doubleentries)) = false;
    end
end
SelElecs.PEEffect(temp_filter_double,:) = [];

clear temp*


%3.3 Determine elecs that show sign effects for both Pred and PE effect (independent of TD, signal)
temp        = struct;
temp.label          = [];
temp.sub          = [];
temp.inputsignal    = [];
temp.TD_signeffect  = [];
temp.anatcat        = [];
temp.detailedlabel  = [];
temp.coords_elec    = [];
temp.index_ssub     = [];

for i_inputdatatyp = 1:length(FuncInput_DataType)
    for i_Predelec = 1:height(SelElecs.PredEffect)
        %Compare pred vs. PE elec via label & sub
        filt_elec_botheffects = [];
        filt_elec_botheffects(:,1) = ...
            strcmp(SelElecs.PredEffect{i_Predelec,1}, SelElecs.PEEffect{:,1});
        filt_elec_botheffects(:,2) = ...
            strcmp(SelElecs.PredEffect{i_Predelec,2}, SelElecs.PEEffect{:,2});
        filt_elec_botheffects = all(filt_elec_botheffects,2);
        
        if any(filt_elec_botheffects)
            %Select 1st common elec (for stable parameters like coords)
            entry1 = find(filt_elec_botheffects);
            entry1 = entry1(1);
            
            %Read out variable components (TD, Signal)
            temp.label = [temp.label; SelElecs.PredEffect{i_Predelec,1}];
            temp.sub = [temp.sub; SelElecs.PredEffect{i_Predelec,2}];
            temp.inputsignal = [temp.inputsignal; ...
                strcat(SelElecs.PredEffect{i_Predelec,3}, '-', ...
                strjoin(SelElecs.PEEffect{filt_elec_botheffects,3}))];
            
            temp.TD_signeffect = [temp.TD_signeffect; ...
                strcat(SelElecs.PredEffect{i_Predelec,4}, '-', ...
                strjoin(SelElecs.PEEffect{filt_elec_botheffects,4}))];
            
            temp.anatcat = [temp.anatcat; SelElecs.PredEffect{i_Predelec,5}];
            temp.detailedlabel = [temp.detailedlabel; SelElecs.PredEffect{i_Predelec,6}];
            temp.coords_elec = [temp.coords_elec; SelElecs.PredEffect{i_Predelec,7}];
            temp.index_ssub     = [temp.index_ssub; SelElecs.PredEffect{i_Predelec,10}];
            
        end
    end
end

varnames = ...
    {'Electrode Label', 'Subject Label', 'Input Signal', 'TD of effect', ...
    'AnatCat Label', 'Detailed Label', 'Elec Coords', ...
    'Elec Index (Ssub elecs)'};
SelElecs.PredPEEffect = table(...
    temp.label, ...
    temp.sub, ...
    temp.inputsignal, ...
    temp.TD_signeffect, ...
    temp.anatcat, ...
    temp.detailedlabel, ...
    temp.coords_elec, ...
    temp.index_ssub, ...
    'VariableNames', varnames);
SelElecs.PredPEEffect.Properties.Description = ...
    ['Information for electrodes showing both a sign. prediciton and complex prediction error effects selected at p<' num2str(param.pval_plotting) ' (' FDR_label ')'];

clear temp*

%Delete double entries due to multiple effects (TD, input signal)
temp_filter_double = false(height(SelElecs.PredPEEffect), height(SelElecs.PredPEEffect));
for i_elec = 1:height(SelElecs.PredPEEffect)
    temp_filter_double(:,i_elec) = ... %find double entries with same elc and sub
        strcmp(SelElecs.PredPEEffect{i_elec,1}, SelElecs.PredPEEffect{:,1}) & ...
        strcmp(SelElecs.PredPEEffect{i_elec,2}, SelElecs.PredPEEffect{:,2});
    temp_filter_double(i_elec, i_elec) = false; %ensure that its not current elec
    %For each double entry, compare input signal & TD. If they are different, append label to show double info
    temp_index_double = find(temp_filter_double(:,i_elec));
    for i_doubleelec = 1:length(temp_index_double)
        if strcmp(SelElecs.PredPEEffect{i_elec,3}, SelElecs.PredPEEffect{temp_index_double(i_doubleelec),3}) == false
            if contains(char(SelElecs.PredPEEffect{i_elec,3}), char(SelElecs.PredPEEffect{temp_index_double(i_doubleelec),3})) == false
                SelElecs.PredPEEffect{i_elec,3} = ...
                    strcat(SelElecs.PredPEEffect{i_elec,3}, '/', ...
                    SelElecs.PredPEEffect{temp_index_double(i_doubleelec),3});
            end
        end
        if strcmp(SelElecs.PredPEEffect{i_elec,4}, SelElecs.PredPEEffect{temp_index_double(i_doubleelec),4}) == false
            if contains(char(SelElecs.PredPEEffect{i_elec,4}), char(SelElecs.PredPEEffect{temp_index_double(i_doubleelec),4})) == false
                SelElecs.PredPEEffect{i_elec,4} = ...
                    strcat(SelElecs.PredPEEffect{i_elec,4}, '/', ...
                    SelElecs.PredPEEffect{temp_index_double(i_doubleelec),4});
            end
        end
    end
end
temp_filter_double = any(temp_filter_double,2);
%find 1st entry of double elec to keep
temp_index_double = find(temp_filter_double);
temp_list_doubleeleclabel = [];
for i_doubleentries = 1:sum(temp_filter_double)
    if strcmp(temp_list_doubleeleclabel, SelElecs.PredPEEffect{temp_index_double(i_doubleentries),1}) == false
        temp_list_doubleeleclabel = ...
            [temp_list_doubleeleclabel, SelElecs.PredPEEffect{temp_index_double(i_doubleentries),1}];
        temp_filter_double(temp_index_double(i_doubleentries)) = false;
    end
end
SelElecs.PredPEEffect(temp_filter_double,:) = [];

clear temp*


%% 4) Plot pairings between sign. Pred - Pred  electrodes
%Plot electrode pairings between electrodes with sign effects (across subjects, TD),
%seperately for each prediciton-effect-type.

SelElecs.Pairs = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_poststepFigs == 1
    %Frontal - temporal
    SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    %Frontal - Parietal
    SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in parietal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %Parietal - Temporal
    SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'SupParLob') %If current electrode is Parietal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %% Plot figure showing electrodes and connections between all anatomical
    %(across subjects)
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [1 3];
    sgtitle(['All subs (n = ' num2str(length(subs)) ') - Pred to Pred (' folderlabel_inputsignal ')'])
    
    effect_comp = 'Pred_Pred';
    cmap  = distinguishable_colors(3);
    clims = [1 3];
    textcolor_rgb   = [0 0 0];
    
    %Frontal - Temporal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            if ~isnan(SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                %Note: the indices refer to the respective effect-specific  across-subject electrode listing
                %(e.g., 5 = 5th electrode in SelElecs.PredEffect)
                if iscolumn(unique(SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'VentralT'))) = 2;
    
    %Note: Problem: elec_coords is mixed (source + target) and ordered,
    %while elec_pairings is not mixed and not ordered. Both var refer to
    %the same elecs and thus should have the same orientation
    elec_pairings = SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_Pred.Frontal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 1, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Frontal - Parietal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            if ~isnan(SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                if iscolumn(unique(SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'SupParLob'))) = 3;

    elec_pairings = SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
    
    sp_title = {['{\color{blue}Frontal} to {\color{green}Parietal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_Pred.Frontal_Parietal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 2, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Parietal - Temporal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            if ~isnan(SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                if iscolumn(unique(SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'SupParLob'))) = 3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PredEffect{curr_elec_ind,5}, 'VentralT'))) = 2;

    elec_pairings = SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
    
    sp_title = {['{\color{green}Parietal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_Pred.Parietal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 3, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
end

filename     = ['Surf1H_' ...
    'Allsubn' num2str(length(subs)) '_' ...
    'PairsSignElec_PredPred_'  ...
    'AllTD_' folderlabel_inputsignal '.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version
close;


%% 5) Plot pairings between sign. PE - PE  electrodes
%Plot electrode pairings between electrodes with sign effects (across subjects, TD),
%seperately for each prediciton-effect-type.

if plot_poststepFigs == 1
    %Frontal - temporal
    SelElecs.Pairs.PE_PE.Frontal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_PE.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_PE.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ... %Same subject
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'AntPFC'))); %and specific location
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    %Frontal - Parietal
    SelElecs.Pairs.PE_PE.Frontal_Parietal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in parietal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_PE.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_PE.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %Parietal - Temporal
    SelElecs.Pairs.PE_PE.Parietal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'SupParLob') %If current electrode is Parietal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_PE.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_PE.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %% Plot figure showing electrodes and connections between all anatomical
    %(across subjects)
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [1 3];
    sgtitle(['All subs (n = ' num2str(length(subs)) ') - PE to PE (' folderlabel_inputsignal ')'])
    
    effect_comp = 'PE_PE';
    cmap  = distinguishable_colors(3);
    clims = [1 3];
    textcolor_rgb   = [0 0 0];
    
    %Frontal - Temporal
    %To do: Copy workaround to other plots and check result
    curr_elec_ind = [];
    for i_sub = 1:length(subs) 
        if any(strcmp(subs(i_sub),SelElecs.PEEffect{:,2}))
            if ~isnan(SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                if iscolumn(unique(SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    curr_elec_ind = sort(curr_elec_ind); %Important to subsequently order elecs
    
    elec_coords = SelElecs.PEEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'VentralT'))) = 2;

    elec_pairings = SelElecs.Pairs.PE_PE.Frontal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
       
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.PE_PE.Frontal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 1, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Frontal - Parietal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PEEffect{:,2}))
            if ~isnan(SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                if iscolumn(unique(SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    elec_coords = SelElecs.PEEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'AntPFC'))) = 1;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'SupParLob'))) = 3;

    elec_pairings = SelElecs.Pairs.PE_PE.Frontal_Parietal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{green}Parietal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.PE_PE.Frontal_Parietal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 2, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Parietal - Temporal
    curr_elec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PEEffect{:,2}))
            if ~isnan(SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                if iscolumn(unique(SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}))
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})];
                else %Workaround because output is row vector when only 1 pair found
                    curr_elec_ind = [curr_elec_ind; ...
                        unique(SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})'];
                end
            end
        end
    end
    elec_coords = SelElecs.PEEffect{curr_elec_ind,7};
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    for i_elec = 1:size(elec_coords,1) %No electrode labels
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_elec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = ones(size(elec_coords,1),1)*3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'SupParLob'))) = 3;
    elec_typeindex(~cellfun('isempty',strfind(SelElecs.PEEffect{curr_elec_ind,5}, 'VentralT'))) = 2;

    elec_pairings = SelElecs.Pairs.PE_PE.Parietal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
        for i_column = 1:size(elec_pairings,2)
            elec_pairings_ordered(i_row,i_column) = ...
                find(curr_elec_ind == elec_pairings(i_row,i_column));
        end
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{green}Parietal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.PE_PE.Parietal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 3, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
end

filename     = ['Surf1H_' ...
    'Allsubn' num2str(length(subs)) '_' ...
    'PairsSignElec_PEPE_'  ...
    'AllTD_' folderlabel_inputsignal '.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version
close;


%% 5) Plot pairings between sign. Pred - PE  electrodes

if plot_poststepFigs == 1
    %Frontal - temporal
    SelElecs.Pairs.Pred_PE.Frontal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_PE.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_PE.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    %Frontal - Parietal
    SelElecs.Pairs.Pred_PE.Frontal_Parietal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in parietal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_PE.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_PE.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %Parietal - Temporal
    SelElecs.Pairs.Pred_PE.Parietal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PredEffect)
        %Determine current elec information
        temp_label  = SelElecs.PredEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PredEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PredEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'SupParLob') %If current electrode is Parietal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PEEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.Pred_PE.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PEEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.Pred_PE.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %% Plot figure showing electrodes and connections between all anatomical
    %(across subjects)
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [1 3];
    sgtitle(['All subs (n = ' num2str(length(subs)) ') - Pred to PE (' folderlabel_inputsignal ')'])
    
    effect_comp = 'Pred_PE';
    cmap  = distinguishable_colors(3);
    clims = [1 3];
    textcolor_rgb   = [0 0 0];
    
    %Frontal - Temporal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_Predelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PEEffect{curr_PEelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_Predelec_ind),1); ones(length(curr_PEelec_ind),1)*2];

    elec_pairings = SelElecs.Pairs.Pred_PE.Frontal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,2));
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred to PE elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_PE.Frontal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 1, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Frontal - Parietal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_Predelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PEEffect{curr_PEelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_Predelec_ind),1); ones(length(curr_PEelec_ind),1)*3];

    elec_pairings = SelElecs.Pairs.Pred_PE.Frontal_Parietal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,2));
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{green}Parietal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred to PE elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_PE.Frontal_Parietal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 2, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Parietal - Temporal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PredEffect{curr_Predelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PEEffect{curr_PEelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_Predelec_ind),1)*3; ones(length(curr_PEelec_ind),1)*2];

    elec_pairings = SelElecs.Pairs.Pred_PE.Parietal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,2));
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{green}Parietal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' Pred to PE elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_PE.Parietal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 3, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
end

filename     = ['Surf1H_' ...
    'Allsubn' num2str(length(subs)) '_' ...
    'PairsSignElec_PredPE_'  ...
    'AllTD_' folderlabel_inputsignal '.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version
close;


%% 5) Plot pairings between sign. PE - Pred  electrodes

if plot_poststepFigs == 1
    %Frontal - temporal
    SelElecs.Pairs.PE_Pred.Frontal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_Pred.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_Pred.Frontal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    %Frontal - Parietal
    SelElecs.Pairs.PE_Pred.Frontal_Parietal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'AntPFC') %If current electrode is frontal
            
            temp_index_otherelecs = ...%Find other sign. elecs in parietal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_Pred.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_Pred.Frontal_Parietal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'AntPFC')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'SupParLob')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %Parietal - Temporal
    SelElecs.Pairs.PE_Pred.Parietal_Temporal.Pairs_perelec = struct;
    for i_elec = 1:height(SelElecs.PEEffect)
        %Determine current elec information
        temp_label  = SelElecs.PEEffect{i_elec,1}; %Elec label
        temp_sub    = SelElecs.PEEffect{i_elec,2}; %Sub
        temp_anatcat = SelElecs.PEEffect{i_elec,5}; %Anatomical location
        if strfind(temp_anatcat{1}, 'SupParLob') %If current electrode is Parietal
            
            temp_index_otherelecs = ...%Find other sign. elecs in temporal areas
                find(...
                ~strcmp(temp_label, SelElecs.PredEffect{:,1}) & ...
                ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
            if ~isempty(strfind(temp_anatcat{1}, ',')) %Create label for subfield
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' extractBefore(temp_anatcat{1}, ',')];
            else
                temp_fieldname = ...
                    [temp_label{1} '_' temp_sub{1} '_' temp_anatcat{1}];
            end
            
            if ~isempty(temp_index_otherelecs) %Store pair information in new struct
                SelElecs.Pairs.PE_Pred.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = ...
                    SelElecs.PredEffect(temp_index_otherelecs,:);
            else
                SelElecs.Pairs.PE_Pred.Parietal_Temporal.Pairs_perelec.(temp_fieldname) = [];
            end
        end
    end
    %Determine number of possible pairings for each subject
    SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings = [];
    for i_sub = 1:length(subs)
        sub = subs(i_sub);
        sign_elecs_currsub1 = ...
            find(strcmp(sub, SelElecs.PEEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PEEffect{:,5}, 'SupParLob')));
        sign_elecs_currsub2 = ...
            find(strcmp(sub, SelElecs.PredEffect{:,2}) & ...
            ~cellfun('isempty',strfind(SelElecs.PredEffect{:,5}, 'VentralT')));
        %Determine number of possible pairings for each subject
        [p,q] = meshgrid(sign_elecs_currsub1, sign_elecs_currsub2);
        temp_pairs1 = [p(:) q(:)];
        if size(temp_pairs1,1) > 1
            %delete same-same pairings
            temp_pairs2 = [];
            for i_pairs = 1:length(temp_pairs1)
                if temp_pairs1(i_pairs,1) ~= temp_pairs1(i_pairs,2)
                    temp_pairs2 = [temp_pairs2; temp_pairs1(i_pairs,:)];
                end
            end
            %delete doubles
            i_doubles = false(length(temp_pairs2),1);
            i_compared = false(length(temp_pairs2),1);
            for i_pairs = 1:length(temp_pairs2)
                temp_doubles = find(all(temp_pairs2(i_pairs,:) == fliplr(temp_pairs2),2));
                i_compared(i_pairs) = true;
                
                if i_compared(temp_doubles) == false
                    i_doubles(temp_doubles) = true;
                end
            end
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                size(temp_pairs2(~i_doubles,:),1);
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs2(~i_doubles,:);
        elseif size(temp_pairs1,1) == 1 && temp_pairs1(1) ~= temp_pairs1(2)
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                1;
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                temp_pairs1;
        else
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Num_UniquePairings_persub(1,i_sub) = ...
                0;
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub} = ...
                [nan nan];
        end
        SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings = ...
            [SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings; ...
            SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}];
    end
    clear temp* sign*
    
    
    %% Plot figure showing electrodes and connections between all anatomical
    %(across subjects)
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    DimSubplot = [1 3];
    sgtitle(['All subs (n = ' num2str(length(subs)) ') - PE to Pred (' folderlabel_inputsignal ')'])
    
    effect_comp = 'PE_Pred';
    cmap  = distinguishable_colors(3);
    clims = [1 3];
    textcolor_rgb   = [0 0 0];
    
    %Frontal - Temporal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PEEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PEEffect{curr_PEelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PredEffect{curr_Predelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_PEelec_ind),1); ones(length(curr_Predelec_ind),1)*2];

    elec_pairings = SelElecs.Pairs.PE_Pred.Frontal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,2));
    end

    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' PE to Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.PE_Pred.Frontal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 1, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Frontal - Parietal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub})
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PEEffect{curr_PEelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PredEffect{curr_Predelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_PEelec_ind),1); ones(length(curr_Predelec_ind),1)*3];
 
    elec_pairings = SelElecs.Pairs.PE_Pred.Frontal_Parietal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,2));
    end

    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{blue}Frontal} to {\color{green}Parietal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' PE to Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.PE_Pred.Frontal_Parietal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 2, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
    %Parietal - Temporal
    curr_Predelec_ind = [];
    curr_PEelec_ind = [];
    for i_sub = 1:length(subs)
        if any(strcmp(subs(i_sub),SelElecs.PredEffect{:,2}))
            %Determine pairings or, if not present, sign. elecs
            if ~isnan(SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub})
                curr_PEelec_ind = [curr_PEelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}(:,1))];
                curr_Predelec_ind = [curr_Predelec_ind; ...
                    unique(SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings_persub{i_sub}(:,2))];
            end
        end
    end
    elec_coords = SelElecs.PEEffect{curr_PEelec_ind,7};
    elec_coords = [elec_coords; SelElecs.PredEffect{curr_Predelec_ind,7}];
    elec_coords(:,1) = abs(elec_coords(:,1)) * -1;%Project all electrodes on one hemisphere
    elec_labels = {};
    i_entry     = 0;
    for i_elec = 1:length(curr_PEelec_ind)
        elec_labels{i_elec} = cell2mat(SelElecs.PEEffect{curr_PEelec_ind(i_elec),1});
        i_entry = i_entry + 1;
    end
    for i_elec = 1:length(curr_Predelec_ind)
        elec_labels{i_elec + i_entry} = cell2mat(SelElecs.PredEffect{curr_Predelec_ind(i_elec),1});
    end
    elec_data = ones(size(elec_coords,1),1);
    elec_signindex = true(size(elec_coords,1),1);
    elec_typeindex = [ones(length(curr_PEelec_ind),1)*3; ones(length(curr_Predelec_ind),1)*2];

    elec_pairings = SelElecs.Pairs.PE_Pred.Parietal_Temporal.Ind_UniquePairings;
    elec_pairings = elec_pairings(~isnan(elec_pairings(:,1)),:);
    elec_pairings_ordered = zeros(size(elec_pairings));
    for i_row = 1:size(elec_pairings,1)
            elec_pairings_ordered(i_row,1) = ...
                find(curr_PEelec_ind == elec_pairings(i_row,1));
            elec_pairings_ordered(i_row,2) = ...
                find(curr_Predelec_ind == elec_pairings(i_row,2));
    end
    
    chanSize = ones(1,length(elec_data))*3; %electrode size (arbitrary)
    sp_title = {['{\color{green}Parietal} to {\color{red}Temporal}  {\color{black}electrodes} '], ...
        [ num2str(size(elec_coords,1)) ' PE to Pred elecs - ' ...
        num2str(sum(SelElecs.Pairs.Pred_PE.Parietal_Temporal.Num_UniquePairings_persub)) ...
        ' elec pairs']};
    
    if ~isempty(elec_data)
        NASTD_ECoG_Plot_SubplotSignElecsSurf_ColorAnat_LH_Pairs...
            (elec_coords, elec_labels, elec_typeindex, elec_signindex, ...
            elec_pairings_ordered, effect_comp, ...
            chanSize, clims, cmap, textcolor_rgb, ...
            DimSubplot, 3, [0 -0.02 0 0], [], sp_title);
    end
    clear elec*
    
end

filename     = ['Surf1H_' ...
    'Allsubn' num2str(length(subs)) '_' ...
    'PairsSignElec_PEPred_'  ...
    'AllTD_' folderlabel_inputsignal '.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version
close;

end