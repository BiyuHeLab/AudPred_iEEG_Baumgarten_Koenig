function NASTD_ECoG_StimCorr_CompareResElecs ...
    (subs, ...
    FuncInput_DataType, FuncInput_ToneDur_text, ...
    param, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Compare the number and location of sign. electrodes determined by
%1) Single tone processing analysis
%2) Identical sequence procesing analysis vs. Similar sequence processing analysis
%3) Prediction effect
%4) History Integration effect
%in order to check how they overlap in terms of number/identity and spatial
%location.

% set(0, 'DefaultFigureVisible', 'on');

%% 0.1) Specify vars, paths, and setup fieldtrip
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

path_fig = ([paths_NASTD_ECoG.ECoGdata_StimCorr ... 
    '/Allsub_n' num2str(length(subs)) '/Figs/CompareResElecs/p' ...
    num2str(param.pval_plotting) '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

%% 1) Load stimulus correlation data and aggregate relevant info across subjects
clear ElecLabels Effect Effect_allsub
ElecLabels_allsub.elec_coords       = [];
ElecLabels_allsub.elec_sublabels    = [];
ElecLabels_allsub.elec_labels       = [];
ElecLabels_allsub.elec_anatlabels   = [];
usedElecs_chanposIndex              = [];
Effect_allsub.sub_index             = [];

for i_sub = 1:length(subs)
    
    tic
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --'])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    %1.1 Load ECoG preproc data for channel labels and position
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
    load(loadfile_ECoGpreprocdata);    
    %Determine basic parameters
    SampleFreq  = DataClean_AllTrials.fsample;
    
    %1.2 Load stimulus correlation data and select current effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
    load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
        FuncInput_ToneDur_text 'sTD.mat'], ...
        'corr_ttest', 'SensorLabels');
    ElecLabels{i_sub}.StimCorr  = SensorLabels;
    Effect{i_sub}.StimCorr_t    = corr_ttest.t;
    Effect{i_sub}.StimCorr_pval = corr_ttest.p;   
    clear corr_ttest SensorLabels
    
    %1.3 Load tone processing data
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
    load([path_inputdata sub '_ToneProc_' FuncInput_DataType '_' ...
        FuncInput_ToneDur_text 'sTD.mat'], ...
        'stat_N100Ampvs0_BLperTone', 'SensorLabels');
    ElecLabels{i_sub}.ToneProc  = SensorLabels;
    Effect{i_sub}.ToneProc_t    = stat_N100Ampvs0_BLperTone.t;
    Effect{i_sub}.ToneProc_pval = stat_N100Ampvs0_BLperTone.p;   
    clear stat_N100Ampvs0_BLperTone SensorLabels
    
    %1.4A Load prediction effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];
    load([path_inputdata sub '_PredEffectsCluster_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'sTD.mat'], ...
      'PredEffect' , 'labels_loadedData');
    ElecLabels{i_sub}.Pred          = labels_loadedData;
    for i_elec = 1:length(PredEffect.clusterstat)
        Effect{i_sub}.Pred_clusterstat(i_elec,1)  = ...
            PredEffect.clusterstat{i_elec}.maxStatSumAbs;
        Effect{i_sub}.Pred_cluster_pval(i_elec,1) = ...
            min(PredEffect.clusterstat{i_elec}.cluster_pval);
    end  
    clear PredEffect labels_loadedData
    %1.4B Load simple prediction error effect
    load([path_inputdata sub '_PredEffectsCluster_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'sTD.mat'], ...
      'SimplePredErrEffect' , 'labels_loadedData');
    ElecLabels{i_sub}.SimplePredErr = labels_loadedData;
    for i_elec = 1:length(SimplePredErrEffect.clusterstat)
        Effect{i_sub}.SimplePredErr_clusterstat(i_elec,1)  = ...
            SimplePredErrEffect.clusterstat{i_elec}.maxStatSumAbs;
        Effect{i_sub}.SimplePredErr_cluster_pval(i_elec,1) = ...
            min(SimplePredErrEffect.clusterstat{i_elec}.cluster_pval);
    end    
    clear SimplePredErrEffect labels_loadedData
    %1.4C Load complex prediction error effect
    load([path_inputdata sub '_PredEffectsCluster_' FuncInput_DataType ...
        '_' FuncInput_ToneDur_text 'sTD.mat'], ...
      'ComplexPredErrEffect' , 'labels_loadedData');
    ElecLabels{i_sub}.ComplexPredErr = labels_loadedData;
    for i_elec = 1:length(ComplexPredErrEffect.clusterstat)
        Effect{i_sub}.ComplexPredErr_clusterstat(i_elec,1)  = ...
            ComplexPredErrEffect.clusterstat{i_elec}.maxStatSumAbs;
        Effect{i_sub}.ComplexPredErr_cluster_pval(i_elec,1) = ...
            min(ComplexPredErrEffect.clusterstat{i_elec}.cluster_pval);
    end
    clear ComplexPredErrEffect labels_loadedData
    
    %1.5 Load SHI effect
    path_inputdata = [paths_NASTD_ECoG.ECoGdata_HisTrack sub ...
        '/ExpvsShuff/50msTW/'];
    SHI_LabelToneDur = [num2str(str2double(FuncInput_ToneDur_text) * 1000) 'msTD'];
    if strcmp(FuncInput_DataType, 'LP35Hz') || strcmp(FuncInput_DataType, 'HighGamma')
        load([path_inputdata sub '_' FuncInput_DataType '_' SHI_LabelToneDur ...
            '_HisTrack_ExpvsShuffKprimeCombFolds_unbalanced.mat']); %filename: ExpKprime_data
        ElecLabels{i_sub}.SHI = labels_loadedData;
        for i_TW = 1:length(Kprime_data.Exp.avgFolds.Kprime_Avg)
            Effect{i_sub}.SHI_Kprime(:,i_TW)    = ...
                Kprime_data.Exp.avgFolds.Kprime_Avg{i_TW}';
            Effect{i_sub}.SHI_pval(:,i_TW)      = ...
                Kprime_data.Exp.avgFolds.pval_ExpvsShuff{i_TW}';
        end
        clear Kprime_data labels_loadedData        
    end

    %Check if all loaded data has the same amount of analyzed electrodes
    nSensors_subfield = [];
    for i_subfields = 1:length(fieldnames(ElecLabels{i_sub}))
        curr_fieldnames = fieldnames(ElecLabels{i_sub});
        nSensors_subfield = [nSensors_subfield ...
            length(ElecLabels{i_sub}.(curr_fieldnames{i_subfields}))];
    end
    if length(unique(nSensors_subfield)) > 1
        disp(['Inconsistent number of electrodes analyzed for ' sub])
        pause
    end
    nSensors = length(ElecLabels{i_sub}.StimCorr);    
    
    %1.6 Read electrode labels, the corresponding coordinates, 
    %and anatomical labels for all analyzed electrodes 
    %and aggregate them across subjects   
    for i_elec = 1:length(ElecLabels{i_sub}.StimCorr)
        usedElecs_chanposIndex(1,i_elec) = ...
            find(strcmp(ElecLabels{i_sub}.StimCorr{i_elec}, ...
            DataClean_AllTrials.elec.label));
    end    
    ElecLabels{i_sub}.elec_coords       = ...
        DataClean_AllTrials.elec.chanpos(usedElecs_chanposIndex,:);    
    ElecLabels_allsub.elec_coords       = ...
        [ElecLabels_allsub.elec_coords; ElecLabels{i_sub}.elec_coords];
    ElecLabels_allsub.elec_anatlabels   = ...
        [ElecLabels_allsub.elec_anatlabels; ...
        DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];    
    for i_elec = 1:length(usedElecs_chanposIndex)
        ElecLabels_allsub.elec_labels{end+1,1}     = ...
            [ElecLabels{i_sub}.StimCorr{i_elec}];
        ElecLabels_allsub.elec_sublabels{end+1,1}  = ...
            [ElecLabels{i_sub}.StimCorr{i_elec} ' ' sub]; 
    end   
    
    %1.7 Store effect results across all subjects
    Effect_allsub.sub_index = ...
        [Effect_allsub.sub_index; ...
        ones(length(ElecLabels{i_sub}.elec_coords),1)*i_sub];
    
    if i_sub == 1
        Effect_allsub.StimCorr_t        = Effect{i_sub}.StimCorr_t;
        Effect_allsub.StimCorr_pval     = Effect{i_sub}.StimCorr_pval;
        Effect_allsub.ToneProc_t        = Effect{i_sub}.ToneProc_t;
        Effect_allsub.ToneProc_pval     = Effect{i_sub}.ToneProc_pval;
        
        Effect_allsub.Pred_clusterstat              = ...
            Effect{i_sub}.Pred_clusterstat;
        Effect_allsub.Pred_cluster_pval             = ...
            Effect{i_sub}.Pred_cluster_pval;
        Effect_allsub.SimplePredErr_clusterstat     = ...
            Effect{i_sub}.SimplePredErr_clusterstat;
        Effect_allsub.SimplePredErr_cluster_pval    = ...
            Effect{i_sub}.SimplePredErr_cluster_pval;
        Effect_allsub.ComplexPredErr_clusterstat    = ...
            Effect{i_sub}.ComplexPredErr_clusterstat;
        Effect_allsub.ComplexPredErr_cluster_pval   = ...
            Effect{i_sub}.ComplexPredErr_cluster_pval;
        
        Effect_allsub.SHI_Kprime        = Effect{i_sub}.SHI_Kprime;
        Effect_allsub.SHI_pval          = Effect{i_sub}.SHI_pval;        
        
    else
        Effect_allsub.StimCorr_t        = ...
            [Effect_allsub.StimCorr_t; Effect{i_sub}.StimCorr_t];
        Effect_allsub.StimCorr_pval     = ...
            [Effect_allsub.StimCorr_pval; Effect{i_sub}.StimCorr_pval];
        Effect_allsub.ToneProc_t        = ...
            [Effect_allsub.ToneProc_t; Effect{i_sub}.ToneProc_t];
        Effect_allsub.ToneProc_pval     = ...
            [Effect_allsub.ToneProc_pval; Effect{i_sub}.ToneProc_pval];
    
        Effect_allsub.Pred_clusterstat              = ...
            [Effect_allsub.Pred_clusterstat; ...
            Effect{i_sub}.Pred_clusterstat];
        Effect_allsub.Pred_cluster_pval             = ...
            [Effect_allsub.Pred_cluster_pval; ...
            Effect{i_sub}.Pred_cluster_pval];
        Effect_allsub.SimplePredErr_clusterstat     = ...
            [Effect_allsub.SimplePredErr_clusterstat; ...
            Effect{i_sub}.SimplePredErr_clusterstat];
        Effect_allsub.SimplePredErr_cluster_pval    = ...
            [Effect_allsub.SimplePredErr_cluster_pval; ...
            Effect{i_sub}.SimplePredErr_cluster_pval];
        Effect_allsub.ComplexPredErr_clusterstat    = ...
            [Effect_allsub.ComplexPredErr_clusterstat; ...
            Effect{i_sub}.ComplexPredErr_clusterstat];
        Effect_allsub.ComplexPredErr_cluster_pval   = ...
            [Effect_allsub.ComplexPredErr_cluster_pval; ...
            Effect{i_sub}.ComplexPredErr_cluster_pval];  
        
        Effect_allsub.SHI_Kprime    = ...
            [Effect_allsub.SHI_Kprime; Effect{i_sub}.SHI_Kprime];
        Effect_allsub.SHI_pval      = ...
            [Effect_allsub.SHI_pval; Effect{i_sub}.SHI_pval];  
    end
    
    %Cleanup
    usedElecs_chanposIndex = [];
    clear DataClean_AllTrials
    disp(['-- done loading & setting up data struct in ' ...
        num2str(toc) ' sec --'])    
end

%% 2) Determine sign. electrode index and labels for each effect
clear SignElecs
SignElecs       = struct;
curr_fieldnames = fieldnames(Effect_allsub);

for i_subfield = 1:length(curr_fieldnames)    
    curr_fieldname  = curr_fieldnames{i_subfield};
    if ~contains(curr_fieldname, '_pval')
        continue
    end    
    SignElecs.(curr_fieldname) = struct;
    
    if size(Effect_allsub.(curr_fieldname),2) > 1
        for i_TW = 1:size(Effect_allsub.(curr_fieldname),2)
            SignElecs.(curr_fieldname).array(:,i_TW)    = ...
                Effect_allsub.(curr_fieldname)(:,i_TW) < param.pval_plotting;
            SignElecs.(curr_fieldname).index{i_TW}    = ...
                find(SignElecs.(curr_fieldname).array(:,i_TW));
            SignElecs.(curr_fieldname).numelecs(:,i_TW) = ...
                length(SignElecs.(curr_fieldname).index{i_TW});
            for i_sub = unique(Effect_allsub.sub_index)'
                SignElecs.(curr_fieldname).numelecs_persub(i_TW,i_sub) = ...
                    sum(SignElecs.(curr_fieldname).array(...
                    find(Effect_allsub.sub_index == i_sub),i_TW));
            end

            SignElecs.(curr_fieldname).labels{i_TW}       = ...
                ElecLabels_allsub.elec_labels(SignElecs.(curr_fieldname).index{i_TW});
            SignElecs.(curr_fieldname).sublabels{i_TW}    = ...
                ElecLabels_allsub.elec_sublabels(SignElecs.(curr_fieldname).index{i_TW});
            SignElecs.(curr_fieldname).anatlabels{i_TW}   = ...
                ElecLabels_allsub.elec_anatlabels(SignElecs.(curr_fieldname).index{i_TW});
            for i_elec = 1:length(SignElecs.(curr_fieldname).labels{i_TW})
                SignElecs.(curr_fieldname).fulllabels{i_TW}{i_elec,1} = ...
                [SignElecs.(curr_fieldname).sublabels{i_TW}{i_elec}, ' ', ...
                SignElecs.(curr_fieldname).anatlabels{i_TW}{i_elec}];            
            end             
        end
    else
        SignElecs.(curr_fieldname).array    = ...
            Effect_allsub.(curr_fieldname) < param.pval_plotting;
        SignElecs.(curr_fieldname).index    = ...
            find(SignElecs.(curr_fieldname).array);
        SignElecs.(curr_fieldname).numelecs = ...
            length(SignElecs.(curr_fieldname).index);
        SignElecs.(curr_fieldname).numelecs_persub = [];
        for i_sub = unique(Effect_allsub.sub_index)'
            find(Effect_allsub.sub_index == i_sub);
            SignElecs.(curr_fieldname).numelecs_persub(i_sub) = ...
                sum(SignElecs.(curr_fieldname).array(...
                find(Effect_allsub.sub_index == i_sub)));
        end

        SignElecs.(curr_fieldname).labels       = ...
            ElecLabels_allsub.elec_labels(SignElecs.(curr_fieldname).index);
        SignElecs.(curr_fieldname).sublabels    = ...
            ElecLabels_allsub.elec_sublabels(SignElecs.(curr_fieldname).index);
        SignElecs.(curr_fieldname).anatlabels   = ...
            ElecLabels_allsub.elec_anatlabels(SignElecs.(curr_fieldname).index);
        SignElecs.(curr_fieldname).fulllabels   = [];
        for i_elec = 1:length(SignElecs.(curr_fieldname).labels)
            SignElecs.(curr_fieldname).fulllabels{i_elec} = ...
            [SignElecs.(curr_fieldname).sublabels{i_elec}, ' ', ...
            SignElecs.(curr_fieldname).anatlabels{i_elec}];            
        end  
    end
end

%% 3. Plot comparison between sign. elecs for each subject
%Plot venn diagram & barchart with sign. elecs per effect rel. to all elecs & with
%overlapping elecs

%Plot surface plot differently marking effects

%% 4. Plot comparison between sign. elecs across all subjects
%Plot barchart with sign. elecs per effect rel. to all elecs & with
%overlapping elecs

%Plot surface plot differently marking effects
