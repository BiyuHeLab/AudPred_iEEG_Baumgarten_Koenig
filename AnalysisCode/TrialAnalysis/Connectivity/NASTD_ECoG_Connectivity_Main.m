%% Project: NASTD_ECoG
%Compute Granger Causality (GC) directed connectivity measures to
%determine:
%1) the interplay between frontal and temporal prediction effect electrodes (time-wise GC in the HP2-LP30Hz band)
%2) the interplay between frontal and temporal complex prediciton error (PE) effect eletrodes (frequency-wise GC in the high gamma band)
%3) the interplay between frontal prediction and temporal PE effect electrodes
%4) the interplay between global SHI and prediciton effect electrodes

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add project base path
NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
ToneDur_text = {'0.2' '0.4'};

plot_poststepFigs = 1;
save_poststepFigs = 1;

%% 1) Select and plot electrodes for which GC will be analyzed

% param.pval_plotting     = 0.01; %Pval thresh for plotting
% param.pval_FDR          = 0.05;
% param.FDRcorrect        = 0;
% param.plot_SubplotperTW = 1; %Common plot across TW
% param.ElecSelect        = 'All'; %StimCorr, All
% param.SamplesTW         = 25; %25 = 50ms TW
% param.Label_TW          = num2str(round(param.SamplesTW/512,2));
% param.ToneIndex         = 33;
% 
% InputDataType           = {'HP05toLP30Hz'};
% InputEffectType         = {'PredEffect', 'ComplexPredErrEffect'};
% ToneDur_text            = {'0.2' '0.4'};
% subs                    = sub_list(vars.validSubjs);
% 
% %1.1 Highlight same-subject electrodes with selected p-value for either prediciton OR complex prediction error effect
% for i_effect = 1:length(InputEffectType)
%     for i_inputData = 1:length(InputDataType)    
%         NASTD_ECoG_Connectivity_PlotSignElec_AllSubTD...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)
%     end
% end

% %1.2 Highlight same-subject electrodes with selected p-value 
% %with prediciton (for HP2toLP30Hz) AND complex prediction error (for HighGamma_LogAmp) effect
% NASTD_ECoG_Connectivity_PlotSignElecPred2PE_AllSubTD...
%     (subs, ...
%     ToneDur_text, ...
%     param,...
%     save_poststepFigs, paths_NASTD_ECoG)
 
% %1.3 Create output file specifiying selected electrodes
subs                    = sub_list(vars.validSubjs);
plot_poststepFigs       = 0;

param.pval_plotting     = 0.05; %Pval thresh 
param.FDRcorrect        = 0;
param.pval_FDR          = 0.05; 
param.ElecSelect        = 'All'; %StimCorr, All
InputDataType           = {'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

%Effect-based p-val thresholded electrode selection
SelElecs = NASTD_ECoG_Connectivity_ReadOutGCElecs_AllSubTD...
    (subs, ...
    InputDataType, ToneDur_text, ...
    param, plot_poststepFigs, ...
    paths_NASTD_ECoG);
 
% %%Select all existing electrodes (independent of effect or threshold)
% SelElecs = NASTD_ECoG_Connectivity_ReadOutAllElecs_AllSubTD...
%     (subs, ...
%     plot_poststepFigs, ...
%     paths_NASTD_ECoG);

%1.5 Create Plot showing electrode connections that are the basis for GC calculation
%Across subjects and TD, color-coded for region, sign-coded for prediction-effect-type.
%Plot separately for each electrode-selection based on different prediction-effect-types.

NASTD_ECoG_Connectivity_PlotSignElecConnectionsforGC_AllSubTD...
    (subs, ...
    InputDataType, ToneDur_text, ...
    param, plot_poststepFigs, ...
    paths_NASTD_ECoG);

NASTD_ECoG_Connectivity_PlotSignElecCon_AllSubTDLobes... %For Pred-Pred & PE-PE, all lobes together with optional connection lines
    (subs, ...
    InputDataType, ToneDur_text, ...
    param, plot_poststepFigs, ...
    paths_NASTD_ECoG);

%% 2) Compute Granger Causality (GC) for specific electrode combinations

%Hypotheses: 
    %Spatial: 
        %Frontal P -> Temporal P > Chance
        %Temporal PE -> Frontal PE > Chance (CAVE: few frontal PE)
        %Frontal P -> Temporal PE > Frontal P <- Temporal PE
    %Temporal: 
        %P -> PE electrodes stronger in earlier TW compared to late TW
        %P -> PE electrodes stronger in earlier TW compared to P <- PE electrodes
    %Spectral: 
        %P -> PE: low frequencies (alpha/beta)
        %P <- PE: high gamma
        %P -> P: low frequencies (alpha/beta)
        %PE -> PE: high gamma
        
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/');
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.01uncorr.mat') %p < 0.01 thresh,uncorr elec selection    
load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr.mat') %p < 0.05 thresh,uncorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.025FDRcorr.mat') %p < 0.025 thresh, FDRcorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_AllElecs.mat') % no thresh, all elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr_HighGamma.mat') %p < 0.05 thresh,uncorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.01uncorr_HighGamma.mat') %p < 0.05 thresh,uncorr elec selection    

param.GC            = [];
param.GC.fs         = 512;
param.GC.downsample = 0; %Cave: Downsampling lead to problems with tone sample selection
    if param.GC.downsample == 0
        param.GC.newfs  = param.GC.fs;
    else
        param.GC.newfs  = param.GC.fs/param.GC.downsample;
    end
param.GC.nvars      = 2;
param.GC.regmode    = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default
param.GC.maxmorder  = 50;   % maximum model order for model order estimation, rule of thumb = number of samples per input data snippet, but not > 100
% param.GC.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)

param.GC.tstat      = 'F'; % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
param.GC.alpha      = 0.05;   % significance level for significance test
param.GC.mhtc       = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

param.GC.ElecPairEffect   = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

% param.GC.tone_aggregation = {32, 33, 34};
% param.GC.tone_aggregation = {31, 33, 34};
% param.GC.tone_aggregation = {1, 6, 11, 16, 21, 26, 31};
% param.GC.tone_aggregation = {2, 7, 12, 17, 22, 27, 32};
% param.GC.tone_aggregation = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

param.GC.tone_aggregation_label = cellfun(@(c)[c],param.GC.tone_aggregation);
param.GC.tone_aggregation_label = num2str((param.GC.tone_aggregation_label),'_%d');
param.GC.tone_aggregation_label = erase(param.GC.tone_aggregation_label, ' ');

param.GC.epochdur_ms = 'full'; %'full', 200, 100
if strcmp(param.GC.epochdur_ms, 'full')
    param.GC.epochdur_ms_label = 'fullTW';
else
   param.GC.epochdur_ms_label = [num2str(param.GC.epochdur_ms) 'msTW'];
end

for i_sub = vars.validSubjs
        
        tic    
        GCdata{i_sub} = ... %Compute GC for aggregated tones
            NASTD_ECoG_Connectivity_CalculateGC_aggrtone ...
            (sub_list, i_sub, ...
            ToneDur_text, SelElecs, ...
            param, paths_NASTD_ECoG);
%         GCdata{i_sub} = ... %Compute GC for each single tone 
%             NASTD_ECoG_Connectivity_CalculateGC_pertone ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             param, paths_NASTD_ECoG);
        disp(['-- GC computation for sub: ' sub_list{i_sub} ' finished after ' num2str(round(toc/60),2) ' min --'])
   
end

%save output
path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
label_outputdata = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
    param.GC.ElecPairRegion '_' ...
    param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
    param.GC.epochdur_ms_label '_ensemnorm'];
save([path_outputdata label_outputdata], 'GCdata','-v7.3')
% load([path_outputdata label_outputdata], 'GCdata')

%Plot GC results
AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
    {'AntPFC'; ...
    'VentralT'; ...
    'SupParLob'};
%     AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
%         {'AntPFC_IFG'; 'AntPFC_PrecentralG'; ...
%         'VentralT_STG'; 'VentralT_MTG'; ...
%         'SupParLob_PostcentralG'; 'SupParLob_SupramarginalG'};
% AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
%     {'AntPFC'; 'AntPFC_IFG'; 'AntPFC_PrecentralG'; ...
%     'VentralT'; 'VentralT_STG'; 'VentralT_MTG'; ...
%     'SupParLob'; 'SupParLob_PostcentralG'; 'SupParLob_SupramarginalG'; ...
%     'OccipitalL'};

    %Single subject per electrode pairing
    %     NASTD_ECoG_Connectivity_PlotGC_persub_aggrtone ...
    %         (sub_list(vars.validSubjs), ...
    %         ToneDur_text, ...
    %         save_poststepFigs, ...
    %         param, paths_NASTD_ECoG);
    % NASTD_ECoG_Connectivity_PlotGC_persub_pertone ...
    %     (sub_list(vars.validSubjs), ...
    %     ToneDur_text, InputDataType, label_ElecPairSel, ...
    %     save_poststepFigs, ...
    %     paths_NASTD_ECoG);
       
    %Single subject aggregated GC results
%     NASTD_ECoG_Connectivity_PlotGC_persub_aggrtoneanat ...
%         (sub_list(vars.validSubjs), ...
%         ToneDur_text, SelElecs, AnatReg_CatLabels, ...
%         save_poststepFigs, ...
%         param, paths_NASTD_ECoG);

    %Plot group-level aggregated GC results
%     NASTD_ECoG_Connectivity_PlotGC_allsub_aggrtone ...
%         (sub_list(vars.validSubjs), ...
%         ToneDur_text, SelElecs, AnatReg_CatLabels, ...
%         save_poststepFigs, ...
%         param, paths_NASTD_ECoG);
    NASTD_ECoG_Connectivity_PlotGC_allsub_aggrtone_median ...
        (sub_list(vars.validSubjs), ...
        ToneDur_text, SelElecs, AnatReg_CatLabels, ...
        save_poststepFigs, ...
        param, paths_NASTD_ECoG);
    
    %Plot GC output with statistics
    if strcmp(param.GC.tone_aggregation_label,'_32_33_34')
        NASTD_ECoG_Connectivity_PlotGC_statp32p33p34 ...
            (sub_list(vars.validSubjs), ...
            ToneDur_text, SelElecs, AnatReg_CatLabels, ...
            save_poststepFigs, ...
            param, paths_NASTD_ECoG);
    else
        NASTD_ECoG_Connectivity_PlotGC_statp1top31 ...
            (sub_list(vars.validSubjs), ...
            ToneDur_text, SelElecs, AnatReg_CatLabels, ...
            save_poststepFigs, ...
            param, paths_NASTD_ECoG);
    end
   
    


%% 3) Compute Granger Causality (GC) during p34 for different PE effects 
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/');
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.01uncorr.mat') %p < 0.01 thresh,uncorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr.mat') %p < 0.05 thresh,uncorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.025FDRcorr.mat') %p < 0.025 thresh, FDRcorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_AllElecs.mat') % no thresh, all elec selection    
load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr_HighGamma.mat') %p < 0.05 thresh,uncorr elec selection    
% load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.01uncorr_HighGamma.mat') %p < 0.05 thresh,uncorr elec selection    

param.GC.fs         = 512;
param.GC.downsample = 0; %Cave: Downsampling lead to problems with tone sample selection
    if param.GC.downsample == 0
        param.GC.newfs  = param.GC.fs;
    else
        param.GC.newfs  = param.GC.fs/param.GC.downsample;
    end
param.GC.nvars      = 2;
param.GC.regmode    = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default
param.GC.maxmorder  = 50;   % maximum model order for model order estimation, rule of thumb = number of samples per input data snippet, but not > 100
param.GC.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)

param.GC.tstat      = 'F'; % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
param.GC.alpha      = 0.05;   % significance level for significance test
param.GC.mhtc       = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

param.GC.ElecPairEffect   = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

param.GC.epochdur_ms = 'full'; %'full', 200, 100
if strcmp(param.GC.epochdur_ms, 'full')
    param.GC.epochdur_ms_label = 'fullTW';
else
   param.GC.epochdur_ms_label = [num2str(param.GC.epochdur_ms) 'msTW'];
end

for i_sub = vars.validSubjs
        
    tic
    GCdata{i_sub} = ...
        NASTD_ECoG_Connectivity_CalculateGC_PEeffect ...
        (sub_list, i_sub, ...
        ToneDur_text, SelElecs, ...
        param, paths_NASTD_ECoG);
    
    disp(['-- GC PE computation for sub: ' sub_list{i_sub} ' finished after ' num2str(round(toc/60),2) ' min --'])
        
end

%save output
path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/PEeffect/']);
if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
label_outputdata = ['GCdataPEeffect_n' num2str(length(vars.validSubjs)) '_' ...
    param.GC.ElecPairRegion '_' ...
    param.GC.InputDataType{1} '_' ...
    param.GC.epochdur_ms_label '_ensemnorm'];
save([path_outputdata label_outputdata], 'GCdata')
% load([path_outputdata label_outputdata], 'GCdata')

%Plot GC results
AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
    {'AntPFC'; ...
    'VentralT'; ...
    'SupParLob'};

NASTD_ECoG_Connectivity_PlotGC_statPEeffect ...
    (sub_list(vars.validSubjs), ...
    ToneDur_text, SelElecs, AnatReg_CatLabels, ...
    save_poststepFigs, ...
    param, paths_NASTD_ECoG);
