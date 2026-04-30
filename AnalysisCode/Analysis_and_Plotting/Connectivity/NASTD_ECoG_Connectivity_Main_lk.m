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

plot_poststepFigs = 0;
save_poststepFigs = 1;

%% 1) Select and plot electrodes for which GC will be analyzed

param.pval_plotting     = 0.01; %Pval thresh for plotting
param.pval_FDR          = 0.05;
param.FDRcorrect        = 0;
param.plot_SubplotperTW = 1; %Common plot across TW
param.ElecSelect        = 'All'; %StimCorr, All
param.SamplesTW         = 25; %25 = 50ms TW
param.Label_TW          = num2str(round(param.SamplesTW/512,2));
param.ToneIndex         = 33;

InputDataType           = {'HP05toLP30Hz'};
InputEffectType         = {'PredEffect', 'ComplexPredErrEffect'};
ToneDur_text            = {'0.2' '0.4'};
subs                    = sub_list(vars.validSubjs);

%1.1 Highlight same-subject electrodes with selected p-value for either prediciton OR complex prediction error effect
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
InputDataType           = {'HP05toLP30Hz'};%, 'HighGamma_LogAmp'}; 

%Effect-based p-val thresholded electrode selection
% SelElecs = NASTD_ECoG_Connectivity_ReadOutGCElecs_AllSubTD...
%     (subs, ...
%     InputDataType, ToneDur_text, ...
%     param, plot_poststepFigs, ...
%     paths_NASTD_ECoG);
 
% %%Select all existing electrodes (independent of effect or threshold)
% SelElecs = NASTD_ECoG_Connectivity_ReadOutAllElecs_AllSubTD...
%     (subs, ...
%     plot_poststepFigs, ...
%     paths_NASTD_ECoG);

%1.5 Create Plot showing electrode connections that are the basis for GC calculation
%Across subjects and TD, color-coded for region, sign-coded for prediction-effect-type.
%Plot separately for each electrode-selection based on different prediction-effect-types.

% NASTD_ECoG_Connectivity_PlotSignElecConnectionsforGC_AllSubTD...
%     (subs, ...
%     InputDataType, ToneDur_text, ...
%     param, plot_poststepFigs, ...
%     paths_NASTD_ECoG);
% 
% NASTD_ECoG_Connectivity_PlotSignElecCon_AllSubTDLobes... %For Pred-Pred & PE-PE, all lobes together with optional connection lines
%     (subs, ...
%     InputDataType, ToneDur_text, ...
%     param, plot_poststepFigs, ...
%     paths_NASTD_ECoG);

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

% Computing how many significant electrodes per lobe (frontal, parietal, and temporal) and prediction vs. PE effect for Table S2
% Define keywords for each lobe
% frontalKeywords = ["AntPFC", "PrecentralG", "IFG"];
% parietalKeywords = ["SupParLob", "PostcentralG", "SupramarginalG"];
% temporalKeywords = ["VentralT", "STG", "MTG"];
% 
% % Extract AnatCat Label column
% predLabels = SelElecs.PredEffect.("AnatCat Label");
% predLabelsStr = string(predLabels);
% predLabelsStr(contains(predLabelsStr, "OccipitalL")) = []; % Remove entries containing "OccipitalL"
% [uniqueLabels1, ~, idx1] = unique(predLabelsStr);
% counts1 = accumarray(idx1, 1);
% PredlabelCountsTable = table(uniqueLabels1, counts1, 'VariableNames', {'Label', 'Count'});
% 
% peLabels = SelElecs.PEEffect.("AnatCat Label");
% peLabelsStr = string(peLabels);
% peLabelsStr(contains(peLabelsStr, "OccipitalL")) = []; % Remove entries containing "OccipitalL"
% [uniqueLabels2, ~, idx2] = unique(peLabelsStr);
% counts2 = accumarray(idx2, 1);
% PElabelCountsTable = table(uniqueLabels2, counts2, 'VariableNames', {'Label', 'Count'});
% 
% countLobes = @(data, keywords) sum(any(contains(data, keywords), 2));
% 
% % Count electrodes per lobe for each effect using the strict word-boundary matching
% frontalPredCount = countLobes(predLabelsStr, frontalKeywords);
% parietalPredCount = countLobes(predLabelsStr, parietalKeywords);
% temporalPredCount = countLobes(predLabelsStr, temporalKeywords);
% 
% frontalPECount = countLobes(peLabelsStr, frontalKeywords);
% parietalPECount = countLobes(peLabelsStr, parietalKeywords);
% temporalPECount = countLobes(peLabelsStr, temporalKeywords);
% 
% % Create the output table
% LobeSummaryTable = table(["Frontal"; "Parietal"; "Temporal"], ...
%                          [frontalPredCount; parietalPredCount; temporalPredCount], ...
%                          [frontalPECount; parietalPECount; temporalPECount], ...
%                          'VariableNames', {'Lobe', 'PredCount', 'PECount'});
% writetable(LobeSummaryTable, ['/isilon/LFMI/VMdrive/Lua/NASTD/Figures/ElectrodeCountsPerLobe.csv']);

%% Compute GC
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
%param.GC.ElecPairEffect = {'Pred_Pred', 'PE_PE'};
%param.GC.ElecPairEffect = {'Pred_PE', 'PE_Pred'};
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

% param.GC.tone_aggregation = {32, 33, 34};
% param.GC.tone_aggregation = {31, 33, 34};
%param.GC.tone_aggregation = {1, 6, 11, 16, 21, 26, 31};
%param.GC.tone_aggregation = {2, 7, 12, 17, 22, 27, 32};
% param.GC.tone_aggregation = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
% param.GC.tone_aggregation = {30, 31, 32};
%param.GC.tone_aggregation = {30, 31, 32, 33, 34};
param.GC.tone_aggregation = {34};

param.GC.tone_aggregation_label = cellfun(@(c)[c],param.GC.tone_aggregation);
param.GC.tone_aggregation_label = num2str((param.GC.tone_aggregation_label),'_%d');
param.GC.tone_aggregation_label = erase(param.GC.tone_aggregation_label, ' ');

param.GC.epochdur_ms = 'full'; %'full', 200, 100
if strcmp(param.GC.epochdur_ms, 'full')
    param.GC.epochdur_ms_label = 'fullTW';
else
   param.GC.epochdur_ms_label = [num2str(param.GC.epochdur_ms) 'msTW'];
end

% Across-trial GC estimates
% parfor i_sub = vars.validSubjs
%         
%         tic    
%         GCdata{i_sub} = ... %Compute GC for aggregated tones
%             NASTD_ECoG_Connectivity_CalculateGC_aggrtone ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             param, paths_NASTD_ECoG);
%         GCdata{i_sub} = ... %Compute GC for each single tone 
%             NASTD_ECoG_Connectivity_CalculateGC_pertone ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             param, paths_NASTD_ECoG);
%         disp(['-- GC computation for sub: ' sub_list{i_sub} ' finished after ' num2str(round(toc/60),2) ' min --'])
%    
% end

%save output
% path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
% if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
% label_outputdata = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
%     param.GC.ElecPairRegion '_' ...
%     param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
%     param.GC.epochdur_ms_label '_ensemnorm'];
% %save([path_outputdata label_outputdata], 'GCdata','-v7.3')
% load([path_outputdata label_outputdata], 'GCdata')

% Across-trial GC estimates for selected trial IDs
FTPL_ratings = load([paths_NASTD_ECoG.Analysis_Behavior 'FTPLratings/Allsub_n8/Data/', 'Trialwise_FTPL_allsubs.mat']);
FTPL_ratings = FTPL_ratings.Trialwise_FTPL;
rows_all = strcmp(FTPL_ratings.ToneDur, 'all'); %Take ratings across both TD conditions
FTPL_all = FTPL_ratings(rows_all, :);
sub_list_FTPL = sub_list(2:9);

% Averaging over trialwise data for low-FTPL trials and for high-FTPL
% trials
low_FTPL_all = {};
high_FTPL_all = {};
for i_sub = 1:length(sub_list_FTPL)
    current_subID = sub_list_FTPL(i_sub);
    
    % Find all rows in FTPL_all for this subject
    subj_rows = strcmp(FTPL_all.SubID, current_subID);
    
    % Get the FTPL ratings for this subject
    ratings_subj = FTPL_all.FTPLrating(subj_rows);
    median_subj = median(ratings_subj, 'omitnan');
    
    % Create logical indices for low and high FTPL relative to
    % subject-specific median of ratings
    low_idx = ratings_subj < median_subj;
    low_idx = FTPL_all.TrialNum(low_idx); % TrialIDs for low FTPL
    low_FTPL_all{i_sub} = low_idx;
    
    high_idx = ratings_subj > median_subj;
    high_idx = FTPL_all.TrialNum(high_idx);
    high_FTPL_all{i_sub} = high_idx;
end

path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
% parfor i_sub = vars.validSubjs(2:9)
%         tic    
%         GCdata_lowFTPL{i_sub} = ... %Compute GC for aggregated tones
%             NASTD_ECoG_Connectivity_CalculateGC_aggrtone_SelTrials ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             low_FTPL_all{i_sub-1}, ...
%             param, paths_NASTD_ECoG);
%         disp(['-- GC computation for sub: ' sub_list{i_sub} ' for LOW FTPL finished after ' num2str(round(toc/60),2) ' min --'])
%         GCdata_highFTPL{i_sub} = ... %Compute GC for aggregated tones
%             NASTD_ECoG_Connectivity_CalculateGC_aggrtone_SelTrials ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             high_FTPL_all{i_sub-1}, ...
%             param, paths_NASTD_ECoG);
%         disp(['-- GC computation for sub: ' sub_list{i_sub} ' for LOW FTPL finished after ' num2str(round(toc/60),2) ' min --'])
% end
label_outputdata1 = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
            param.GC.ElecPairRegion '_' ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
            param.GC.epochdur_ms_label '_lowFTPL'];
%save([path_outputdata label_outputdata1], 'GCdata_lowFTPL','-v7.3')
GCdata_lowFTPL=load([path_outputdata label_outputdata1]);

label_outputdata2 = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
            param.GC.ElecPairRegion '_' ...
            param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
            param.GC.epochdur_ms_label '_highFTPL'];
%save([path_outputdata label_outputdata2], 'GCdata_highFTPL','-v7.3')
GCdata_highFTPL = load([path_outputdata label_outputdata2]);

% Trialwise GC estimates
% parfor i_sub = vars.validSubjs
%         
%         tic    
%         GCdata{i_sub} = ... %Compute GC for each single tone and for each trial; used for control analysis checking how GC values differ for different trial-wise FTPL ratings; only run for tone 34
%             NASTD_ECoG_Connectivity_CalculateGC_aggrtone_trialwise ...
%             (sub_list, i_sub, ...
%             ToneDur_text, SelElecs, ...
%             param, paths_NASTD_ECoG);
%         disp(['-- GC computation for sub: ' sub_list{i_sub} ' finished after ' num2str(round(toc/60),2) ' min --'])
%    
% end
%save output
% path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
% if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
% label_outputdata = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
%     param.GC.ElecPairRegion '_' ...
%     param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
%     param.GC.epochdur_ms_label '_trialwise_ensemnorm'];
% %save([path_outputdata label_outputdata], 'GCdata','-v7.3')
% load([path_outputdata label_outputdata], 'GCdata')

%% Compare GC values for low vs. high FTPL ratings for parietal -> frontal during p34 using trialwise data

%First make a GCdata for low FTPL ratings, average across trials,NASTD_ECoG_Connectivity_stat save GC.
%Do the same for high FTP ratings, average across trials, save GC.
%Pass both averaged GC through anatomical averaging script (NASTD_ECoG_Connectivity_PlotGC_statp32p33p34.m in Thomas' folder)
% This will output 10 (source anat regions) x 10 (target anat regions) x 2
% (TDs) arrays for the low-FTPL-GC and the high-FTPL-GC. Then we can select
% parietal to frontal by indexing source anat region = 7 and target anat
% region = 1 for each of them, and do the stats.
% GCdata_FTPL = GCdata(2:9);
% numSubjects = length(GCdata_FTPL);
% FTPL_ratings = load([paths_NASTD_ECoG.Analysis_Behavior 'FTPLratings/Allsub_n8/Data/', 'Trialwise_FTPL_allsubs.mat']);
% FTPL_ratings = FTPL_ratings.Trialwise_FTPL;
% rows_all = strcmp(FTPL_ratings.ToneDur, 'all'); %Take ratings across both TD conditions
% FTPL_all = FTPL_ratings(rows_all, :);
% all_trialIDs = cell(numSubjects, 1);
% for i_sub = 1:numSubjects
%     % Combine trial IDs across tone durations (1x2 cell)
%     all_trialIDs{i_sub} = [GCdata{i_sub}.info.trialIDs{1}; GCdata{i_sub}.info.trialIDs{2}];
% end
% sub_list_FTPL = sub_list(2:9);
% 
% % Averaging over trialwise data for low-FTPL trials and for high-FTPL
% % trials
% for i_sub = 1:numSubjects
%     current_subID = sub_list_FTPL(i_sub);
%     
%     % Find all rows in FTPL_all for this subject
%     subj_rows = strcmp(FTPL_all.SubID, current_subID);
%     
%     % Get the FTPL ratings for this subject
%     ratings_subj = FTPL_all.FTPLrating(subj_rows);
%     median_subj = median(ratings_subj, 'omitnan');
%     trialIDs_subj = all_trialIDs{i_sub};
%     
%     % Create logical indices for low and high FTPL relative to
%     % subject-specific median of ratings
%     low_idx = ratings_subj < median_subj;
%     low_idx = FTPL_all.TrialNum(low_idx); % TrialIDs for low FTPL
%     
%     high_idx = ratings_subj > median_subj;
%     high_idx = FTPL_all.TrialNum(high_idx);
%     
%     % Find the Trial IDs for the GC data
%     trialIDs_GC_TD1 = GCdata_FTPL{i_sub}.info.trialIDs{1}; %Trial IDs of GC data for TD1
%     trialIDs_GC_TD2 = GCdata_FTPL{i_sub}.info.trialIDs{2}; %Trial IDs for GC data for TD2
%     
%     % Identify corresponding trials in GC data for low and high FTPL
%     % ratings
%     [~,match_idx_low_TD1] = ismember(low_idx, trialIDs_GC_TD1); % Finds corresponding trials in the GC data for TD1
%     idx_low_in_GC_TD1 = match_idx_low_TD1(match_idx_low_TD1 > 0);
%     [~, match_idx_low_TD2] = ismember(low_idx, trialIDs_GC_TD2);  % Finds corresponding trials in the GC data for TD2
%     idx_low_in_GC_TD2 = match_idx_low_TD2(match_idx_low_TD2 > 0);
%     
%     [~, match_idx_high_TD1] = ismember(high_idx, trialIDs_GC_TD1); % Finds corresponding trials in the GC data for TD1
%     idx_high_in_GC_TD1 = match_idx_high_TD1(match_idx_high_TD1 > 0);
%     [~, match_idx_high_TD2] = ismember(high_idx, trialIDs_GC_TD2); % Finds corresponding trials in the GC data for TD2
%     idx_high_in_GC_TD2 = match_idx_high_TD2(match_idx_high_TD2 > 0);
%     
%     for i_effect = 1:length(GCdata_FTPL{i_sub}.temporalGC)
%         for i_TD = 1:length(GCdata_FTPL{i_sub}.temporalGC{i_effect})
%             if i_TD == 1
%                 idx_low_in_GC = idx_low_in_GC_TD1;
%                 idx_high_in_GC = idx_high_in_GC_TD1;
%             else
%                 idx_low_in_GC = idx_low_in_GC_TD2;
%                 idx_high_in_GC = idx_high_in_GC_TD2;
%             end
%             
%            % === TEMPORAL GC ===
%             temp_struct = GCdata_FTPL{i_sub}.temporalGC{i_effect}{i_TD};
%             GCdata_average_lowFTPL{i_sub}.temporalGC{i_effect}{i_TD}.source2target = ...
%                 mean(temp_struct.source2target(:,:,:,idx_low_in_GC), 4, 'omitnan');
%             GCdata_average_lowFTPL{i_sub}.temporalGC{i_effect}{i_TD}.target2source = ...
%                 mean(temp_struct.target2source(:,:,:,idx_low_in_GC), 4, 'omitnan');
% 
%             GCdata_average_highFTPL{i_sub}.temporalGC{i_effect}{i_TD}.source2target = ...
%                 mean(temp_struct.source2target(:,:,:,idx_high_in_GC), 4, 'omitnan');
%             GCdata_average_highFTPL{i_sub}.temporalGC{i_effect}{i_TD}.target2source = ...
%                 mean(temp_struct.target2source(:,:,:,idx_high_in_GC), 4, 'omitnan');
% 
%             % === PVAL TEMPORAL GC ===
%             pval_struct = GCdata_FTPL{i_sub}.pval_temporalGC{i_effect}{i_TD};
%             GCdata_average_lowFTPL{i_sub}.pval_temporalGC{i_effect}{i_TD}.source2target = ...
%                 mean(pval_struct.source2target(:,:,:,idx_low_in_GC), 4, 'omitnan');
%             GCdata_average_lowFTPL{i_sub}.pval_temporalGC{i_effect}{i_TD}.target2source = ...
%                 mean(pval_struct.target2source(:,:,:,idx_low_in_GC), 4, 'omitnan');
% 
%             GCdata_average_highFTPL{i_sub}.pval_temporalGC{i_effect}{i_TD}.source2target = ...
%                 mean(pval_struct.source2target(:,:,:,idx_high_in_GC), 4, 'omitnan');
%             GCdata_average_highFTPL{i_sub}.pval_temporalGC{i_effect}{i_TD}.target2source = ...
%                 mean(pval_struct.target2source(:,:,:,idx_high_in_GC), 4, 'omitnan');
% 
%             % === SPECTRAL GC ===
%             spec_struct = GCdata_FTPL{i_sub}.spectralGC{i_effect}{i_TD};
%             GCdata_average_lowFTPL{i_sub}.spectralGC{i_effect}{i_TD}.source2target = ...
%                 mean(spec_struct.source2target(:,:,:,:,idx_low_in_GC), 5, 'omitnan');
%             GCdata_average_lowFTPL{i_sub}.spectralGC{i_effect}{i_TD}.target2source = ...
%                 mean(spec_struct.target2source(:,:,:,:,idx_low_in_GC), 5, 'omitnan');
% 
%             GCdata_average_highFTPL{i_sub}.spectralGC{i_effect}{i_TD}.source2target = ...
%                 mean(spec_struct.source2target(:,:,:,:,idx_high_in_GC), 5, 'omitnan');
%             GCdata_average_highFTPL{i_sub}.spectralGC{i_effect}{i_TD}.target2source = ...
%                 mean(spec_struct.target2source(:,:,:,:,idx_high_in_GC), 5, 'omitnan');
%         end
%     end
%     GCdata_average_lowFTPL{i_sub}.info = GCdata_FTPL{i_sub}.info;
%     GCdata_average_lowFTPL{i_sub}.label_allelecs = GCdata_FTPL{i_sub}.label_allelecs;
%     GCdata_average_lowFTPL{i_sub}.ind_pairedelecs_from = GCdata_FTPL{i_sub}.ind_pairedelecs_from;
%     GCdata_average_lowFTPL{i_sub}.ind_pairedelecs_to = GCdata_FTPL{i_sub}.ind_pairedelecs_to;
%     
%     GCdata_average_highFTPL{i_sub}.info = GCdata_FTPL{i_sub}.info;
%     GCdata_average_highFTPL{i_sub}.label_allelecs = GCdata_FTPL{i_sub}.label_allelecs;
%     GCdata_average_highFTPL{i_sub}.ind_pairedelecs_from = GCdata_FTPL{i_sub}.ind_pairedelecs_from;
%     GCdata_average_highFTPL{i_sub}.ind_pairedelecs_to = GCdata_FTPL{i_sub}.ind_pairedelecs_to;
% end

% AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
%     {'AntPFC'; ...
%     'VentralT'; ...
%     'SupParLob'};
% 
% Gavg_tempGC_peranatreg_lowFTPL = NASTD_ECoG_Connectivity_ObtainGCByAnatRegions ...
%     (sub_list_FTPL, ...
%     GCdata_average_lowFTPL, ...
%     ToneDur_text, ...
%     SelElecs, AnatReg_CatLabels, ...
%     save_poststepFigs, ...
%     param, paths_NASTD_ECoG); %This will have the format, for each field, of i_effect cells, with format nSourceElecs x nTargetElecs x nTD
% 
% Gavg_tempGC_peranatreg_highFTPL = NASTD_ECoG_Connectivity_ObtainGCByAnatRegions ...
%     (sub_list_FTPL, ...
%     GCdata_average_highFTPL, ...
%     ToneDur_text, ...
%     SelElecs, AnatReg_CatLabels, ...
%     save_poststepFigs, ...
%     param, paths_NASTD_ECoG);  %This will have the format, for each field, of i_effect cells, with format nSourceElecs x nTargetElecs x nTD


% Compute GC low vs high FTPL for Parietal > Frontal
% 
% index_of_interest_source = 7; %Parietal
% index_of_interest_target = 1; %Frontal
% 
% GC_Par_to_Front_lowFTPL = nan(numSubjects, 1);  % one entry per subject
% GC_Par_to_Front_highFTPL = nan(numSubjects, 1);  % one entry per subject
% 
% GC_Par_to_Front_lowFTPL_all = {}; % Pool across 4 different effect types, average across TDs
% GC_Par_to_Front_highFTPL_all = {};
% 
% for i_sub = 1:numSubjects
%     all_GC_rows_low = [];
%     all_GC_rows_high = [];
% 
%     for i_effect = 1:length(GCdata_FTPL{i_sub}.temporalGC)
%         for i_TD = 1:length(GCdata_FTPL{i_sub}.temporalGC{i_effect})
%             % Get data and subject indices for this effect and TD
%             GCvals_low = Gavg_tempGC_peranatreg_lowFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
%             subjIDs_low = Gavg_tempGC_peranatreg_lowFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
% 
%             if ~isempty(GCvals_low) && ~isempty(subjIDs_low)
%                 % Find rows belonging to this subject
%                 rows_low = find(subjIDs_low == i_sub);
%                 if ~isempty(rows_low)
%                     all_GC_rows_low = [all_GC_rows_low; GCvals_low(rows_low, :)];
%                 end
%             end
%             
%             % Get data and subject indices for this effect and TD
%             GCvals_high = Gavg_tempGC_peranatreg_highFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
%             subjIDs_high = Gavg_tempGC_peranatreg_highFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
% 
%             if ~isempty(GCvals_high) && ~isempty(subjIDs_high)
%                 % Find rows belonging to this subject
%                 rows_high = find(subjIDs_high == i_sub);
%                 if ~isempty(rows_high)
%                     all_GC_rows_high = [all_GC_rows_high; GCvals_high(rows_high, :)];
%                 end
%             end
%         end
%     end
% 
%     if ~isempty(all_GC_rows_low)
%         GC_Par_to_Front_lowFTPL(i_sub,1) = mean(all_GC_rows_low, 1, 'omitnan');
%     else
%         GC_Par_to_Front_lowFTPL(i_sub,1) = NaN;  % or [] if you prefer
%     end
%     
%     if ~isempty(all_GC_rows_high)
%         GC_Par_to_Front_highFTPL(i_sub,1) = mean(all_GC_rows_high, 1, 'omitnan');
%     else
%         GC_Par_to_Front_highFTPL(i_sub,1) = NaN;  % or [] if you prefer
%     end
% end
% 
% % Do stats
% [~, pval, ~, stats] = ttest(GC_Par_to_Front_lowFTPL, GC_Par_to_Front_highFTPL);
% [pval, h, stats] = signrank(GC_Par_to_Front_lowFTPL, GC_Par_to_Front_highFTPL);
% 
% % Plot
% % Combine data
% % all_data = [GC_Par_to_Front_lowFTPL(:); GC_Par_to_Front_highFTPL(:)];
% % group = [ones(length(GC_Par_to_Front_lowFTPL),1); 2*ones(length(GC_Par_to_Front_highFTPL),1)];
% all_data = [GC_Par_to_Front_lowFTPL(:); mean(GC_Par_to_Front_lowFTPL); GC_Par_to_Front_highFTPL(:); mean(GC_Par_to_Front_highFTPL)];
% group = [ones(length(GC_Par_to_Front_lowFTPL)+1,1); 2*ones(length(GC_Par_to_Front_highFTPL)+1,1)];
% 
% % Create boxplot
% figure; hold on;
% h = boxplot(all_data, group, 'Labels', {'Low FTPL', 'High FTPL'}, ...
%     'Colors', 'k', 'Widths', 0.5);
% 
% % Change box colors
% boxColors = [0.2 0.6 1; 1 0.4 0.4];  % Blue, Red
% boxes = findobj(gca, 'Tag', 'Box');
% for j = 1:length(boxes)
%     patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), ...
%         boxColors(3-j,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
% end
% 
% % Overlay paired data
% for i = 1:length(GC_Par_to_Front_lowFTPL)
%     plot([1, 2], [GC_Par_to_Front_lowFTPL(i), GC_Par_to_Front_highFTPL(i)], '-o', ...
%         'Color', [0.6 0.6 0.6], ...
%         'MarkerFaceColor', 'k', ...
%         'MarkerEdgeColor', 'none', ...
%         'LineWidth', 1.2);
% end
% 
% plot([1, 2], [mean(GC_Par_to_Front_lowFTPL), mean(GC_Par_to_Front_highFTPL)], '-o', ...
%         'Color', [0.6 0.6 0.6], ...
%         'MarkerFaceColor', 'k', ...
%         'MarkerEdgeColor', 'none', ...
%         'LineWidth', 1.2);
% 
% % Significance annotation
% [~, pval] = ttest(GC_Par_to_Front_lowFTPL, GC_Par_to_Front_highFTPL);
% y_max = max([GC_Par_to_Front_lowFTPL(:); GC_Par_to_Front_highFTPL(:)]) + 0.01;
% plot([1, 2], [y_max, y_max], 'k', 'LineWidth', 1.5)
% 
% if pval < 0.001
%     sig_label = '***';
% elseif pval < 0.01
%     sig_label = '**';
% elseif pval < 0.05
%     sig_label = '*';
% else
%     sig_label = 'n.s.';
% end
% text(1.5, y_max + 0.005, sig_label, ...
%     'FontSize', 16, ...
%     'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center');
% 
% % Formatting
% title('GC during p34 Parietal > Frontal');
% ylabel('Average GC');
% ylim([min([GC_Par_to_Front_lowFTPL; GC_Par_to_Front_highFTPL]) - 0.01, y_max + 0.02]);
% set(gca, 'FontSize', 12);
% box on;
% 
% filename = ['boxplot_compare_GC_lowhigh_FTPL_Par_to_Front.png'];
% figfile = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figures/' filename];
% saveas(gcf, figfile, 'png');
% 
% 
% %% Computing stats for FTPL low vs high GC comparison across pairs from different effect-type pairings using TRIALWISE data
% % Compute GC low vs high FTPL for Parietal > Frontal
% 
% index_of_interest_source = 7; %Parietal
% index_of_interest_target = 1; %Frontal
% effectTypeNames = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};
% 
% GC_Par_to_Front_HighMinusLow = struct();  % Final output
% 
% for i_effect = 1:length(GCdata_FTPL{1}.temporalGC)
%     effect_field = effectTypeNames{i_effect};
%     for i_TD = 1:length(GCdata_FTPL{1}.temporalGC{1})
%         if ~isfield(GC_Par_to_Front_HighMinusLow, effect_field)
%             GC_Par_to_Front_HighMinusLow.(effect_field) = struct();
%         end
%         if ~isfield(GC_Par_to_Front_HighMinusLow.(effect_field), sprintf('TD_%d', i_TD))
%             GC_Par_to_Front_HighMinusLow.(effect_field).(sprintf('TD_%d', i_TD)) = [];
%         end
% 
%         pooled_diffs = [];  % temporary holder for this TD across all subjects
% 
%         for i_sub = 1:numSubjects
%             % LOW
%             GCvals_low = Gavg_tempGC_peranatreg_lowFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
%             subjIDs_low = Gavg_tempGC_peranatreg_lowFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
% 
%             if ~isempty(GCvals_low) && ~isempty(subjIDs_low)
%                 rows_low = (subjIDs_low == i_sub);
%                 GC_low = GCvals_low(rows_low, :);
%             else
%                 GC_low = NaN;
%             end
% 
%             % HIGH
%             GCvals_high = Gavg_tempGC_peranatreg_highFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
%             subjIDs_high = Gavg_tempGC_peranatreg_highFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
% 
%             if ~isempty(GCvals_high) && ~isempty(subjIDs_high)
%                 rows_high = (subjIDs_high == i_sub);
%                 GC_high = GCvals_high(rows_high, :);
%             else
%                 GC_high = NaN;
%             end
% 
%            if ~isempty(GC_low) && ~isempty(GC_high)
%                 n_rows = min(size(GC_low,1), size(GC_high,1));
%                 GC_diff = GC_high(1:n_rows, :) - GC_low(1:n_rows, :);  % pairwise subtraction
%                 pooled_diffs = [pooled_diffs; GC_diff];  % append
%            end
%         end
%         GC_Par_to_Front_HighMinusLow.(effect_field).(sprintf('TD_%d', i_TD)) = pooled_diffs;
%     end
% end
% 
% % Average across TDs
% GC_all = cell(1, length(effectTypeNames));
% pvals = zeros(1, length(effectTypeNames));
% numEffects = numel(effectTypeNames);
% 
% % Prepare data
% for i = 1:length(effectTypeNames)
%     data_TD1 = GC_Par_to_Front_HighMinusLow.(effectTypeNames{i}).TD_1;
%     data_TD2 = GC_Par_to_Front_HighMinusLow.(effectTypeNames{i}).TD_2;
%     
%     % Average across TDs
%     avgData = mean([data_TD1(:), data_TD2(:)], 2, 'omitnan');
%     GC_all{i} = avgData;
%     
%     % t-test against zero
%     [~, p] = ttest(avgData);
%     pvals(i) = p;
% end
% 
% % Plot
% effectNames_plot = {'Pred --> Pred', 'PE --> PE', 'Pred --> PE', 'PE --> Pred'};
% figure; hold on; box on;
% colors = lines(numEffects);  % distinct colors for each cloud
% rng(0);  % for reproducible jitter
% 
% for i = 1:numEffects
%     x = i + (rand(size(GC_all{i})) - 0.5) * 0.4;  % add jitter
%     scatter(x, GC_all{i}, 25, ...
%         'MarkerFaceColor', colors(i,:), ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceAlpha', 0.6, ...
%         'MarkerEdgeAlpha', 0.6);
% end
% 
% % Formatting
% xlim([0.5, numEffects + 0.5]);
% xticks(1:numEffects);
% xticklabels(effectNames_plot);
% xtickangle(20);
% ylabel('GC Difference (High FTPL − Low FTPL)');
% yline(0, '--k', 'LineWidth', 1);
% 
% % Add exact p-values above each scatter cloud
% yl = ylim;
% extra_space = 0.1 * range(yl);  % 10% extra space
% ylim([yl(1) + 2 * extra_space, yl(2) + extra_space]);
% 
% % Update y-limit variable for use in placing the p-values
% yl = ylim;
% text_y = yl(2) - 0.05 * range(yl);
% for i = 1:numEffects
%     text(i, text_y, sprintf('p = %.3g', pvals(i)), ...
%         'HorizontalAlignment', 'center', ...
%         'FontSize', 11, 'FontWeight', 'bold');
% end
% 
% % Save figure
% set(gcf, 'Color', 'w');
% 
% 
% filename = ['boxplot_compare_GC_lowhigh_FTPL_Par_to_Front_allpairs.png'];
% figfile = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figures/' filename];
% saveas(gcf, figfile, 'png');
% 

%% Computing stats for FTPL low vs high GC comparison across pairs from different effect-type pairings using AVERAGED data
% Compute GC low vs high FTPL for Parietal > Frontal
AnatReg_CatLabels = ... %Select which anatomical regions should be plotted
    {'AntPFC'; ...
    'VentralT'; ...
    'SupParLob'};
sub_list_FTPL = sub_list(2:9);

Gavg_tempGC_peranatreg_lowFTPL = NASTD_ECoG_Connectivity_ObtainGCByAnatRegions ...
    (sub_list_FTPL, ...
    GCdata_lowFTPL.GCdata_lowFTPL(2:9), ...
    ToneDur_text, ...
    SelElecs, AnatReg_CatLabels, ...
    save_poststepFigs, ...
    param, paths_NASTD_ECoG); %This will have the format, for each field, of i_effect cells, with format nSourceElecs x nTargetElecs x nTD

Gavg_tempGC_peranatreg_highFTPL = NASTD_ECoG_Connectivity_ObtainGCByAnatRegions ...
    (sub_list_FTPL, ...
    GCdata_highFTPL.GCdata_highFTPL(2:9), ...
    ToneDur_text, ...
    SelElecs, AnatReg_CatLabels, ...
    save_poststepFigs, ...
    param, paths_NASTD_ECoG);  %This will have the format, for each field, of i_effect cells, with format nSourceElecs x nTargetElecs x nTD

index_of_interest_source = 7; %Parietal
index_of_interest_target = 1; %Frontal
effectTypeNames = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};

GC_Par_to_Front_HighMinusLow = struct();  % Final output

for i_effect = 1:length(GCdata_lowFTPL.GCdata_lowFTPL{2}.temporalGC)
    effect_field = effectTypeNames{i_effect};
    for i_TD = 1:length(GCdata_lowFTPL.GCdata_lowFTPL{2}.temporalGC{1})
        if ~isfield(GC_Par_to_Front_HighMinusLow, effect_field)
            GC_Par_to_Front_HighMinusLow.(effect_field) = struct();
        end
        if ~isfield(GC_Par_to_Front_HighMinusLow.(effect_field), sprintf('TD_%d', i_TD))
            GC_Par_to_Front_HighMinusLow.(effect_field).(sprintf('TD_%d', i_TD)) = [];
        end

        pooled_diffs = [];  % temporary holder for this TD across all subjects

        for i_sub = 1:8
            % LOW
            GCvals_low = Gavg_tempGC_peranatreg_lowFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
            subjIDs_low = Gavg_tempGC_peranatreg_lowFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};

            if ~isempty(GCvals_low) && ~isempty(subjIDs_low)
                rows_low = (subjIDs_low == i_sub);
                GC_low = GCvals_low(rows_low, :);
            else
                GC_low = NaN;
            end

            % HIGH
            GCvals_high = Gavg_tempGC_peranatreg_highFTPL.source2target{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};
            subjIDs_high = Gavg_tempGC_peranatreg_highFTPL.subject_index{i_effect}{index_of_interest_source, index_of_interest_target, i_TD};

            if ~isempty(GCvals_high) && ~isempty(subjIDs_high)
                rows_high = (subjIDs_high == i_sub);
                GC_high = GCvals_high(rows_high, :);
            else
                GC_high = NaN;
            end

           if ~isempty(GC_low) && ~isempty(GC_high)
                n_rows = min(size(GC_low,1), size(GC_high,1));
                GC_diff = GC_high(1:n_rows, :) - GC_low(1:n_rows, :);  % pairwise subtraction
                pooled_diffs = [pooled_diffs; GC_diff];  % append
           end
        end
        GC_Par_to_Front_HighMinusLow.(effect_field).(sprintf('TD_%d', i_TD)) = pooled_diffs;
    end
end

% Average across TDs
GC_all = cell(1, length(effectTypeNames));
pvals = zeros(1, length(effectTypeNames));
numEffects = numel(effectTypeNames);

% Prepare data
for i = 1:length(effectTypeNames)
    data_TD1 = GC_Par_to_Front_HighMinusLow.(effectTypeNames{i}).TD_1;
    data_TD2 = GC_Par_to_Front_HighMinusLow.(effectTypeNames{i}).TD_2;
    
    % Average across TDs
    minLen = min(numel(data_TD1), numel(data_TD2));
    avgData = mean([data_TD1(1:minLen), data_TD2(1:minLen)], 2, 'omitnan');
    GC_all{i} = avgData;
    
    % t-test against zero
    [~, p] = ttest(avgData);
    pvals(i) = p;
end

% Plot
effectNames_plot = {'Pred --> Pred', 'PE --> PE', 'Pred --> PE', 'PE --> Pred'};
figure; hold on; box on;
colors = lines(numEffects);  
rng(0);  

% Prepare data for boxplot
allData = [];
groupData = [];
for i = 1:numEffects
    allData = [allData; GC_all{i}(:)];
    groupData = [groupData; repmat(i, numel(GC_all{i}), 1)];
end

% Overlay scatter points (lower alpha)
for i = 1:numEffects
    x = i + (rand(size(GC_all{i})) - 0.5) * 0.4;
    scatter(x, GC_all{i}, 25, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.3, ...
        'MarkerEdgeAlpha', 0.3);
    
    % Mean line
    m = mean(GC_all{i}, 'omitnan');
    plot([i-0.3, i+0.3], [m, m], 'Color', colors(i,:), 'LineWidth', 3);
end

% Formatting
xlim([0.5, numEffects + 0.5]);
xticks(1:numEffects);
xticklabels(effectNames_plot);
xtickangle(20);
ylabel('GC Difference (High FTPL > Low FTPL)');
yline(0, '--k', 'LineWidth', 1);

% Add exact p-values or significance stars
ylim([-0.01, 0.02]);
yl = ylim;
text_y = yl(2) - 0.05 * range(yl);
for i = 1:numEffects
    if pvals(i) < 0.001
        sigLabel = '***';
    elseif pvals(i) < 0.01
        sigLabel = '**';
    elseif pvals(i) < 0.05
        sigLabel = '*';
    else
        sigLabel = sprintf('p=%.3g', pvals(i));
    end
    text(i, text_y, sigLabel, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
end

set(gcf, 'Color', 'w');

filename = ['boxplot_compare_GC_lowhigh_FTPL_Par_to_Front_allpairs_avgData.png'];
figfile = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figures/' filename];
saveas(gcf, figfile, 'png');

%% Plot GC results
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
if strcmp(param.GC.tone_aggregation_label,'_30_31_32_33_34')
    NASTD_ECoG_Connectivity_PlotGC_statp30top32_vs_p33p34_lk ...
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
%%%load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr_HighGamma.mat') %p < 0.05 thresh,uncorr elec selection    
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


%% 4) Plot and do statistics for electrode combinations between regions of the auditory hierarchy

addpath('/isilon/LFMI/VMdrive/Lua/Temp2A/Temp2A_fMRI_pilot/toolboxes/spm12');
load('/isilon/LFMI/VMdrive/Thomas//NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Analysis/Connectivity/ElecSelect/Allsub_n9/SelElecs_p0.05uncorr.mat') %p < 0.05 thresh,uncorr elec selection    
atlas_path = 'brodmann.nii'; 

% Determine if any of the electrodes are located in one of the auditory
% regions of interest
% Classify each electrode as part of a Brodmann's area
BA_areas_Pred = assign_BA_from_MNI(SelElecs.PredEffect.("Elec Coords"), atlas_path);
BA_areas_PE = assign_BA_from_MNI(SelElecs.PEEffect.("Elec Coords"), atlas_path);

% Define ROI map: Brodmann areas → auditory/cognitive labels
roi_map = containers.Map( ...
    [41, 42, 44, 45, 8, 40, 6, 22], ...
    {'Auditory Cortex', ...
     'Auditory Cortex', ...
     'IFC Dorsal', ...
     'IFC Ventral', ...
     'dlPFC', ...
     'IPL', ...
     'PMC', ...
     'STS'});

% For Prediction electrodes
roi_labels_pred = cell(size(BA_areas_Pred));
for i = 1:length(BA_areas_Pred)
    roi_labels_pred{i} = get_roi_label(BA_areas_Pred(i), roi_map);
end

% For Prediction Error electrodes
roi_labels_PE = cell(size(BA_areas_PE));
for i = 1:length(BA_areas_PE)
    roi_labels_PE{i} = get_roi_label(BA_areas_PE(i), roi_map);
end

% Combine everything into one table
SelElecs.PredEffect.BAlabel = roi_labels_pred;
SelElecs.PEEffect.BAlabel = roi_labels_PE;

% Also need to update SelElecs.Pairs for stats
% === Update Pred_Pred with PredEffect ===
predPairFields = fieldnames(SelElecs.Pairs.Pred_Pred.AllRegions.Pairs_perelec);
predEffect = SelElecs.PredEffect;

for i = 1:length(predPairFields)
    pairName = predPairFields{i};
    targetElecs = SelElecs.Pairs.Pred_Pred.AllRegions.Pairs_perelec.(pairName);
    BAlabels = cell(height(targetElecs), 1);
    for j = 1:height(targetElecs)
        elecLabel = targetElecs.("Electrode Label"){j};
        subjLabel = targetElecs.("Subject Label"){j};

        % Find matching row in PredEffect
        matchIdx = strcmp(predEffect{:,1}, elecLabel) & strcmp(predEffect{:,2}, subjLabel);

        if any(matchIdx)
            baLabel = predEffect{matchIdx,12}{1};  % Get BAlabel (12th column)
        else
            baLabel = 'NA';
        end
        
        BAlabels{j,1} = baLabel;
    end
    targetElecs.BAlabel = BAlabels;  % Add as new column
    SelElecs.Pairs.Pred_Pred.AllRegions.Pairs_perelec.(pairName) = targetElecs;  % Update
end

% === Update PE_PE with PEEffect ===
pePairFields = fieldnames(SelElecs.Pairs.PE_PE.AllRegions.Pairs_perelec);
peEffect = SelElecs.PEEffect;

for i = 1:length(pePairFields)
    pairName = pePairFields{i};
    targetElecs = SelElecs.Pairs.PE_PE.AllRegions.Pairs_perelec.(pairName);
    BAlabels = cell(height(targetElecs), 1);
    for j = 1:height(targetElecs)
        elecLabel = targetElecs.("Electrode Label"){j};
        subjLabel = targetElecs.("Subject Label"){j};

        % Find matching row in PEEffect
        matchIdx = strcmp(peEffect{:,1}, elecLabel) & strcmp(peEffect{:,2}, subjLabel);

        if any(matchIdx)
            baLabel = peEffect{matchIdx,12}{1};  % Get BAlabel (12th column)
        else
            baLabel = 'NA';
        end
        
       BAlabels{j,1} = baLabel;
    end
    targetElecs.BAlabel = BAlabels;  % Add as new column
    SelElecs.Pairs.PE_PE.AllRegions.Pairs_perelec.(pairName) = targetElecs;  % Update
end

% Output a surface project of Brodmanns areas with only the desired ROIs
% addpath '/isilon/LFMI/VMdrive/Lua/surfice_atlas-master'
% desired_ROIs = [41, 42, 44, 45, 8, 40, 6, 22];
% nii_nii2atlas_rois('brodmann.nii', 'mylut.lut', desired_ROIs);

%% Extract GC values for each pair 
subs                    = sub_list(vars.validSubjs);
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

%param.GC.ElecPairEffect   = {'Pred_Pred','PE_PE','Pred_PE','PE_Pred'};
param.GC.ElecPairEffect   = {'Pred_Pred','PE_PE'};
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

param.GC.tone_aggregation = {1, 6, 11, 16, 21, 26, 31};
param.GC.tone_aggregation = {2, 7, 12, 17, 22, 27, 32};
param.GC.tone_aggregation = {30, 31, 32, 33, 34};

param.GC.tone_aggregation_label = cellfun(@(c)[c],param.GC.tone_aggregation);
param.GC.tone_aggregation_label = num2str((param.GC.tone_aggregation_label),'_%d');
param.GC.tone_aggregation_label = erase(param.GC.tone_aggregation_label, ' ');

param.GC.epochdur_ms = 'full'; %'full', 200, 100
if strcmp(param.GC.epochdur_ms, 'full')
    param.GC.epochdur_ms_label = 'fullTW';
else
   param.GC.epochdur_ms_label = [num2str(param.GC.epochdur_ms) 'msTW'];
end
param.GC.ElecPairRegion   = 'AllRegions';
param.GC.InputDataType    = {'Broadband'}; %{'Broadband', 'HP05toLP30Hz', 'HighGamma_LogAmp'}; 

% Load in GC data
path_outputdata = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/']);
if (~exist(path_outputdata, 'dir')); mkdir(path_outputdata); end
label_outputdata = ['GCdata_n' num2str(length(vars.validSubjs)) '_' ...
    param.GC.ElecPairRegion '_' ...
    param.GC.InputDataType{1} param.GC.tone_aggregation_label '_' ...
    param.GC.epochdur_ms_label '_ensemnorm'];
load([path_outputdata label_outputdata], 'GCdata')

%Add in BA label information to GCdata
% for i_sub = vars.validSubjs
%     for i_effect = 1:2
%         % Assign labels that correspond to ind_pairedelecs_from
%         sub = subs(i_sub);
%         if i_effect == 1
%             effecttype = 'Pred';
%         else 
%             effecttype = 'PE';
%         end
%         GCdata{i_sub}.BAlabels{i_effect} = MNI_coords_all.BAlabel(strcmp(MNI_coords_all.("Subject Label"), sub) & strcmp(MNI_coords_all.EffectType, effecttype),:);
%     end
% end

% Now the GC results also have the BA labels for auditory cortex

%Plot GC results
AnatReg_CatLabels = {'Auditory Cortex', ...  %Select which anatomical regions should be plotted
    'IFC Dorsal', ...
    'IFC Ventral', ...
    'dlPFC', ...
    'IPL', ...
    'PMC', ...
    'STS'};

%Extract GC per anat regions
NASTD_ECoG_Connectivity_BAregions_lk ...
        (sub_list(vars.validSubjs), GCdata,...
        ToneDur_text, SelElecs, AnatReg_CatLabels, ...
        save_poststepFigs, ...
        param, paths_NASTD_ECoG);
    
% GCdata = cell(1, 9);  % Initialize output
% 
% for i = 1:9
%     % Start from GCdata_T1 and overwrite the field we want to average
%     GCdata{i} = GCdata_T1{i};
% 
%     avg_tgc = cell(1, 4);  % Correct shape: 1×4
% 
%     for j = 1:4
%         avg_tgc{j} = cell(1, 2);  % Each is a 1×2 cell
%         
%         for k = 1:2
%             % Extract source2target matrices
%             mat1 = GCdata_T1{i}.temporalGC{j}{k}.source2target;
%             mat2 = GCdata_T2{i}.temporalGC{j}{k}.source2target;
% 
%             % Compute average
%             avg_mat = (mat1 + mat2) / 2;
% 
%             % Copy one struct and overwrite source2target with the average
%             avg_struct = GCdata_T1{i}.temporalGC{j}{k};
%             avg_struct.source2target = avg_mat;
% 
%             % Save back
%             avg_tgc{j}{k} = avg_struct;
%         end
%     end
% 
%     % Assign averaged temporalGC to output
%     GCdata{i}.temporalGC = avg_tgc;
% end