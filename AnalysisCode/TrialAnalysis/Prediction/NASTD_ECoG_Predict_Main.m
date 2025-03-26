%% Project: NASTD_ECoG
%Compute prediction and prediciton error effects
%1) prediction effect (linear regression for neural activity during p33 on p*34)
%2) prediction error effect (linear regression for neural activity during p34 on p34)

%1. Load in pre-processed data
%2. Prepare input data (ERF-range and amplitude envelopes)
%3. Compute Prediction/Prediction error effects
%3a: Windowed approach (Neural data averaged across 50 ms time windows)
%3b: Sample-based approach (Effects computed for each sample and
%corrected across samples)
%4. Plot output
%4a: PredEffect ERF (p33 per p*34)
%4b: PredErrorEffect ERF (p34 per p34-p*34/p33-p34 pitch distance)

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
save_poststepFigs = 0;

%% 1. Load in pre-processed data
for i_sub = vars.validSubjs
    
    sub = sub_list{i_sub};
    disp(['Processing subject: ' sub])
    
    %     prediction_dir = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/'];
    %     if (~exist(prediction_dir, 'dir')); mkdir(prediction_dir); end
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
    % 0.4) Specify analysis info
    %     Selected_Elecs = ...
    %         setdiff(DataClean_AllTrials.cfg.info_elec.selected.index4EDF, ...
    %         subs_PreProcSettings.(sub).rejectedChan_index);
    SampleFreq = DataClean_AllTrials.fsample;
    
    %% 2) Prepare input data
    %2.1) Select only clean/valid trials
    DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);
    
    %2.2) Select only valid electrodes (Grid+Strip, clean, MNI-coordinates present)
    Index_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.index4EDF, ...
        subs_PreProcSettings.(sub).rejectedChan_index); %via index
    Label_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.Label, ...
        subs_PreProcSettings.(sub).rejectedChan_label); %via label
    %Check via label
    if ~isempty(find((strcmp(...
            sort(DataClean_CleanTrials.label(Index_Valid_Elecs)), ...
            sort(Label_Valid_Elecs))) == 0))
        disp(['Mismatch in electrode selection']);
    end
    
    %    cd([paths_NASTD_ECoG.FieldTrip 'utilities/']) %In case ft_selectdata not found
    cfg         = [];
    cfg.channel = (Label_Valid_Elecs);
    DataClean_CleanTrialsElecs = ft_selectdata(cfg, DataClean_CleanTrials);
    
    %2.3) Select only trials with specific Tone Duration (TD)
    for i_TD = 1:length(ToneDur_text)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_Preproc_SelTrialsperTD...
            (ToneDur_text{i_TD}, DataClean_CleanTrialsElecs);
        
        disp(['Total trial number for TD ' ToneDur_text{i_TD} 's: ' ...
            num2str(length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial))])
        
        %Adjust info subfield
        DataClean_CleanTrialsElecs_perTD{i_TD}.info.elec = ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_elec;
        DataClean_CleanTrialsElecs_perTD{i_TD}.info.trigger = ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_trigger;
        DataClean_CleanTrialsElecs_perTD{i_TD}.info.ref = ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_ref;
        
        DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous = ...
            rmfield(DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous, 'info_elec');
        DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous = ...
            rmfield(DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous, 'info_trigger');
        DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous = ...
            rmfield(DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous, 'info_ref');
    end
    
    %Cleanup
    clear DataClean_CleanTrialsElecs DataClean_CleanTrials DataClean_AllTrials
    
    %% 3) Prepare input signals: slow arrhythmic activity & gamma amplitude envelope
    
    %3.0A Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus (i.e., tone sequence start) baseline
    %Note: Prestimulus baseline = Fixation button shown prior to sequence start for ~0.7s
    %End trial x-1 to start trial x = 0.001s
    %Start trial x to Onset Fixation Dot = 0.5s
    %Onset fixation dot to Sequence start = ~0.7s
    
    BL_win = [-0.5 0];
    
    for i_TD = 1:length(ToneDur_text)
        for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
            
            for i_elec = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.label)
                DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) = ...
                    DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) - ...
                    nanmean(DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,...
                    find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == BL_win(1)): ...
                    find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == BL_win(2))));
            end
            
        end
    end
    
    %3.0B Remove trials where some electrodes have only NaN samples (NY798)
    for i_TD = 1:length(ToneDur_text)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
            (sub, ToneDur_text{i_TD}, DataClean_CleanTrialsElecs_perTD{i_TD});
    end
    
    %3.1 Apply Filter for Slow ERF-range activity and deal with NaNs
%     for i_TD = 1:length(ToneDur_text)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.LP35Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.LP35Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_LP35Hz...
%             (sub, ToneDur_text{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end
%     for i_TD = 1:length(ToneDur_text)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.HP01toLP30Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP01toLP30Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_HP01toLP30Hz...
%             (sub, ToneDur_text{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.HP05toLP30Hz, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    end
%     for i_TD = 1:length(ToneDur_text)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.HP1toLP30Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP1toLP30Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_HP1toLP30Hz...
%             (sub, ToneDur_text{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end
%     for i_TD = 1:length(ToneDur_text)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.HP2toLP30Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP2toLP30Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_HP2toLP30Hz...
%             (sub, ToneDur_text{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end
    
    %     %3.2 Compute Frequency-band specific amplitude envelope and deal with NaNs
    %     FrequencyBands = [8 12; 15 30; 30, 70; 70, 150]; %Hz [HP,LP]
    %     FrequencyBand_Labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
    %
    %     for i_TD = 1:length(ToneDur_text)
    %         for i_freqbands = 1:length(FrequencyBands)
    %
    %             FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
    %
    %             [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
    %                 DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
    %                 NASTD_ECoG_FiltNaNinterp_AmpEnvel...
    %                 (sub, ToneDur_text{i_TD}, ...
    %                 FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel,...
    %                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
    %                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    %
    %         end
    %     end
    
    %To do: Add a sample*channel plot showing how many trials are going
    %present for each samples. This will sho whow much data is going into
    %stat and avg
    
    
    %% 4) Compute Prediction & Prediction Error effects
    %4.0 Parameters:
    %     InputDataType = {'HP2toLP30Hz', 'HP05toLP30Hz', 'LP35Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
    InputDataType = {'HP05toLP30Hz'};
    
    param.PredSeqrange = 33; %final tone defining range for tone pitches based on which prediction should be computed (i.e., 1-33)
    param.ToneIndex = 33; %tone index indicating for which tone ERF/ERP should be computed, which is then used in regression analysis of prediction effect
    
    param.SamplesTW = 25; %25 = 50ms TW
    param.Label_TW = num2str(round(param.SamplesTW/SampleFreq,2));
    
    param.alpha = 0.025; %i.e., 2-sided test at p < 0.05
    param.clusteralpha = 0.05;
    param.numreps = 1000;%Number of repetitions to compute Montecarlo-based cluster statistics
    
    % 4.1 Compute prediction & prediction error effect (Windowed, uncorrected across time)
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(InputDataType)
            
            NASTD_ECoG_Predict_CompPred_TW_Subs...
                (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
            
        end
    end
    
    % 4.2 Compute prediction & prediction error effect (sample-wise, cluster-corrected across time)
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(InputDataType)
            
            NASTD_ECoG_Predict_CompPred_Sample_Subs...
                (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
            
        end
    end
    
    %4.3A Plot prediction effect (neural signals per predp34 during p33) and mark sign. TWs and samples
    for i_TD = 1:length(ToneDur_text)
        for i_inputData = 1:length(InputDataType)
            NASTD_ECoG_Predict_ComparePredEffects_TWvsSample...
                (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                save_poststepFigs, paths_NASTD_ECoG);
        end
    end
    
    %4.3B Plot prediction error effect (neural signals per predp34-p34 distance during p34) and mark sign. TWs and samples
    InputEffectType = {'SimplePredErrEffect','ComplexPredErrEffect'};
    for i_effect = 1:length(InputEffectType)
        for i_TD = 1:length(ToneDur_text)
            for i_inputData = 1:length(InputDataType)
                NASTD_ECoG_Predict_ComparePredErrEffects_TWvsSample...
                    (sub, InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                    param,...
                    DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                    save_poststepFigs, paths_NASTD_ECoG);
            end
        end
    end
    
    %4.3C Plot Summmary plot showing electrode location on surface brain,
    %and ERF timelines + sign. samples for sign. electrodes
    InputEffectType = {'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};
    for i_effect = 1:length(InputEffectType)
        for i_TD = 1:length(ToneDur_text)
            for i_inputData = 1:length(InputDataType)
                NASTD_ECoG_Predict_PlotEffectSummary_Sub...
                    (sub, InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                    param,...
                    DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                    save_poststepFigs, paths_NASTD_ECoG);
            end
        end
    end
end

%% 5) Aggregate prediction effects and FDR-Correct for multiple comparisons across electrodes

%5.0) input parameters
% InputDataType = {'HP2toLP30Hz', 'HP05toLP30Hz', 'LP35Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
InputDataType = {'HP05toLP30Hz'};
InputEffectType =  {'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};

SampleFreq      = 512;

param.alpha     = 0.025;
param.alpha_FDR = 0.05;
param.SamplesTW = 25; %25 = 50ms TW
param.Label_TW  = num2str(round(param.SamplesTW/SampleFreq,2));

%5.1 Compute FDR correction, aggregate effects and plot summary figure
%showing sign. elecs across freqs
SignElecIndices = NASTD_ECoG_Predict_PredEffectTable...
    (vars.validSubjs, sub_list, InputDataType, InputEffectType, ToneDur_text, ...
    param,...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);


%% 6) Plot electrodes showing sign. prediction effect
param.pval_plotting     = 0.025; %Pval thresh for plotting
param.FDRcorrect        = 0;
param.pval_FDR          = 0.05;
param.plot_SubplotperTW = 1; %Common plot across TW
param.ElecSelect        = 'All'; %StimCorr, All
%  set(0, 'DefaultFigureVisible', 'off') %show surface plots

param.SamplesTW         = 25; %25 = 50ms TW
param.Label_TW          = num2str(round(param.SamplesTW/512,2));
param.ToneIndex         = 33;

InputDataType           = {'HP05toLP30Hz', 'HighGamma_LogAmp'};
%     {'LP35Hz', ...
%     'HP05toLP30Hz', 'HP1toLP30Hz',  'HP2toLP30Hz', ...
%     'Alpha_LogAmp', 'Beta_LogAmp', ...
%     'Gamma_LogAmp', 'HighGamma_LogAmp'};
InputEffectType         = {'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};
ToneDur_text            = {'0.2' '0.4'};

% %6.1 Single-subject Volume for time-windowed results
% for i_pval = 1:length(param.pval_plotting)
%     for i_sub = vars.validSubjs
%         for i_effect = 1:length(InputEffectType)
%             for i_inputData = 1:length(InputDataType)
%                 for i_TD = 1:length(ToneDur_text) %tone duration condition
%
%                     NASTD_ECoG_Predict_PlotSignElec_Sub...
%                         (sub_list{i_sub}, ...
%                         InputDataType{i_inputData}, InputEffectType{i_effect}, ToneDur_text{i_TD}, ...
%                         param,...
%                         save_poststepFigs, paths_NASTD_ECoG)
%
%                 end
%             end
%         end
%     end
% end

% %6.2 All subjects on one volume for time-windowed results
% subs = sub_list(vars.validSubjs);
%
% for i_pval = 1:length(param.pval_plotting) %p-value for plotting
%     for i_effect = 1:length(InputEffectType)
%         for i_inputData = 1:length(InputDataType)
%             for i_TD = 1:length(ToneDur_text) %tone duration condition
%
%                 NASTD_ECoG_Predict_PlotSignElec_AllSub...
%                     (subs, ...
%                     InputDataType{i_inputData}, InputEffectType{i_effect}, ToneDur_text{i_TD}, ...
%                     param,...
%                     save_poststepFigs, paths_NASTD_ECoG)
%
%             end
%         end
%     end
% end

%6.3.1 All subjects on one volume for sample-wise results per TD
subs = sub_list(vars.validSubjs);
% subs = sub_list(2:9);

for i_effect = 1:length(InputEffectType)
    for i_inputData = 1:length(InputDataType)
        for i_TD = 1:length(ToneDur_text) %tone duration condition
            
            %%Plot specific parameter (i.e., clusterstat, p-val, onset, offset, effect duration)
            %             for i_param2plot = 1:length(param2plot) %1 plot per parameter
            %                 param.param2plot = param2plot{i_param2plot};
            %                 NASTD_ECoG_Predict_PlotSignClusterElec_AllSub...
            %                     (subs, ...
            %                     InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
            %                     param,...
            %                     save_poststepFigs, paths_NASTD_ECoG)
            %             end
            
            %Plot all parameters on a summary plot
            NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSub...
                (subs, ...
                InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                save_poststepFigs, paths_NASTD_ECoG)
            
            
            %Plot summary figure showing 
            %1) traces per p*34 condition during p33 for all sign. elecs
            %2) traces per p*34 condition during p1-p33 for all sign. elecs
            %3) surface plot overview
            %for all subjects, with traces grouped per anat region
            if strcmp(InputEffectType{i_effect}, 'PredEffect')
                %Note: Makes only sense for pred effect, since p*34 distinction
                %is not useful for PE
                NASTD_ECoG_Predict_PlotERF4SignClusterElec_AllSub...
                    (subs, ...
                    InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                    param,...
                    save_poststepFigs, paths_NASTD_ECoG)
            end            
            
        end
    end
end

%6.3.2 All subjects on one volume for sample-wise results aggregated across TD
subs = sub_list(vars.validSubjs);
% subs = sub_list(2:9);

for i_effect = 1:length(InputEffectType)
    for i_inputData = 1:length(InputDataType)
        
        NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSubTD...
            (subs, ...
            InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
            param,...
            save_poststepFigs, paths_NASTD_ECoG)
%         
%         %figure output saved as movie
%         NASTD_ECoG_Predict_PlotFDRSignClusterElec_AllSubTD_Mov...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)
        
        %Plot summary figure showing
        %1) traces per p*34 condition during p33 for all sign. elecs
        %2) traces per p*34 condition during p1-p33 for all sign. elecs
        %3) surface plot overview
        %for all subjects, with traces grouped per anat region
        if strcmp(InputEffectType{i_effect}, 'PredEffect')
            %Note: Makes only sense for pred effect, since p*34 distinction
            %is not useful for PE
            NASTD_ECoG_Predict_PlotERF4SignClusterElec_AllSubTD...
                (subs, ...
                InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
                param,...
                save_poststepFigs, paths_NASTD_ECoG)
        end
    end
end

%6.4 All subjects on one volume for sample-wise results per TD with FDR-correction per lobe/anatomical parcel
subs = sub_list(vars.validSubjs);
% subs = sub_list(2:9);

for i_effect = 1:length(InputEffectType)
    for i_inputData = 1:length(InputDataType)
        
        NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_RegionalFDR...
            (subs, ...
            InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
            param,...
            save_poststepFigs, paths_NASTD_ECoG)
        
%To do: ERF for regionalFDR corrected prediction effects
        
    end
end
