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

plot_poststepFigs = 0;
save_poststepFigs = 0;

FTPL_ratings = load([paths_NASTD_ECoG.Analysis_Behavior 'FTPLratings/Allsub_n8/Data/', 'Trialwise_FTPL_allsubs.mat']);
FTPL_ratings = FTPL_ratings.Trialwise_FTPL;

%% 1. Load in pre-processed data
DataAll = cell(length(sub_list),1); 
for i_sub = vars.validSubjs
    
    plot_poststepFigs = 0;

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
    %         FrequencyBands = [8 12; 15 30; 30, 70; 70, 150]; %Hz [HP,LP]
    %         FrequencyBand_Labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
    FrequencyBands = [70, 150; 1 2]; %Hz [HP,LP]
    FrequencyBand_Labels = {'HighGamma'};
    %         FrequencyBands = [8 12; 15 30]; %Hz [HP,LP]
    %         FrequencyBand_Labels = {'Alpha', 'Beta'};
    %         FrequencyBands = [70, 150; 30, 70];
    %         FrequencyBand_Labels = {'HighGamma','Gamma'};
    for i_TD = 1:length(ToneDur_text)
        for i_freqbands = 1:(length(FrequencyBands) -1)
            
            FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
            
            [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
                DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, ToneDur_text{i_TD}, ...
                FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
            
        end
    end
    DataAll{i_sub} = DataClean_CleanTrialsElecs_perTD;

        %% July 2025: New analysis running prediciton and PE analyses at repetition rates to control for entrainment effects
%         FrequencyBands = [4.5 5.5; 2 3]; %Hz [5 Hz = 200 ms TD, 2.5 Hz = 400 ms TD]
%         FrequencyBand_Labels = {'RR_200ms_TD', 'RR_400ms_TD'}; % RR = repetition rate
%         for i_TD = 1:length(ToneDur_text)
%                 FieldLabel = [FrequencyBand_Labels{i_TD} '_LogAmp'];
%     
%                 [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
%                     DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
%                     NASTD_ECoG_FiltNaNinterp_AmpEnvel...
%                     (sub, ToneDur_text{i_TD}, ...
%                     FrequencyBands(i_TD,1), FrequencyBands(i_TD,2), FieldLabel,...
%                     DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                     plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%         end

    
    %To do: Add a sample*channel plot showing how many trials are going
    %present for each samples. This will sho whow much data is going into
    %stat and avg
    
    
    %% December 2025: Compute power spectra for amplitude envelopes of four main frequency bands to show they fluctuate at slow time scales
%     fs = DataClean_CleanTrialsElecs_perTD{1}.fsample;
%     nBands = length(FrequencyBands);
% 
%     for i_TD = 1:length(ToneDur_text)
%         for i_freqbands = 1:nBands
% 
%             FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
%             filtData = DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel);
% 
%             nTrials = numel(filtData);
%             nCh = size(filtData{1},1);
% 
%             % --- Determine minimum samples across trials ---
%             minSamples = min(cellfun(@(x) size(x,2), filtData));
% 
%             % --- Preallocate ---
%             Pxx_all = [];
% 
%             for ch = 1:nCh
%                 for t = 1:nTrials
% 
%                     X = filtData{t}(ch, 1:minSamples);
%                     X = X - mean(X);   % remove DC
% 
%                     % --- FFT-based PSD (full epoch) ---
%                     N = length(X);
%                     Xf = fft(X);
%                     Pxx = (abs(Xf).^2) / (fs * N);
% 
%                     % One-sided spectrum
%                     Pxx = Pxx(1:floor(N/2)+1);
%                     f = (0:floor(N/2)) * fs / N;
% 
%                     Pxx_all(:, end+1) = Pxx; %#ok<SAGROW>
%                 end
%             end
% 
%             % --- Average across trials and channels ---
%             Pxx_mean = mean(Pxx_all, 2);
% 
%             PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).freq = f;
%             PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).Pxx = Pxx_mean;
% 
%         end
%     end

    
    %% 4) Compute Prediction & Prediction Error effects
    %4.0 Parameters:
    %     InputDataType = {'HP2toLP30Hz', 'HP05toLP30Hz', 'LP35Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
    % InputDataType = {'HP05toLP30Hz'};
    %InputDataType = {'Alpha_LogAmp','Beta_LogAmp'};
    InputDataType = {'HP05toLP30Hz','HighGamma_LogAmp'};
    %InputDataType = {'RR_200ms_TD_LogAmp', 'RR_400ms_TD_LogAmp'};
    
    param.PredSeqrange = 33; %final tone defining range for tone pitches based on which prediction should be computed (i.e., 1-33)
    param.ToneIndex = 33; %tone index indicating for which tone ERF/ERP should be computed, which is then used in regression analysis of prediction effect
    
    SampleFreq      = 512;
    param.SamplesTW = 25; %25 = 50ms TW
    param.Label_TW = num2str(round(param.SamplesTW/SampleFreq,2));
    
    param.alpha = 0.025; %i.e., 2-sided test at p < 0.05
    param.clusteralpha = 0.05;
    param.numreps = 1000;%Number of repetitions to compute Montecarlo-based cluster statistics
    
    plot_poststepFigs = 0;
    
    % 4.1 Compute prediction & prediction error effect (Windowed, uncorrected across time)
%     for i_TD = 1:length(ToneDur_text)
%         for i_inputData = 1:length(InputDataType)
%             
%             NASTD_ECoG_Predict_CompPred_TW_Subs...
%                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                 param,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%             
%         end
%     end

     % 4.1.B Compute prediction and prediction error effect (windowed,
    % uncorrected across time) for RR control analysis
%     for i_TD = 1:length(ToneDur_text)
%         NASTD_ECoG_Predict_CompPred_TW_Subs...
%             (sub, InputDataType{i_TD}, ToneDur_text{i_TD}, ...
%             param,...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end
    
        % 4.1.C Compute prediction & prediction error effect (Windowed, uncorrected across time)
        % to predict FTLP for control analysis
%         disp(['Processing subject # ' num2str(i_sub) '/9'])
%     for i_TD = 1:length(ToneDur_text)
%         for i_inputData = 1:length(InputDataType)
%             
%             NASTD_ECoG_Predict_CompPred_TW_Subs_FTPL...
%                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                 param,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG, FTPL_ratings(strcmp(FTPL_ratings.SubID, sub) & strcmp(FTPL_ratings.ToneDur, ToneDur_text{i_TD}), :));
%             
%         end
%     end
    
    % 4.2 Compute prediction & prediction error effect (sample-wise, cluster-corrected across time)
%     for i_TD = 1:length(ToneDur_text)
%         for i_inputData = 1:length(InputDataType)
%             
%             NASTD_ECoG_Predict_CompPred_Sample_Subs...
%                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                 param,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%             
%             % Version by LK December 2024 without cluster correction
% %              NASTD_ECoG_Predict_CompPred_Sample_Subs_NoClusterCorr_lk...
% %                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
% %                 param,...
% %                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
% %                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%             
%         end
%     end
    
    % 4.2.B Compute prediction and prediction error effect (sample-wise,
    % for control analysis looking at 'entrainment' and repetition rate of
%     % each TD)
%     for i_TD = 1:length(ToneDur_text)
%         % With cluster correction
%         NASTD_ECoG_Predict_CompPred_Sample_Subs...
%                 (sub, InputDataType{i_TD}, ToneDur_text{i_TD}, ...
%                 param,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%         
%         % Version by LK December 2024 without cluster correction
% %         NASTD_ECoG_Predict_CompPred_Sample_Subs_NoClusterCorr_lk...
% %             (sub, InputDataType{i_TD}, ToneDur_text{i_TD}, ...
% %             param,...
% %             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
% %             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end

     % 4.2.C Compute prediction and prediction error effect (sample-wise,
    % to predict FTPL for control analysis
%     for i_TD = 1:length(ToneDur_text)
%         for i_inputData = 1:length(InputDataType)
%    
%             NASTD_ECoG_Predict_CompPred_Sample_Subs_FTPL...
%                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                 param,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG, FTPL_ratings(strcmp(FTPL_ratings.SubID, sub) & strcmp(FTPL_ratings.ToneDur, ToneDur_text{i_TD}), :));
% 
% %               %No Cluster correction
% %             NASTD_ECoG_Predict_CompPred_Sample_Subs_NoClusterCorr_FTPL...
% %                 (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
% %                 param,...
% %                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
% %                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG, FTPL_ratings(strcmp(FTPL_ratings.SubID, sub) & strcmp(FTPL_ratings.ToneDur, ToneDur_text{i_TD}), :));
% %             
%         end
%     end
    
    
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
%     
%     %4.3B Plot prediction error effect (neural signals per predp34-p34 distance during p34) and mark sign. TWs and samples
%     InputEffectType = {'SimplePredErrEffect','ComplexPredErrEffect'};
%     for i_effect = 1:length(InputEffectType)
%         for i_TD = 1:length(ToneDur_text)
%             for i_inputData = 1:length(InputDataType)
%                 NASTD_ECoG_Predict_ComparePredErrEffects_TWvsSample...
%                     (sub, InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                     param,...
%                     DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                     save_poststepFigs, paths_NASTD_ECoG);
%             end
%         end
%     end
%     
%     %4.3C Plot Summmary plot showing electrode location on surface brain,
%     %and ERF timelines + sign. samples for sign. electrodes
%     InputEffectType = {'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};
%     for i_effect = 1:length(InputEffectType)
%         for i_TD = 1:length(ToneDur_text)
%             for i_inputData = 1:length(InputDataType)
%                 NASTD_ECoG_Predict_PlotEffectSummary_Sub...
%                     (sub, InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                     param,...
%                     DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                     save_poststepFigs, paths_NASTD_ECoG);
%             end
%         end
%     end
end

%% Reviewer request: plot power spectra of amplitude envelope signals

% Save the power spectra
% save('/isilon/LFMI/VMdrive/Lua/NASTD/PowerSpectrum_Subject.mat', 'PowerSpectrum_Subject');
% 
% figure;
% 
% nSubjects = length(fieldnames(PowerSpectrum_Subject)); % assumes each subject is a field
% 
% for i_TD = 1:length(ToneDur_text)
%     subplot(1, length(ToneDur_text), i_TD);
%     hold on
% 
%     for i_freqbands = 1:nBands
%         % Initialize to accumulate power across subjects
%         Pxx_all = [];
%         f_all = [];
%         
%         subjectNames = fieldnames(PowerSpectrum_Subject);
%         for i_sub = 1:nSubjects
%             sub = subjectNames{i_sub};
%             
%             f   = PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).freq;
%             Pxx = PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).Pxx;
% 
%             % Store for averaging
%             if isempty(Pxx_all)
%                 Pxx_all = zeros(nSubjects, length(Pxx));
%                 f_all   = f;  % assume frequency vector is same for all subjects
%             end
%             
%             Pxx_all(i_sub, :) = Pxx;
%         end
% 
%         % Average across subjects
%         meanPxx = mean(Pxx_all, 1);
% 
%         % Plot (skip DC component)
%         loglog(f_all(2:end), meanPxx(2:end), 'LineWidth', 2);
%     end
% 
%     xlabel('Frequency (Hz)')
%     ylabel('Power')
%     title(sprintf('TD = %s', ToneDur_text{i_TD}))
%     xlim([0.1 30])   % log-log: avoid 0
%     set(gca, 'FontSize', 14, 'Box', 'off')
% 
%     if i_TD == length(ToneDur_text)
%         legend(FrequencyBand_Labels, 'Location', 'northeast')
%     end
% end
% set(gca,'FontSize',14);
% saveas(gcf, fullfile('/isilon/LFMI/VMdrive/Lua/NASTD/', 'Group_PSD_LogPower_AmpEnvel.png'));   % PNG


%% 5) Aggregate prediction effects and FDR-Correct for multiple comparisons across electrodes

%5.0) input parameters
% InputDataType = {'HP2toLP30Hz', 'HP05toLP30Hz', 'LP35Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
%InputDataType = {'Alpha_LogAmp','Beta_LogAmp'};
InputDataType = {'HP05toLP30Hz','HighGamma_LogAmp'};
%InputDataType = {'RR_200ms_TD_LogAmp', 'RR_400ms_TD_LogAmp'};
%InputEffectType =  {'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};
InputEffectType =  {'PredEffect', 'ComplexPredErrEffect'};
%InputEffectType = {'PredErrFTPL'}; %This is for the FTPL control analysis
%InputEffectType = {'ComplexPredErrEffect'};

SampleFreq      = 512;

param.alpha     = 0.025;
param.alpha_FDR = 0.05;
param.SamplesTW = 25; %25 = 50ms TW
param.Label_TW  = num2str(round(param.SamplesTW/SampleFreq,2));
plot_poststepFigs = 0;

%5.1 Compute FDR correction, aggregate effects and plot summary figure
%showing sign. elecs across freqs
SignElecIndices = NASTD_ECoG_Predict_PredEffectTable_lk...
    (vars.validSubjs, sub_list, InputDataType, InputEffectType, ToneDur_text, ...
    param,...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);


% 5.1.B If running on RR control analysis, need to run only one inputdatatype per
% TD (since RR 200 ms is for TD 1 and RR 400 ms is for TD 2), so use this:
% SignElecIndices = NASTD_ECoG_Predict_PredEffectTable_RR_lk...
%     (vars.validSubjs, sub_list, InputDataType, InputEffectType, ToneDur_text, ...
%     param,...
%     plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
% 
% 
% % 5.1.C If running on FTPL control analysis
% SignElecIndices = NASTD_ECoG_Predict_PredEffectTable_FTPL_lk...
%     (vars.validSubjs(2:9), sub_list(2:9), InputDataType, InputEffectType, ToneDur_text, ...
%     param,...
%     plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);

%% 6) Plot electrodes showing sign. prediction effect
param.pval_plotting     = 0.025; %Pval thresh for plotting
param.FDRcorrect        = 1;
param.pval_FDR          = 0.05;
param.plot_SubplotperTW = 1; %Common plot across TW
param.ElecSelect        = 'All'; %StimCorr, All
%  set(0, 'DefaultFigureVisible', 'off') %show surface plots

param.SamplesTW         = 25; %25 = 50ms TW
param.Label_TW          = num2str(round(param.SamplesTW/512,2));
param.ToneIndex         = 33;

%InputDataType           = {'Alpha_LogAmp', 'Beta_LogAmp'}; %, 'HighGamma_LogAmp'};
%InputDataType           = {'RR_200ms_TD_LogAmp', 'RR_400ms_TD_LogAmp'};
InputDataType = {'HP05toLP30Hz','HighGamma_LogAmp'};
%InputDataType = {'Gamma_LogAmp'};
%     {'LP35Hz', ...
%     'HP05toLP30Hz', 'HP1toLP30Hz',  'HP2toLP30Hz', ...
%     'Alpha_LogAmp', 'Beta_LogAmp', ...
%     'Gamma_LogAmp', 'HighGamma_LogAmp'};
InputEffectType         = {'PredEffect', 'ComplexPredErrEffect'}; %{'PredEffect', 'SimplePredErrEffect', 'ComplexPredErrEffect'};
InputEffectType = {'PredEffect'};
%InputEffectType = {'PredErrFTPL'}; %This is for the FTPL control analysis
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
subs = sub_list(vars.validSubjs);

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

% 6.2.B All subjects on one volume for time-windowed results for RR control
% analysis
% for i_pval = 1:length(param.pval_plotting) %p-value for plotting
%     for i_effect = 1:length(InputEffectType)
%         for i_inputData = 1:length(InputDataType)
%                 NASTD_ECoG_Predict_PlotSignElec_AllSub...
%                     (subs, ...
%                     InputDataType{i_inputData}, InputEffectType{i_effect}, ToneDur_text{i_inputData}, ...
%                     param,...
%                     save_poststepFigs, paths_NASTD_ECoG)
%         end
%     end
% end

%6.3.1 All subjects on one volume for sample-wise results per TD
subs = sub_list(vars.validSubjs);
% subs = sub_list(2:9);

% for i_effect = 1:length(InputEffectType)
%     for i_inputData = 1:length(InputDataType)
%         for i_TD = 1:length(ToneDur_text) %tone duration condition
%             
%             %%Plot specific parameter (i.e., clusterstat, p-val, onset, offset, effect duration)
%             %             for i_param2plot = 1:length(param2plot) %1 plot per parameter
%             %                 param.param2plot = param2plot{i_param2plot};
%             %                 NASTD_ECoG_Predict_PlotSignClusterElec_AllSub...
%             %                     (subs, ...
%             %                     InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%             %                     param,...
%             %                     save_poststepFigs, paths_NASTD_ECoG)
%             %             end
%             
%             %Plot all parameters on a summary plot
% %             NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSub...
% %                 (subs, ...
% %                 InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
% %                 param,...
% %                 save_poststepFigs, paths_NASTD_ECoG)
%             
%              %Plot all parameters on a summary plot (no cluster corr)
%             NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSub_lk...
%                 (subs, ...
%                 InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                 param,...
%                 save_poststepFigs, paths_NASTD_ECoG)
%             
%             
%             %Plot summary figure showing 
%             %1) traces per p*34 condition during p33 for all sign. elecs
%             %2) traces per p*34 condition during p1-p33 for all sign. elecs
%             %3) surface plot overview
%             %for all subjects, with traces grouped per anat region
%             if strcmp(InputEffectType{i_effect}, 'PredEffect')
%                 %Note: Makes only sense for pred effect, since p*34 distinction
%                 %is not useful for PE
%                 NASTD_ECoG_Predict_PlotERF4SignClusterElec_AllSub...
%                     (subs, ...
%                     InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
%                     param,...
%                     save_poststepFigs, paths_NASTD_ECoG)
%             end            
%             
%         end
%     end
% end

% 6.3.1.B For RR control analysis
% for i_effect = 1:length(InputEffectType)
%     for i_inputData = 1:length(InputDataType)
%         NASTD_ECoG_Predict_PlotSignClusterElec_AllSub_FDR_Uncorr...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text{i_inputData}, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)
%     end
% end
            
%6.3.2 All subjects on one volume for sample-wise results aggregated across TD
subs = sub_list(vars.validSubjs);
% subs = sub_list(2:9);

Onset_times_all = {}; %3 conds of interest: pred, simple PE and complex PE, each one will contain an array with format electrodes w/ sig cluster x
save_poststepFigs = 1;

for i_effect = 1:length(InputEffectType)
    for i_inputData = 1:length(InputDataType)
        
        % Approach for computing the onset time as the first significant
        % cluster
%         Onset_times_ms = NASTD_ECoG_Predict_PlotSignClusterElec_SubplotParam_AllSubTD_lk...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)
%         
%         Onset_times_all{i_effect} = Onset_times_ms;
        
        % Plot cluster-corrected in small spheres and FDR-corrected in
        % large yellow spheres
        NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_FDR_Uncorr...
            (subs, ...
            InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
            param,...
            save_poststepFigs, paths_NASTD_ECoG)
        
        % Plot on both HEMIs
%         NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_FDR_Uncorr_Both...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)

        % Plot cluster-corrected in small spheres and FDR-corrected in
        % large yellow spheres COLOR CODE BY SUBJECT ID
%         NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_BySubj_FDR...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)
        
%         
        % Approach for computing the onset time as the first FDR-corrected
        % significant p-value
%         Onset_times_ms = NASTD_ECoG_Predict_PlotSignClusterElec_OnsetTimes...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG);
%         
%         Onset_times_all{i_effect} = Onset_times_ms;
%         
%         %figure output saved as movie
%         NASTD_ECoG_Predict_PlotFDRSignClusterElec_AllSubTD_Mov...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG)

        % For FTPL control analysis only
%         NASTD_ECoG_Predict_PlotSignClusterElec_AllSubTD_FDR_Uncorr_FTPL...
%             (subs, ...
%             InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
%             param,...
%             save_poststepFigs, paths_NASTD_ECoG);
        
        %Plot summary figure showing
        %1) traces per p*34 condition during p33 for all sign. elecs
        %2) traces per p*34 condition during p1-p33 for all sign. elecs
        %3) surface plot overview
        %for all subjects, with traces grouped per anat region
        %% March 2026 - plotting p33 and p34 4per p*34 and p32 and evoked p32
        if strcmp(InputEffectType{i_effect}, 'PredEffect')
            %Note: Makes only sense for pred effect, since p*34 distinction
            %is not useful for PE
            NASTD_ECoG_Predict_PlotERF4SignClusterElec_AllSubTD_lk...
                (subs, ...
                InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
                param,...
                save_poststepFigs, paths_NASTD_ECoG)
        end
        
        %% To extract Long p34 and plot
        if strcmp(InputEffectType{i_effect}, 'PredEffect')
            %Note: Makes only sense for pred effect, since p*34 distinction
            %is not useful for PE
            NASTD_ECoG_Predict_PlotERF4SignClusterElec_longp34_lk...
                (subs, ...
                InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
                param,...
                save_poststepFigs, paths_NASTD_ECoG)
        end
        
         %% To run regression of neural activity onto tone pitch across sequence and for p1 pitch tracking
        if strcmp(InputEffectType{i_effect}, 'PredEffect')
            %Note: Makes only sense for pred effect, since p*34 distinction
            %is not useful for PE
            NASTD_ECoG_Predict_PlotERF4SignClusterElec_pitchregression_lk...
                (subs, ...
                InputEffectType{i_effect}, InputDataType{i_inputData}, ToneDur_text, ...
                param,...
                save_poststepFigs, paths_NASTD_ECoG)
        end
    end
end

%% Bar plots 
% Create figure
% data1_col1 = Onset_times_all{1}(:,1);  % Prediction TD1
% data1_col2 = Onset_times_all{1}(:,2);  % Prediction TD2
% data3_col1 = Onset_times_all{2}(:,1);  % Prediction Error TD1 (complex)
% data3_col2 = Onset_times_all{2}(:,2);  % Prediction Error TD2 (complex)
% means = [nanmean(data1_col1), nanmean(data1_col2), nanmean(data3_col1), nanmean(data3_col2)];
% sem = [nanstd(data1_col1)/length(data1_col1(~isnan(data1_col1))), nanstd(data1_col2)/length(data1_col2(~isnan(data1_col2))), nanstd(data3_col1)/length(data3_col1(~isnan(data3_col1))), nanstd(data3_col2)/length(data3_col2(~isnan(data3_col2)))];
% T = array2table([means; sem], 'VariableNames', {'Pred_TD1', 'Pred_TD2', 'PE_TD1', 'PE_TD2'}, 'RowNames', {'Means', 'SEM'});
% figures_dir = '/isilon/LFMI/VMdrive/Lua/NASTD/Figures/';
% writetable(T, [figures_dir, 'OnsetTimePerTD_', InputDataType{i_inputData}, '_FDRmethod.csv'], 'WriteRowNames', true);
% 
% d1 = data1_col1(~isnan(data1_col1));
% d2 = data1_col2(~isnan(data1_col2));
% d3 = data3_col1(~isnan(data3_col1));
% d4 = data3_col2(~isnan(data3_col2));
% maxLength = max([length(d1), length(d2), length(d3), length(d4)]);
% 
% % Pad shorter vectors with NaN
% d1_padded = [d1; NaN(maxLength - length(d1), 1)];
% d2_padded = [d2; NaN(maxLength - length(d2), 1)];
% d3_padded = [d3; NaN(maxLength - length(d3), 1)];
% d4_padded = [d4; NaN(maxLength - length(d4), 1)];
% 
% % Create the table
% T2 = table(d1_padded, d2_padded, d3_padded, d4_padded, 'VariableNames', {'d1', 'd2', 'd3', 'd4'});
% writetable(T2, [figures_dir, 'OnsetTimePerTD_', InputDataType{i_inputData}, 'alldat_FDRmethod.csv'], 'WriteRowNames', true);
% 
% T2_imported = readtable([figures_dir, 'OnsetTimePerTD_', InputDataType{i_inputData}, 'alldat_FDRmethod.csv'], 'ReadRowNames', false);
% 
% % Run ANOVA to see if theres an effect
% data = [data1_col1; data1_col2; data3_col1; data3_col2];
% timepoint = [repmat('TD1', length(data1_col1), 1); repmat('TD2', length(data1_col2), 1); repmat('TD1', length(data3_col1), 1); repmat('TD2', length(data3_col2), 1)];
% condition = [repmat('Pred', length(data1_col1), 1); repmat('Pred', length(data1_col2), 1); repmat('PE  ', length(data3_col1), 1); repmat('PE  ', length(data3_col2), 1)];
% tbl = table(data, timepoint, condition, 'VariableNames', {'Data', 'Timepoint', 'Condition'});
% tbl = tbl(~isnan(tbl.Data), :);
% lm = fitlm(tbl, 'Data ~ Timepoint*Condition');
% anova_results = anova(lm);
% disp(anova_results);
% 
% %Wilcoxon
% [p1,h1,stats1] = ranksum(data1_col1(~isnan(data1_col1)),data3_col1(~isnan(data3_col1)));
% [p2,h2,stats2] = ranksum(data1_col2(~isnan(data1_col2)),data3_col2(~isnan(data3_col2)));
% 
% T2_imported = table2array(T2_imported);
% means = [nanmean(T2_imported(:,1)), nanmean(T2_imported(:,2)), nanmean(T2_imported(:,3)), nanmean(T2_imported(:,4))];
% sem = [nanstd(T2_imported(:,1))/sum(~isnan(T2_imported(:,1))), nanstd(T2_imported(:,2))/sum(~isnan(T2_imported(:,2))), nanstd(T2_imported(:,3))/sum(~isnan(T2_imported(:,3))), nanstd(T2_imported(:,4))/sum(~isnan(T2_imported(:,4)))];
% 
% % Group bar plot
% figure;
% colors = {[0.5, 0.7, 1], [1, 0.6, 0], [0.5, 0.7, 1], [1, 0.6, 0]};
% offsets = {0, 0.5, 0, 0.5};
% hold on;  
% b1 = bar(1, means(1), 'FaceColor', colors{1}, 'BarWidth', 0.4);
% b2 = bar(2, means(2), 'FaceColor', colors{2}, 'BarWidth', 0.4);
% b3 = bar(3, means(3), 'FaceColor', colors{3}, 'BarWidth', 0.4);
% b4 = bar(4, means(4), 'FaceColor', colors{4}, 'BarWidth', 0.4);
% set(b1, 'XData', 1);
% set(b2, 'XData', 1.5);  % Move bar 2 closer to bar 1
% set(b3, 'XData', 3);
% set(b4, 'XData', 3.5);  % Move bar 4 closer to bar 3
% % errorbar(1, means(1), sem(1), 'k.', 'LineWidth', 2);
% % errorbar(1.5, means(2), sem(2), 'k.', 'LineWidth', 2);
% % errorbar(3, means(3), sem(3), 'k.', 'LineWidth', 2);
% % errorbar(3.5, means(4), sem(4), 'k.', 'LineWidth', 2);
% jitterAmount = 0.3; 
% for i = 1:4
%     groupdata = T2_imported(:,i);
%     lengthGroup = sum(~isnan(T2_imported(:,i)));
%     xPositions = repelem(i, lengthGroup) - offsets{i} + ...
%              (rand(size(groupdata)) - 0.5) * jitterAmount;
%     scatter(repelem(i, length(groupdata)) - offsets{i}, groupdata, 'filled', ...
%         'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'k', 'XJitter','randn', 'XJitterWidth', jitterAmount);
% end
% xticks([]);
% xticklabels({'', '', '', ''});  % We will manually set labels below the bars
% text(1.25, -0.8, 'Pred', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 18);
% text(3.25, -0.8, 'PE', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 18);
% legend([b1, b2], {'200 ms', '400 ms'}, 'Location', 'northwest', 'FontSize', 16);
% ylabel('Mean Onset Time (ms)', 'FontSize', 18);
% ax = gca; % Get the current axes handle
% ax.FontSize = 16; % Set the font size of the axis tick labels
% ylim([0 300]);
% hold off;


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
