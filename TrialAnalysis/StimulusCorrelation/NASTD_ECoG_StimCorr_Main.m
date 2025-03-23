%% Project: NASTD_ECoG
%Compute Stimulus correlation analyses
%Determine electrodes that show tone sequence tracking to indentical vs.
%similar ton sequences.

%1. Load in pre-processed data
%2. Prepare input data (ERF-range and amplitude envelopes)
%3. Compute Stimulus correlation analyses

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

for i_sub = vars.validSubjs
    
    tic
    disp(['-- Preparing input data for sub: ' sub_list{i_sub} ' --'])
   
%% 1. Set paths and load pre-processed data   
    sub = sub_list{i_sub};
    
   stimcorr_dir = [paths_NASTD_ECoG.ECoGdata_StimCorr sub '/'];
    if (~exist(stimcorr_dir, 'dir')); mkdir(stimcorr_dir); end
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];    
    load(loadfile_ECoGpreprocdata);
    
    SampleFreq = DataClean_AllTrials.fsample;
    
       
%% 2. Prepare input data
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
    

    %% 3. Prepare input signals: slow arrhythmic activity & gamma amplitude envelope
    
    %3.0 Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus baseline
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
    
    %3.1 Remove trials where some electrodes have only NaN samples (NY798)
    for i_TD = 1:length(ToneDur_text)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
            (sub, ToneDur_text{i_TD}, DataClean_CleanTrialsElecs_perTD{i_TD});   
    end    
        
%Old version
%     %3.2 Apply Filter for Slow ERF-range activity and deal with NaNs
%     for i_TD = 1:length(ToneDur_text)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.LP35Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.LP35Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_LP35Hz...
%             (sub, ToneDur_text{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end 
    
%New Version
    %3.2 Apply Filter for Slow ERF-range activity and deal with NaNs
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.HP05toLP30Hz, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    end    
    
    
%     %3.3 Compute Frequency-band specific amplitude envelope and deal with NaNs
%     FrequencyBands = [8 12; 15 30; 30 70; 70 150]; %Hz [HP,LP]
%     FrequencyBand_Labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
%     
%     for i_TD = 1:length(ToneDur_text)
%         for i_freqbands = 1:length(FrequencyBand_Labels)            
%             FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];            
%             [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
%                 NASTD_ECoG_FiltNaNinterp_AmpEnvel...
%                 (sub, ToneDur_text{i_TD}, ...
%                 FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);            
%         end
%     end

    disp(['-- Finished input data for sub: ' sub_list{i_sub} ' after ' ...
        num2str(round(toc/60,2)) 'min --'])
    
        
%% 3. Comute stimulus correlation analyses
    %3.0 Parameters:
    InputDataType = {'HP05toLP30Hz'};
%     InputDataType = {'HP05toLP30Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
    
    %3.1 Compute sequence processing / Stimulus correlation effect (identical and similar)
    %for specific input data and TD    
    param.StimCorr_Seqrange = 33; 
    %final tone defining range for tone pitches based on which 
    %stimulus correlation should be computed (i.e., 1-33)
    param.TW_Samples = 25; %25 = 50ms TW
    param.TW_Label = num2str(round(param.TW_Samples/SampleFreq,2));

    for i_inputData = 1:length(InputDataType)
        for i_TD = 1:length(ToneDur_text)
            NASTD_ECoG_StimCorr_IdentvsSim ...
                (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                save_poststepFigs, paths_NASTD_ECoG);
        end
    end

    %3.2 Compute tone processing effect for specific input data and TD    
    param.N100_TW       = [0.075 0.125]; %should focus on N100-peak around 100 ms
    param.BL_tone       = [0 0.05]; %baseline at tone beginning
    param.pval_plotting = 0.05; %Pval thresh to select sign. electrode for plotting

    for i_inputData = 1:length(InputDataType)
        for i_TD = 1:length(ToneDur_text)
            NASTD_ECoG_StimCorr_ToneProc ...
                (sub, InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
                param,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                save_poststepFigs, paths_NASTD_ECoG);
        end
    end
    
    clear DataClean_CleanTrialsElecs_perTD
    
end

%% 4. Plot aggregated group-data on common surface
param.FDRcorrect    = 0;
param.pval_plotting = 0.05; %Pval thresh for plotting

InputDataType = {'HP05toLP30Hz', 'HighGamma_LogAmp'};
% InputDataType   = {'HP05toLP30Hz','Alpha_LogAmp', 'Beta_LogAmp', 'Gamma_LogAmp', 'HighGamma_LogAmp'};
ToneDur_text    = {'0.2' '0.4'};

%4.1 Compare tone proc vs similar vs identical electrodes
for i_inputData = 1:length(InputDataType)
    for i_TD = 1:length(ToneDur_text)
        
        NASTD_ECoG_StimCorr_CompareResElecs...
            (subs, ...
            InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
            param,...
            save_poststepFigs, paths_NASTD_ECoG)
        
    end
end    

%4.2 All subjects on one volume separately per TD
% subs = sub_list(vars.validSubjs);
subs = sub_list(2:9);

for i_inputData = 1:length(InputDataType)
    for i_TD = 1:length(ToneDur_text)
        
        NASTD_ECoG_StimCorr_PlotSignElec_AllSub...
            (subs, ...
            InputDataType{i_inputData}, ToneDur_text{i_TD}, ...
            param,...
            save_poststepFigs, paths_NASTD_ECoG)
        
    end
end

%4.3 All subjects on one volume aggregated across TD
for i_inputData = 1:length(InputDataType)
    
    NASTD_ECoG_StimCorr_PlotSignElec_AllSubTD...
        (subs, ...
        InputDataType{i_inputData}, ToneDur_text, ...
        param,...
        save_poststepFigs, paths_NASTD_ECoG)
        
end
