%% Project: NASTD_ECoG
%Compute Time-Frequency-Analysis (TFA) across all frequency bands and the entire sequence

%1. Load in pre-processed data
%2. Prepare input data (baseline correction)
%3. Compute TFA across trials per electrode, TD, and subject
%4. Plot output 

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add project base path

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

ToneDur_text = {'0.2' '0.4'};

plot_poststepFigs = 1;
save_poststepFigs = 1;

%% 1. Load in pre-processed data
for i_sub = vars.validSubjs
    
    sub = sub_list{i_sub};
    disp(['Processing subject: ' sub])
    
    prediction_dir = [paths_NASTD_ECoG.ECoGdata_Prediction sub '/'];
    if (~exist(prediction_dir, 'dir')); mkdir(prediction_dir); end
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    tic
    disp(['Loading preprocessed data set for sub: ' sub])
    load(loadfile_ECoGpreprocdata);
    disp(['done loading in ' num2str(toc) ' sec'])
    
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
    end
    
    %Cleanup
    clear DataClean_CleanTrialsElecs DataClean_CleanTrials DataClean_AllTrials
    
    %No BL correction in time domain, later done in frequency domain
    
    %% 3. Compute TFA
    %Decompose the signal in time and frequency using time-resolved 
    %Fourier-based spectral decomposition.

    param.TFA_FOI = {[6:1:35], [30:5:150]}; %, [5:5:150]
    
    %3.1 Compute TFA across all trials
    %To do: Better plotting for topoplot electrode size
    for i_freq = 1:length(param.TFA_FOI)
        for i_TD = 1:length(ToneDur_text)        
            NASTD_ECoG_TFA_TFAallTrials...
                (sub, ToneDur_text{i_TD}, ...
                param.TFA_FOI{i_freq}, ...  
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);        
        end
    end
    
%     %3.2 Linearly correlate prediction performance metric with power values
%     %and check if power changes are indicative of prediction performance   
%     for i_freq = 1:length(param.TFA_FOI)
%         for i_TD = 1:length(ToneDur_text)        
%             NASTD_ECoG_TFA_CorrPredPerfPowChanges...
%                 (sub, ToneDur_text{i_TD}, ...
%                 param.TFA_FOI{i_freq}, ...  
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);        
%         end 
%     end

    %3.3 Outdated: Stat comparison between high vs. low performance trials
    %     for i_freq = 1:length(param.TFA_FOI)
    %         for i_TD = 1:length(ToneDur_text)
    %             NASTD_ECoG_TFA_PowComparisonFTPL...
    %                 (sub, ToneDur_text{i_TD}, ...
    %                 param.TFA_FOI{i_freq}, ...
    %                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
    %                 plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    %         end
    %     end
    
end