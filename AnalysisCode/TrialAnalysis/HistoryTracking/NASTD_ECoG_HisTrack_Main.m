%% Project: NASTD_ECoG
%% Task: Compute History Tracking effects per subject/electrode/tone duration condition

%1. Load in pre-processed data
%2. Prepare input data (ERF-range and amplitude envelopes)
%3a. Compute experimental Kprime values (includes across-fold combination)
%3b. Compute shuffled Kprime values
%3c. Combine experimental and shuffled Kprime computations to 1 Kprime output
%4. Plot Kprime values (thresholded vs. un-thresholded, single-subject vs
%group-level)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%Determine subjects
sub_list = vars.sub_list;

%Load in file with individual preproc infos
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

%Set analysis parameters
ToneDur_label = {'0.2' '0.4'};
inputData_label = { 'HP05toLP30Hz'};
% inputData_label = {'HP05toLP30Hz','Alpha_LogAmp','Beta_LogAmp','Gamma_LogAmp','HighGamma_LogAmp'};
    
fsample     = 512;
SamplesTW   = 25; %25 = 50ms TW
Label_TW    = [num2str(round(SamplesTW/fsample,2)*1000) 'msTW'];
NumFolds    = 5;  %Number folds used for Kprime computation (1 test vs. NumFolds - 1 train folds)
NumRuns     = 5;  %Number runs used for Kprime computation. Estimates averaged across runs
nReps       = 1000; %Number repetitions for shuffled Kprime distribution

%Set plotting parameters
plot_poststepFigs = 1;
save_poststepFigs = 0;

%% 1. Load in pre-processed data
for i_sub = vars.validSubjs
    
    sub = sub_list{i_sub};
    tic
    disp(['Loading and processing subject: ' sub])    
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];%path to indiv preprocessed ECoG data
    
    % 0.3) Load in preproc data
    load(loadfile_ECoGpreprocdata);

%     SampleFreq = DataClean_AllTrials.fsample;
    
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
    for i_TD = 1:length(ToneDur_label)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_Preproc_SelTrialsperTD...
            (ToneDur_label{i_TD}, DataClean_CleanTrialsElecs);
        
        disp(['Total trial number for TD ' ToneDur_label{i_TD} 's: ' ...
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
    
    %2.4A Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus baseline
    BL_win = [-0.5 0];    
    for i_TD = 1:length(ToneDur_label)
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
    
    %2.4B Remove trials where some electrodes have only NaN samples (NY798)
    for i_TD = 1:length(ToneDur_label)
        DataClean_CleanTrialsElecs_perTD{i_TD} = ...
            NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
            (sub, ToneDur_label{i_TD}, DataClean_CleanTrialsElecs_perTD{i_TD});   
    end    
        
    %2.5 Apply Filter for Slow ERF-range activity and deal with NaNs
%     for i_TD = 1:length(ToneDur_label)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.LP35Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.LP35Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_LP35Hz...
%             (sub, ToneDur_label{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             0, 0, paths_NASTD_ECoG);
%     end    
%     for i_TD = 1:length(ToneDur_label)
%         [DataClean_CleanTrialsElecs_perTD{i_TD}.HP2toLP30Hz, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP2toLP30Hz] = ...
%             NASTD_ECoG_FiltNaNinterp_HP2toLP30Hz...
%             (sub, ToneDur_label{i_TD}, ...
%             DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%             plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
%     end  
    for i_TD = 1:length(ToneDur_label)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.HP05toLP30Hz, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, ToneDur_label{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    end
    
%     %2.6 Compute Frequency-band specific amplitude envelope and deal with NaNs
%     FrequencyBands = [70, 150]; %Hz [HP,LP]
%     FrequencyBand_Labels = {'HighGamma'};
% %     FrequencyBands = [8 12; 15 30; 30, 70; 70, 150]; %Hz [HP,LP]
% %     FrequencyBand_Labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};
%     
%     for i_TD = 1:length(ToneDur_label)
%         for i_freqbands = 1:length(FrequencyBand_Labels)
%             
%             FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
%             
%             [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
%                 NASTD_ECoG_FiltNaNinterp_AmpEnvel...
%                 (sub, ToneDur_label{i_TD}, ...
%                 FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel,...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 0, 0, paths_NASTD_ECoG);
%             
%         end
%     end
    
    disp(['Done loading and preprocessing subject: ' sub ' in ' ...
        num2str(toc/60) ' min'])
    
    
%% 3) Compute Kprime-values 
    %3.1 Experimental K' (~1 min per run)    
    plotFig_CompFolds = 1; %plot summary fig showing avg Kprime and min squared residuals across folds 
    plot_SubplotperTW = 1; %Subplot per TW
    param.ElecSelect  = 'All'; %StimCorr, All

    for i_inputData = 1:length(inputData_label)
        for i_TD = 1:length(ToneDur_label)

            tic
            process_info = (['-- Starting computation of EXP Kprime for Sub: ' ...
                sub ' - Input: ' inputData_label{i_inputData} ...
                ' - TD: ' ToneDur_label{i_TD} ' --']);
            disp(process_info)

%             %Compute exp Kprime in 1 run (CAVE: Kprime estimate not overly
%             %reliable)
%             NASTD_ECoG_HisTrack_CompExpKvals...
%                 (sub, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 SamplesTW, Label_TW, NumFolds, ...
%                 paths_NASTD_ECoG, ...
%                 plotFig_CompFolds);
            
            %Compute exp Kprime in multiple runs and average estimates across runs
            NASTD_ECoG_HisTrack_CompExpKvals_MultiRun...
                (sub, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                SamplesTW, Label_TW, NumFolds, NumRuns, ...
                paths_NASTD_ECoG, ...
                plotFig_CompFolds);            

            process_info = (['-- Finished computation of EXP Kprime for Sub: ' ...
                sub ' - Input: ' inputData_label{i_inputData} ...
                ' - TD: ' ToneDur_label{i_TD} ...
                ' after (' num2str(toc/60) ' min)--']);
            disp(process_info)
                        
            %Surface-plot EXP Kprime results
            NASTD_ECoG_HisTrack_PlotExpKvals_SsubperTW...
                (sub, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
                SamplesTW, Label_TW, NumRuns,...
                param, ...
                plot_SubplotperTW, save_poststepFigs, ...
                paths_NASTD_ECoG)            

        end
    end
    
%     %3.2 Shuffled K' - Computed on BigPurple Cluster
%     for i_inputData = 1:length(inputData_label)
%         for i_TD = 1:length(ToneDur_label)
%            
%             tic
%             process_info = (['-- Starting computation of SHUFF Kprime for Sub: ' ...
%                 sub ' - Input: ' inputData_label{i_inputData} ...
%                 ' - TD: ' ToneDur_label{i_TD} ' --']);
%             disp(process_info)           
%        
%             NASTD_ECoG_HisTrack_CompShuffKvals...
%                 (sub, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
%                 DataClean_CleanTrialsElecs_perTD{i_TD}, ...
%                 SamplesTW, Label_TW, NumFolds, nReps,...                
%                 paths_NASTD_ECoG, ...
%                 plotFig_CompFolds);
%             
%             process_info = (['-- Finished computation of SHUFF Kprime for Sub: ' ...
%                 sub ' - Input: ' inputData_label{i_inputData} ...
%                 ' - TD: ' ToneDur_label{i_TD} ...
%                 ' after (' num2str(toc/60) ' min)--']);
%             disp(process_info)            
%   
%         end
%     end

clear DataClean_CleanTrialsElecs_perTD
end

%2.3 Combine Exp and Shuff Kprime (combined across folds) to one outputfile
for i_sub = vars.validSubjs
    for i_inputData = 1:length(inputData_label)
        for i_TD = 1:length(ToneDur_label)
          
            tic
            process_info = (['-- Fusing Exp and Shuffled combined K for Sub: ' ...
                sub_list{i_sub} ' Input: ' inputData_label{i_inputData} ...
                ' TD: ' ToneDur_label{i_TD} ' --']);
            disp(process_info)
       
            NASTD_ECoG_HisTrack_CombineExpShuffK...
                (sub_list{i_sub}, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
                SamplesTW, Label_TW, NumRuns, ...
                plot_poststepFigs, ...
                paths_NASTD_ECoG)
            
            process_info = (['-- Finished fusing Exp and Shuffled combined K (' ...
                num2str(toc) ' sec) for Sub: ' sub_list{i_sub} ...
                ' Input: ' inputData_label{i_inputData} ...
                ' TD: ' ToneDur_label{i_TD} ' --']);
            disp(process_info)             
                     
        end
    end
end


%% 3) Plot exp Kprime values 
save_poststepFigs = 1;
plot_SubplotperTW = 0; %Subplot per TW
param.ElecSelect  = 'All'; %StimCorr, All
if strcmp(param.ElecSelect, 'All')
    vars.validSubjs = 1:9;
elseif strcmp(param.ElecSelect, 'StimCorr')
    vars.validSubjs = 2:9;
end
   
%3.1 Plot unthresholded (all) Exp Kprime values 
%3.1.1 For each subject (per time window)
for i_sub = vars.validSubjs
    for i_inputData = 1:length(inputData_label)
        for i_TD = 1:length(ToneDur_label)       
            
            NASTD_ECoG_HisTrack_PlotExpKvals_SsubperTW...
                (sub_list{i_sub}, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
                SamplesTW, Label_TW, NumRuns, ...
                param, ...
                plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)
            
        end
    end
end

%3.1.2 All subjects together (per time window)
for i_inputData = 1:length(inputData_label)
    for i_TD = 1:length(ToneDur_label)
        
        NASTD_ECoG_HisTrack_PlotExpKvals_AllSubperTW...
            (sub_list(vars.validSubjs), ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
            SamplesTW, Label_TW, NumRuns, ...
            param, ...
            plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)
        
    end
end

%3.2 Plot thresholded (Exp vs Shuff p < 0.05) Exp Kprime values with
%optional FDR correction (takes p-value from thresholding)
pval_threshold  = 0.05;
FDR_correct     = 1;

%3.2.1 For each subject (per and averaged across TW)
for i_sub = vars.validSubjs
    for i_inputData = 1:length(inputData_label)
        for i_TD = 1:length(ToneDur_label)       
                
                NASTD_ECoG_HisTrack_PlotSignExpKvals_SsubperTW...
                    (sub_list{i_sub}, ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
                    SamplesTW, Label_TW,...
                    param, ...
                    FDR_correct, pval_threshold,...
                    plot_SubplotperTW, save_poststepFigs, ...
                    paths_NASTD_ECoG)
                
        end
    end
end

%3.2.2 All subjects together (per and averaged across TW)
for i_inputData = 1:length(inputData_label)
    for i_TD = 1:length(ToneDur_label)
        
        %plot figure
        NASTD_ECoG_HisTrack_PlotSignExpKvals_AllSubperTW...
            (sub_list(vars.validSubjs), ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
            SamplesTW, Label_TW,...
            param, ...
            FDR_correct, pval_threshold,...
            plot_SubplotperTW, save_poststepFigs, ...
            paths_NASTD_ECoG)      
        
%         %save movie (1 frame per TW)
%         NASTD_ECoG_HisTrack_PlotSignExpKvals_AllSubperTW_Mov...
%             (sub_list(vars.validSubjs), ToneDur_label{i_TD}, inputData_label{i_inputData}, ...
%             SamplesTW, Label_TW,...
%             param, ...
%             FDR_correct, pval_threshold,...
%             paths_NASTD_ECoG)

    end
end

%% 4. Compute and plot Kprime ratio across TD
save_poststepFigs = 1;
plot_SubplotperTW = 1; %Subplot per TW
param.ElecSelect  = 'All'; %StimCorr, All
if strcmp(param.ElecSelect, 'All')
    vars.validSubjs = 1:9;
elseif strcmp(param.ElecSelect, 'StimCorr')
    vars.validSubjs = 2:9;
end

%4.1 Plot unthresholded (all electrodes) Kprime ratios across TDs
%4.1.1 For each subject (per time window)
for i_sub = vars.validSubjs 
    for i_inputData = 1:length(inputData_label)
            
            NASTD_ECoG_HisTrack_PlotExpKvalsBothTD_SsubperTW...
                (sub_list{i_sub}, ToneDur_label, inputData_label{i_inputData}, ...
                SamplesTW, Label_TW, NumRuns, ...
                param, ...
                plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)
            
    end
end

%4.1.2 All subjects together (per time window)
for i_inputData = 1:length(inputData_label)%         
        NASTD_ECoG_HisTrack_PlotExpKvalsBothTD_AllSubperTW...
            (sub_list(vars.validSubjs), ToneDur_label, inputData_label{i_inputData}, ...
            SamplesTW, Label_TW, NumRuns, ...
            param, ...
            plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)        
end

%4.2 Plot thresholded (only electrodes where Exp vs Shuff p < 0.05) 
    %Exp Kprime values for both TD
pval_threshold  = 0.05;
FDR_correct     = 0;

%4.2.1 For each subject (per time window)
for i_sub = vars.validSubjs 
    for i_inputData = 1:length(inputData_label)
                
        NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_SsubperTW...
            (sub_list{i_sub}, ToneDur_label, inputData_label{i_inputData}, ...
            SamplesTW, Label_TW, NumRuns, ...
            FDR_correct, pval_threshold,...
            param, ...
            plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)
                
    end
end

%4.2.2 For all subjects (per time window)
for i_inputData = 1:length(inputData_label)
    
    NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AllSubperTW...
        (sub_list(vars.validSubjs), ToneDur_label, inputData_label{i_inputData}, ...
        SamplesTW, Label_TW, NumRuns, ...
        FDR_correct, pval_threshold,...
        param, ...
        plot_SubplotperTW, save_poststepFigs, paths_NASTD_ECoG)
    
end

%4.2.3 For all subjects avg across time windows
for i_inputData = 1:length(inputData_label)
    
    NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AllSubAvgTW...
        (sub_list(vars.validSubjs), ToneDur_label, inputData_label{i_inputData}, ...
        SamplesTW, Label_TW, ...
        FDR_correct, pval_threshold,...
        param, ...
        save_poststepFigs, paths_NASTD_ECoG)
    
end


%% 5) Create anatomical filter to compute and compare Kprime values for anatomical parcellations
param.ElecSelect  = 'All'; %StimCorr, All
if strcmp(param.ElecSelect, 'All')
    vars.validSubjs = 1:9;
elseif strcmp(param.ElecSelect, 'StimCorr')
    vars.validSubjs = 2:9;
end

pval_threshold      = 0.05;
FDR_correct         = 1;

%5.1 For all subjects (per and across time windows)
for i_inputData = 1:length(inputData_label)
    %For each anatomical region: avg kprime per subject, TD, kprime ratio
    
    NASTD_ECoG_HisTrack_PlotSignExpKvalsBothTD_AnatReg_AllSubperTW...
        (sub_list(vars.validSubjs), ToneDur_label, inputData_label{i_inputData}, ...
        SamplesTW, Label_TW, ...
        FDR_correct, pval_threshold,...
        param, ...
        save_poststepFigs, paths_NASTD_ECoG)
    
end
