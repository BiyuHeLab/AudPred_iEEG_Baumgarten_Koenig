function GCdata = ...
    NASTD_ECoG_Connectivity_CalculateGC_PEeffect ...
    (sub_list, i_sub, ...
    ToneDur_text, SelElecs, ...
    param, paths_NASTD_ECoG)

%Aim: Compute temporal and spectral GC estimates for simple and complex prediction errors
%during p34 as a function of p*34-p34 difference
%for all selected electrode pairs in a specific subject.
%Return temporal ans spectral Granger values and p-values.


%% 0.1) Set empty GC data struct and output paths
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/core/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/utils/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/stats/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/gc/');

param.GC.tone_aggregation = {34}; %PE computation only makes sense for p34

%% 1) Load and prepare single-subject data

%Load in preproc single-trial data
sub = sub_list{i_sub};
disp(['Computing CG for subject: ' sub]);

NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
tic
disp(['Loading preprocessed data set for sub: ' sub])
load(loadfile_ECoGpreprocdata, 'DataClean_AllTrials');
disp(['done loading in ' num2str(toc) ' sec'])

%Save basic GC info
GCdata                      = struct;
GCdata.info.sub             = sub;
GCdata.info.InputDataType   = param.GC.InputDataType;
GCdata.info.ElecPairEffect  = param.GC.ElecPairEffect;
GCdata.info.ElecPairRegion  = param.GC.ElecPairRegion;
GCdata.info.SelTones        = param.GC.tone_aggregation;
GCdata.info.epochdur_ms     = param.GC.epochdur_ms;
GCdata.info.param           = param.GC;

%1.1) Select only clean/valid trials
DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);

%1.2) Select only valid electrodes (Grid+Strip, clean, MNI-coordinates present)
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

%1.3) Select only trials with specific Tone Duration (TD)
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

%1.4 Determine TP/samples for each tone start+end
for i_TD = 1:length(ToneDur_text)
    curr_toneDur_text = ToneDur_text{i_TD};
    
    TP_Tone_StartStop{i_TD} = NaN(36,2);
    TP_Tone_StartStop{i_TD}(1,1) = 0; %p1 set as t = 0 in trial definition
    for i_tone = 2:36
        Dist = abs(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} - ((str2num(curr_toneDur_text)*i_tone) - str2num(curr_toneDur_text)));
        minDist = min(Dist);
        i_minDist = find(Dist == minDist);
        TP_Tone_StartStop{i_TD}(i_tone,1) = DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(i_minDist); %Tone start
    end
    
    epochdur_ms = param.GC.epochdur_ms;
    if strcmp(param.GC.epochdur_ms, 'full')
        %Option 1: End of tone X is start of tone X+1
        for i_tone = 1:35
            i_LastSampleTone = find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} ...
                == TP_Tone_StartStop{i_TD}(i_tone + 1,1));
            TP_Tone_StartStop{i_TD}(i_tone,2) = ...
                DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(i_LastSampleTone); %Tone end
        end
    else
        %Option 2: End of tone X is Y s after start of tone X (useful if one
        %does not want to include the entire tone in analysis window)
        for i_tone = 1:34
            Dist = abs(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} - ...
                (TP_Tone_StartStop{i_TD}(i_tone,1) + (epochdur_ms/1000)));
            minDist = min(Dist);
            i_minDist = find(Dist == minDist);
            TP_Tone_StartStop{i_TD}(i_tone,2) = ...
                DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(i_minDist); %Tone end
        end
    end
    
    TP_Tone_StartStop{i_TD} = TP_Tone_StartStop{i_TD}(1:34,:);
    %Check if all tones are of equal length
    %(if not then choose min length by deleting the last sample of longer trials)
    minSeqLength_sec = min(TP_Tone_StartStop{i_TD}(:,2)-TP_Tone_StartStop{i_TD}(:,1));
    for i_tone = 1:34
        if (TP_Tone_StartStop{i_TD}(i_tone,2) - TP_Tone_StartStop{i_TD}(i_tone,1)) > minSeqLength_sec
            TP_Tone_StartStop{i_TD}(i_tone,2) = ...
                DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(...
                find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == TP_Tone_StartStop{i_TD}(i_tone,2))-1);
        end
    end
    %Determine samples corresponding to TP
    Sample_Tone_StartStop{i_TD} = NaN(34,2);
    for i_tone = 1:34
        Sample_Tone_StartStop{i_TD}(i_tone,1) = ...
            find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == TP_Tone_StartStop{i_TD}(i_tone,1));
        Sample_Tone_StartStop{i_TD}(i_tone,2) = ...
            find(TP_Tone_StartStop{i_TD}(i_tone,2) == DataClean_CleanTrialsElecs_perTD{i_TD}.time{1});
    end
    %Read out tone groups
    disp(' -- Computing GC for following tone aggregations: --');
    for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
        disp(['Tone group ' num2str(i_tonegroup) ' = tone ' ...
            num2str(param.GC.tone_aggregation{i_tonegroup}(1)) ' to ' ...
            num2str(param.GC.tone_aggregation{i_tonegroup}(end)) ' (' ...
            num2str(TP_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1)) ' - ' ...
            num2str(TP_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2)) 's | sample '...
            num2str(Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1)) ' - ' ...
            num2str(Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2)) ')']) ;
    end
end

%1.5 Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus (i.e., tone sequence start) baseline
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

%1.6 Distinguish trials with different PE magnitudes

% DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.finalID %label p34
% unique(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_final) pitch p34 in log(Hz)
% unique(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.series_f{1}(34)) pitch p34 in Hz

% DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.predID %label p*34
% unique(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_pred) pitch p*34 in log(Hz)
% unique(exp(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_pred)) pitch p*34 in Hz

%Simple PE: p34-p33 (in log(Hz)) - usually 2 different magnitudes
for i_TD = 1:length(ToneDur_text)
    
    %Compute absolute p34-p33 pitch distance (in log(Hz)) per trial
    AbsPitchDiff_p34p33 = [];
    for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
        AbsPitchDiff_p34p33(1,i_trial) = ...
            round(abs(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_final(i_trial) - ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_pen(i_trial)),4);
    end
    % hist(AbsPitchDiff_p34p33)
    
    %Assign index to different PE magnitudes
    Unique_AbsPitchDiff_p34p33 = unique(AbsPitchDiff_p34p33);
    Index_AbsPitchDiff_p34p33_pertrial{i_TD} = ...
        nan(1,length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial));
    
    if length(Unique_AbsPitchDiff_p34p33) > 2
        disp(['-- More than 2 unique p34p33 amps -- ']) %Case for NY688 because altered paradigm
        for i_unique_diff = 1:length(Unique_AbsPitchDiff_p34p33)
            if i_unique_diff < 3
                Index_AbsPitchDiff_p34p33_pertrial{i_TD}(...
                    AbsPitchDiff_p34p33 == ...
                    Unique_AbsPitchDiff_p34p33(i_unique_diff)) ...
                    = 1;
            elseif i_unique_diff >= 3
                Index_AbsPitchDiff_p34p33_pertrial{i_TD}(...
                    AbsPitchDiff_p34p33 == ...
                    Unique_AbsPitchDiff_p34p33(i_unique_diff)) ...
                    = 2;
            end
        end
    else
        disp(['-- Exactly 2 unique p34p33 amps -- ']) %Standard Case
        for i_unique_diff = 1:length(Unique_AbsPitchDiff_p34p33)
            if i_unique_diff < 2
                Index_AbsPitchDiff_p34p33_pertrial{i_TD}(...
                    AbsPitchDiff_p34p33 == ...
                    Unique_AbsPitchDiff_p34p33(i_unique_diff)) ...
                    = 1;
            elseif i_unique_diff < 3
                Index_AbsPitchDiff_p34p33_pertrial{i_TD}(...
                    AbsPitchDiff_p34p33 == ...
                    Unique_AbsPitchDiff_p34p33(i_unique_diff)) ...
                    = 2;
            end
        end
    end
    numtrials_lowampp34p33 = ...
        length(find(Index_AbsPitchDiff_p34p33_pertrial{i_TD} == 1));
    numtrials_highampp34p33 = ...
        length(find(Index_AbsPitchDiff_p34p33_pertrial{i_TD} == 2));
    
    % hist(Index_AbsPitchDiff_p34p33_pertrial{i_TD})
    disp(['-- Found ' num2str(length(Unique_AbsPitchDiff_p34p33)) ...
        ' different simple PE magnitudes with ' ...
        num2str(numtrials_lowampp34p33) ' low / ' ...
        num2str(numtrials_highampp34p33) ' high amp trials ' ...
        'for ' sub ' --'])
    
    %Add index to subfield so that subsequent trial selection/removal is taken into account
    DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34p33 = ...
        Index_AbsPitchDiff_p34p33_pertrial{i_TD};
    
    clear Unique_AbsPitchDiff_p34p33 AbsPitchDiff_p34p33  Index_AbsPitchDiff_p34p33_pertrial
    
end

%Complex PE: p34 - p*34 (in log(Hz); p*34 discretized) - same computation for PE analysis
for i_TD = 1:length(ToneDur_text)
    
    %Compute absolute p34-p*34 pitch distance (in log(Hz)) per trial
    AbsPitchDiff_p34Predp34 = [];
    for i_trial = 1:length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial)
        AbsPitchDiff_p34Predp34(1,i_trial) = ...
            round(abs(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_final(i_trial) - ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.logf_pred(i_trial)),4);
    end
    % hist(AbsPitchDiff_p34Predp34)
    
    %Assign index to different PE magnitudes
    Unique_AbsPitchDiff_p34Predp34 = unique(AbsPitchDiff_p34Predp34);
    Index_AbsPitchDiff_p34Predp34_pertrial{i_TD} = ...
        nan(1,length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial));
    
    if length(Unique_AbsPitchDiff_p34Predp34) > 5
        disp(['-- More than 5 unique p34Predp34 amps -- ']) %Case for NY688 because altered paradigm
        for i_unique_diff = 1:length(Unique_AbsPitchDiff_p34Predp34)
            if i_unique_diff < 3 %low PE magnitude (0.0578 ,0.2310 ,0.2888 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 1;
            elseif i_unique_diff < 7  %medium PE magnitude (0.4043, 0.4621, 0.5199 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 2;
            elseif i_unique_diff >= 7 %high PE magnitude (0.6354, 0.6931, 0.8664 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 3;
            end
        end
    else
        disp(['-- Exactly 5 unique p34Predp34 amps -- '])%Standard Case
        for i_unique_diff = 1:length(Unique_AbsPitchDiff_p34Predp34)
            if i_unique_diff <= 2 %low PE magnitude (0.1733 & 0.3466 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 1;
            elseif i_unique_diff == 3 %medium PE magnitude (0.5199 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 2;
            elseif i_unique_diff > 3 %high PE magnitude (0.6931 & 0.8664 log(Hz))
                Index_AbsPitchDiff_p34Predp34_pertrial{i_TD}(...
                    AbsPitchDiff_p34Predp34 == ...
                    Unique_AbsPitchDiff_p34Predp34(i_unique_diff)) ...
                    = 3;
            end
        end
    end
    % hist(Index_AbsPitchDiff_p34Predp34_pertrial{i_TD})
    numtrials_lowampp34Predp34 = ...
        length(find(Index_AbsPitchDiff_p34Predp34_pertrial{i_TD} == 1));
    numtrials_medampp34Predp34 = ...
        length(find(Index_AbsPitchDiff_p34Predp34_pertrial{i_TD} == 2));
    numtrials_highampp34Predp34 = ...
        length(find(Index_AbsPitchDiff_p34Predp34_pertrial{i_TD} == 3));
    
    % hist(Index_AbsPitchDiff_p34p33_pertrial{i_TD})
    disp(['-- Found ' num2str(length(Unique_AbsPitchDiff_p34Predp34)) ...
        ' different complex PE magnitudes with ' ...
        num2str(numtrials_lowampp34Predp34) ' low / ' ...
        num2str(numtrials_medampp34Predp34) ' med / ' ...
        num2str(numtrials_highampp34Predp34) ' high amp trials ' ...
        'for ' sub ' --'])
    
    %Add index to subfield so that subsequent trial selection/removal is taken into account
    DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34Predp34 = ...
        Index_AbsPitchDiff_p34Predp34_pertrial{i_TD};
    
    clear Unique_AbsPitchDiff_p34Predp34 Index_AbsPitchDiff_p34Predp34_pertrial Index_AbsPitchDiff_p34Predp34_pertrial
    
end

%1.7 Remove trials where some electrodes have only NaN samples (NY798)
for i_TD = 1:length(ToneDur_text)
    DataClean_CleanTrialsElecs_perTD{i_TD} = ...
        NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
        (sub, ToneDur_text{i_TD}, DataClean_CleanTrialsElecs_perTD{i_TD});
end

%1.8 Apply zero-phase-shift filter to get input data
tic
if strcmp(param.GC.InputDataType, 'Broadband') %1) broadband data
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP150Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP150Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            0, 0, paths_NASTD_ECoG);
    end
elseif strcmp(param.GC.InputDataType, 'HP05toLP30Hz') %2) slow freqs (0.5-30 Hz)
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            0, 0, paths_NASTD_ECoG);
    end
elseif strcmp(param.GC.InputDataType, 'HighGamma_LogAmp') %3) high gamma-band amplitude (70-150Hz)
    FrequencyBands          = [70, 150]; %Hz [HP,LP]
    for i_TD = 1:length(ToneDur_text)
        FieldLabel = param.GC.InputDataType{1};
        [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel...
            (sub, ToneDur_text{i_TD}, ...
            FrequencyBands(1), FrequencyBands(2), FieldLabel,...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            0, 0, paths_NASTD_ECoG);
    end
end
disp(['-- Filtering & NaN correction took ' num2str(round(toc/60),2) ' minutes --'])

%1.9 Read out input data parameters
for i_TD = 1:length(ToneDur_text)
    numTrials(i_TD) = ...
        size(DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata,2);
    
    temp_numSamplesperTrial = [];
    for i_trial = 1:numTrials(i_TD)
        temp_numSamplesperTrial(i_trial) = ...
            length(DataClean_CleanTrialsElecs_perTD{i_TD}.time{i_trial});
    end
    numSamplesperTrial(i_TD) = min(temp_numSamplesperTrial);
end
numChans        = size(DataClean_CleanTrialsElecs_perTD{1}.label,1);
label_allelecs  = DataClean_CleanTrialsElecs_perTD{i_TD}.label;

clear temp*


%% Simple PE GC computation
%% 2) Select and pair electrodes
%Determine preproc data indices of selected paired electrodes
clear validF validA temporalGC temporalGC_pval frequencyGC

for i_effect = 1:length(param.GC.ElecPairEffect)
    tstart = tic;
    
    if SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
            param.GC.ElecPairRegion).Num_UniquePairings_persub(i_sub) > 0
        
        %Determine effect label
        if strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_Pred')
            label_effect1 = 'PredEffect';
            label_effect2 = 'PredEffect';
            title_label1   = 'Pred2Pred';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_PE')
            label_effect1 = 'PEEffect';
            label_effect2 = 'PEEffect';
            title_label1   = 'PE2PE';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_PE')
            label_effect1 = 'PredEffect';
            label_effect2 = 'PEEffect';
            title_label1   = 'Pred2PE';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_Pred')
            label_effect1 = 'PEEffect';
            label_effect2 = 'PredEffect';
            title_label1   = 'PE2Pred';
        end
        
        if strcmp(param.GC.ElecPairRegion, 'AllRegions')
            title_label2   = 'AllReg';
        elseif strcmp(param.GC.ElecPairRegion, 'Frontal_Temporal')
            title_label2   = 'Front2Temp';
        elseif strcmp(param.GC.ElecPairRegion, 'Temporal_Frontal')
            title_label2   = 'Temp2Front';
        end
        
        %Determine source and target elec indeces as listed in preproc data (for data selection)
        
        %Note: In SelElecs input data, the elec pairs are column-separated
        %into source vs. target, with reverse-directional double entries removed.
        %For Pred_Pred/PE_PE, however, this is problematic, since this
        %leads to a missing out of a few connections, because we need to
        %take bidirectional connections into account here. Thus, these
        %effect types require a workaround, where we read out the unique
        %indices and redefine the pairs
        if strcmp(param.GC.ElecPairEffect{i_effect},  'Pred_Pred') || ...
                strcmp(param.GC.ElecPairEffect{i_effect},  'PE_PE')
            ind_pairedelecs_from{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub},...
                'stable'), ...
                10}; %Take all elec indices from both columns
            ind_pairedelecs_to{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub},...
                'stable'), ...
                10}; %Take all elec indices from both columns
            
        elseif strcmp(param.GC.ElecPairEffect{i_effect},  'Pred_PE') || ...
                strcmp(param.GC.ElecPairEffect{i_effect},  'PE_Pred')
            ind_pairedelecs_from{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1),...
                'stable'),...
                10};%Take all indices from first (source) column only
            ind_pairedelecs_to{i_effect} = ...
                SelElecs.(label_effect2)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2),...
                'stable'),...
                10};%Take all indices from second (target) column only
        end
        
        ind_pairedelecs_both    = [ind_pairedelecs_from{i_effect}; ind_pairedelecs_to{i_effect}];
        label_pairedelecs_both  = label_allelecs(ind_pairedelecs_both);
        
        %Sanity Check if labels from pairing and preproc data match
        clear ind2_* unique_ind*
        if strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_Pred') || ...
                strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_PE')
            %Determine elec index as listed in elec pairs data (for label comparison)
            ind2_pairedelecs        = ...
                [unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1)); ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2))];
            unique_ind2_pairedelecs = unique(ind2_pairedelecs, 'stable'); %elec indices for SelElecs
            unique_ind_pairedelecs = unique(ind_pairedelecs_both, 'stable');% elec indices for preproc data
            for i_selelec = 1:length(unique_ind_pairedelecs)
                if strcmp(...
                        label_allelecs(unique_ind_pairedelecs(i_selelec)), ...
                        SelElecs.(label_effect1){unique_ind2_pairedelecs(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_both(i_selelec))])
                    pause
                end
            end
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_PE') || ...
                strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_Pred')
            %Determine elec index as listed in elec pairs data (for label comparison)
            ind2_pairedelecs1 = ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1),...
                'stable');
            ind2_pairedelecs2 = ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2),...
                'stable');
            for i_selelec = 1:length(ind_pairedelecs_from{i_effect})
                if strcmp(...
                        label_allelecs(ind_pairedelecs_from{i_effect}(i_selelec)), ...
                        SelElecs.(label_effect1){ind2_pairedelecs1(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_from{i_effect}(i_selelec))])
                    pause
                end
            end
            for i_selelec = 1:length(ind_pairedelecs_to{i_effect})
                if strcmp(...
                        label_allelecs(ind_pairedelecs_to{i_effect}(i_selelec)), ...
                        SelElecs.(label_effect2){ind2_pairedelecs2(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_to{i_effect}(i_selelec))])
                    pause
                end
            end
        end
        
        %Provide elec labels and display current elec info
        label_sourceelecs = [];
        for i_elec = 1:length(label_allelecs(ind_pairedelecs_from{i_effect}))
            label_sourceelecs = [label_sourceelecs, char(label_allelecs(ind_pairedelecs_from{i_effect}(i_elec))) ', ' ];
        end
        label_targetelecs = [];
        for i_elec = 1:length(label_allelecs(ind_pairedelecs_to{i_effect}))
            label_targetelecs = [label_targetelecs, char(label_allelecs(ind_pairedelecs_to{i_effect}(i_elec))) ', ' ];
        end
        disp(['-- Computing GC for ' sub ' for effect ' param.GC.ElecPairEffect{i_effect} ' (' ...
            param.GC.ElecPairRegion ')'])
        disp(['-- from source elecs (' label_sourceelecs ')'])
        disp(['-- to target elecs (' label_targetelecs ')' ])
        
        
        %% 3) Set up empty structs and select elecs and samples for GC input data
        
        %Since trial-specific analyses start here, we differentiate by
        %different simple/complex PE trials.
        
        for i_TD = 1:length(ToneDur_text)
            for i_simplePEeffect = 1:length(unique(...        %Simple PE trial differentiation          
                    DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34p33))
                                
                %3.1 Bring GC input data into requested format (dim: channels, timepoints, trials;)
                trialData = [];
                i_trialselection_simplePE{i_simplePEeffect} = ...
                    find(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34p33 ...
                    == i_simplePEeffect);
                
                for i_trial = 1:length(i_trialselection_simplePE{i_simplePEeffect})
                    trialData(:,:,i_trial) = ...
                        DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata{i_trialselection_simplePE{i_simplePEeffect}(i_trial)}...
                        (:,1:numSamplesperTrial(i_TD));
                end
                
                %Restrict input data to selected electrodes and time frame
                trialData = trialData(ind_pairedelecs_both, :, :);
                
                %Set structs for Output parameters (dim: source elec * target elec)
                validA{i_effect}{i_TD}{i_simplePEeffect} = ... %Indicates if VAR coefficients matrix is valid
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                validF{i_effect}{i_TD}{i_simplePEeffect} = ...%Indicates if GC computation is valid
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                
                temporalGC{i_effect}{i_TD}{i_simplePEeffect}.source2target = ...%pairwise-conditional time-domain multivariate Granger causalities, source->target
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.source2target = ...%pairwise-conditional frequency-domain/ spectral multivariate Granger causalites, source->target
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC{i_effect}{i_TD}{i_simplePEeffect}.target2source = ...%pairwise-conditional time-domain multivariate Granger causalities, target->source
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.target2source = ...%pairwise-conditional frequency-domain/ spectral multivariate Granger causalites, target->source
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                
                frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.source2target = ...
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    (param.GC.fs)+1, ...
                    size(param.GC.tone_aggregation,2));
                frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.target2source = ...
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    (param.GC.fs)+1, ...
                    size(param.GC.tone_aggregation,2));                
                
                %dim: source elec * target elec * tone aggregation
            
            %% 5) Compute GC (for each tone snippet in each selected electrode pairing)
            %5.1 Preprocess to increase stationarity:
            %From Hardstone et al. (2021) and Chao et al. (2018)
            clear trialData_detrend trialData_tempnorm trialData_ensemble*
                for i_elec = 1:size(trialData,1)
                    for i_trial = 1:size(trialData,3)
                        %1) Detrend data for each trial, electrode, TD
                        trialData_detrend(i_elec,:,i_trial) = ...
                            detrend(trialData(i_elec,:,i_trial));
                        %2) Temporal normalization per trial, electrode, TD
                        %(subtraction of the across-sample mean of each
                        %time series and division by the across-sample standard deviation)
                        %Ensures that all variables have equal weights across the trial.
                        trialData_tempnorm(i_elec,:,i_trial) = ...
                            trialData_detrend(i_elec,:,i_trial) - ...
                            nanmean(trialData_detrend(i_elec,:,i_trial));
                        trialData_tempnorm...
                            (i_elec,:,i_trial) = ...
                            trialData_tempnorm(i_elec,:,i_trial) ./ ...
                            nanstd(trialData_tempnorm(i_elec,:,i_trial));
                    end
                    %3)Ensemble normalization per trial, electrode, TD
                    %(pointwise subtraction of the ensemble mean, calculated
                    %by averaging the measured amplitude values at a given channel at each
                    %time point across trials) and division by the ensemble standard deviation,
                    for i_trial = 1:size(trialData,3)
                        trialData_ensemble(i_elec,:,i_trial) = ...
                            trialData_tempnorm(i_elec,:,i_trial) - ...
                            nanmean(squeeze(trialData_tempnorm(i_elec,:,:)),2)';
                    end
                    for i_trial = 1:size(trialData,3)
                        trialData_ensemble2(i_elec,:,i_trial) = ...
                            trialData_ensemble(i_elec,:,i_trial) ./ ...
                            nanstd(squeeze(trialData(i_elec,:,:)),0,2)';
                    end
                end
            
                %4.1) determine appropriate model order for VAR model (tsdata_to_infocrit.m)
                %(balance maximum lag p to achieve best model fit while avoiding overfitting)
                for i_sourceelec = 1:length(ind_pairedelecs_from{i_effect})
                    for i_targetelec = 1:length(ind_pairedelecs_to{i_effect})
                        
                        %Avoid double-dipping (identical electrode as source + target)
                        if ind_pairedelecs_from{i_effect}(i_sourceelec) ~= ind_pairedelecs_to{i_effect}(i_targetelec)
                            
                            [label_allelecs{ind_pairedelecs_from{i_effect}(i_sourceelec)} ...
                                ' -> ' ...
                                label_allelecs{ind_pairedelecs_to{i_effect}(i_targetelec)}]
                            
                            %Calculate standerdized information criteria up to a priori specified maximum model order.
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                ptic('\n*** Determine model order (tsdata_to_infocrit)\n');
                                [AIC, BIC, ...
                                    moAIC(i_sourceelec, i_targetelec, i_tonegroup), ...
                                    moBIC(i_sourceelec, i_targetelec, i_tonegroup)] = ...
                                    tsdata_to_infocrit(...
                                    trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), ...
                                    :), ...
                                    param.GC.maxmorder, param.GC.regmode);
                                ptoc('*** tsdata_to_infocrit took ');
                                % figure(1); clf;
                                % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/param.GC.fs);
                                % title('Model order estimation');
                            end
                        else
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                moAIC(i_sourceelec, i_targetelec, i_tonegroup) = NaN;
                                moBIC(i_sourceelec, i_targetelec, i_tonegroup) = NaN;
                            end
                        end
                    end
                end
                
                %Select median model order across elec pairs * tone snippets as final model order
                %(similar to Hardston et al., Nat Commun 2021)
                param.GC.SimplePE.morder_pereffectTD(i_effect,i_TD,i_simplePEeffect) = floor(nanmedian(moBIC, 'all'));
                fprintf('\nbest model order (BIC) = %d\n',param.GC.SimplePE.morder_pereffectTD(i_effect,i_TD,i_simplePEeffect));
                
                %     %Optional: Sample autocovariance estimation
                %     G = tsdata_to_autocov(trialData_ensemble2,param.GC.morder_pereffectTD(i_effect,i_TD))
                
                num_Pairs{i_effect}{i_TD}{i_simplePEeffect}           = 0;
                num_InvalidPairs{i_effect}{i_TD}{i_simplePEeffect}    = 0;
                
                for i_sourceelec = 1:length(ind_pairedelecs_from{i_effect})
                    for i_targetelec = 1:length(ind_pairedelecs_to{i_effect})
                        
                        %Avoid double-dipping (identical electrode as source + target)
                        if ind_pairedelecs_from{i_effect}(i_sourceelec) ~= ind_pairedelecs_to{i_effect}(i_targetelec)
                            disp(['-- Comuting Simple PE GC for '...
                                label_pairedelecs_both{i_sourceelec} ' -> ' ...
                                label_pairedelecs_both{length(ind_pairedelecs_from{i_effect}) + i_targetelec}]);
                            
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                AS = struct; %Struct for VAR model (contains: A = VAR coefficients matrix, SIG = residuals covariance matrix)
                                FS = struct; %Struct for pairwise-conditional time-domain multivariate Granger causalities
                                
                                %4.2 Fit VAR model to time series data (for up to selected model order)
                                %(obtain estimates of model parameters which maximize the likelihood
                                %function for both full (i.e., X has influence on Y) and reduced
                                %(i.e., X has no influence on Y) VAR models - i.e., model fitting)
                                ptic('\n*** estimating VAR model (tsdata_to_var)... ');
                                [AS.A, AS.SIG] = ...
                                    tsdata_to_var...
                                    (trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), :), ...
                                    param.GC.SimplePE.morder_pereffectTD(i_effect,i_TD,i_simplePEeffect), param.GC.regmode);
                                
                                %                 assert(~isbad(AS.A),'VAR estimation failed - bailing out');
                                ptoc;
                                %A = VAR coefficients matrix (dim: elec * elec * lag/model order)
                                %SIG = residuals covariance matrix (dim: elec * elec)
                                
                                %Check if VAR model estimation worked out (check for rank
                                %defficient VAR output matrix due to colinearlity in input data
                                %and that VAR spectral radius < 1, indicating stationarity or
                                %long-term memory of input data)
                                info = var_info(AS.A, AS.SIG);
                                %                 assert(~info.error,'VAR error(s) found - bailing out');
                                validA{i_effect}{i_TD}{i_simplePEeffect}(i_sourceelec, i_targetelec, i_tonegroup) = ...
                                    ~isbad(AS.A); %Record potential errors
                                
                                %4.3 Directly compute GC and corresponding p-value from full
                                %(i.e., X has influence on Y) and reduced/null (i.e., X has no influence on Y)
                                %VAR models using novel State-space version (not degraded by
                                %moving average components in input data)
                                ptic('*** computing time domain GC estimates and testing against reduced model (var_to_pwcgc)... ');
                                [FS.F, FS.pval] = ...
                                    var_to_pwcgc(...
                                    AS.A, AS.SIG, ...
                                    trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), :), ...
                                    param.GC.regmode, param.GC.tstat);
                                %F  = Time-domain multivariate Granger causality: for the pairwise-
                                %conditional case, an n x n matrix containing pairwise-conditional causalities;
                                %the first index then references target ("to") variables
                                %and the second source ("from") variables,
                                %with NaNs on the diagonal.
                                
                                %F = pairwise conditional GC estimates (dim: elec*elec)
                                %pval = p-value from comparison full model vs. reduced/null model
                                %(i.e., does the past of source elec contain information
                                %about the future of target elec, above and beyond information
                                %contained in the past of the target elec)
                                ptoc;
                                
                                % Check for failed GC calculation
                                % assert(~isbad(FS.F,false),'GC calculation failed - bailing out');
                                validF{i_effect}{i_TD}{i_simplePEeffect}...
                                    (i_sourceelec,i_targetelec, i_tonegroup) = ...
                                    ~isbad(FS.F,false);
                                
                                %4.4 Calculate spectral pairwise-conditional causalities at given frequency
                                % resolution by state-space method.
                                ptic('*** computing frequency domain GC estimates  (var_to_spwcgc)... ');
                                FS.f = var_to_spwcgc...
                                    (AS.A, AS.SIG, ...
                                    param.GC.newfs);
                                % assert(~isbad(FS.f,false),'spectral GC calculation failed - bailing out');
                                ptoc;
                                
                                %4.5 Copy results in output struct
                                if validF{i_effect}{i_TD}{i_simplePEeffect}(i_sourceelec, i_targetelec, i_tonegroup) == 1
                                    
                                    %Richard (05.10.) 
                                    %F(1,2) = chan2(here i_targetelec) -> chan1(here i_sourceelec)
                                    %F(2,1) = chan1(here i_sourceelec ) -> chan2(here i_targetelec)

                                    %Direction source -> target
                                    temporalGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Source elec, Target elec, Tone)
                                        = FS.F(2,1);                                    
                                    temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Source elec, Target elec, Tone)
                                        = FS.pval(2,1);
                                    frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...%(Source elec, Target elec, Frequency, Tone)
                                        = FS.f(2,1,:);
                                    
                                    %Direction target -> source
                                    temporalGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Target elec, Source elec, Tone)
                                        = FS.F(1,2);
                                    temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Target elec, Source elec, Tone)
                                        = FS.pval(1,2);
                                    frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...%(Target elec, Source elec, Frequency, Tone)
                                        = FS.f(1,2,:);
                                    
                                    num_Pairs{i_effect}{i_TD}{i_simplePEeffect} = ...
                                        num_Pairs{i_effect}{i_TD}{i_simplePEeffect} + 1;
                                    
                                else
                                    temporalGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                        = nan;
                                    temporalGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                        = nan;
                                    
                                    num_Pairs{i_effect}{i_TD}{i_simplePEeffect} = ...
                                        num_Pairs{i_effect}{i_TD}{i_simplePEeffect} + 1;
                                    num_InvalidPairs{i_effect}{i_TD}{i_simplePEeffect} = ...
                                        num_InvalidPairs{i_effect}{i_TD}{i_simplePEeffect} + 1;
                                end
                            end
                        else %in case of double dipping, NaN results but don't count as invalid
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                temporalGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                    = nan;
                                temporalGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                temporalGC_pval{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                frequencyGC{i_effect}{i_TD}{i_simplePEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                    = nan;
                            end
                        end
                    end
                end
                disp(['GC computation for ' ToneDur_text{i_TD} 's TD finished - ' ...
                    num2str(num_Pairs{i_effect}{i_TD}{i_simplePEeffect} - ...
                    num_InvalidPairs{i_effect}{i_TD}{i_simplePEeffect}) '/' ...
                    num2str(num_Pairs{i_effect}{i_TD}{i_simplePEeffect}) ...
                    ' valid (GC computation did not bail out)']);
            end
        end
    end
    
    disp(['-- Computing GC from source elecs to target elecs for ' ...
        sub ' and effect ' ...
        param.GC.ElecPairEffect{i_effect} ' (' ...
        param.GC.ElecPairRegion ') took ' ...
        num2str(round(toc(tstart),0)) ' secs --'])
    
end

%% 5) Copy results in common output file
GCdata.SimplePE.dim                  = 'ElecSelectEffect_TD_SimplePEEffect';
GCdata.SimplePE.temporalGC           = temporalGC;
GCdata.SimplePE.pval_temporalGC      = temporalGC_pval;
GCdata.SimplePE.spectralGC           = frequencyGC;
GCdata.SimplePE.trialindices         = i_trialselection_simplePE;

GCdata.SimplePE.info.validVAR        = validA;
GCdata.SimplePE.info.validGC         = validF;
GCdata.SimplePE.info.numPairs        = num_Pairs;
GCdata.SimplePE.info.numInvalidPairs = num_InvalidPairs;
GCdata.SimplePE.info.morder          = param.GC.SimplePE.morder_pereffectTD;

GCdata.elecs.label_allelecs         = label_allelecs;
GCdata.elecs.ind_pairedelecs_from   = ind_pairedelecs_from;
GCdata.elecs.ind_pairedelecs_to     = ind_pairedelecs_to;


%% Complex PE GC computation
%% 2) Select and pair electrodes
%Determine preproc data indices of selected paired electrodes
clear validF validA temporalGC temporalGC_pval frequencyGC

for i_effect = 1:length(param.GC.ElecPairEffect)
    tstart = tic;
    
    if SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
            param.GC.ElecPairRegion).Num_UniquePairings_persub(i_sub) > 0
        
        %Determine effect label
        if strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_Pred')
            label_effect1 = 'PredEffect';
            label_effect2 = 'PredEffect';
            title_label1   = 'Pred2Pred';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_PE')
            label_effect1 = 'PEEffect';
            label_effect2 = 'PEEffect';
            title_label1   = 'PE2PE';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_PE')
            label_effect1 = 'PredEffect';
            label_effect2 = 'PEEffect';
            title_label1   = 'Pred2PE';
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_Pred')
            label_effect1 = 'PEEffect';
            label_effect2 = 'PredEffect';
            title_label1   = 'PE2Pred';
        end
        
        if strcmp(param.GC.ElecPairRegion, 'AllRegions')
            title_label2   = 'AllReg';
        elseif strcmp(param.GC.ElecPairRegion, 'Frontal_Temporal')
            title_label2   = 'Front2Temp';
        elseif strcmp(param.GC.ElecPairRegion, 'Temporal_Frontal')
            title_label2   = 'Temp2Front';
        end
        
        %Determine source and target elec indeces as listed in preproc data (for data selection)
        
        %Note: In SelElecs input data, the elec pairs are column-separated
        %into source vs. target, with reverse-directional double entries removed.
        %For Pred_Pred/PE_PE, however, this is problematic, since this
        %leads to a missing out of a few connections, because we need to
        %take bidirectional connections into account here. Thus, these
        %effect types require a workaround, where we read out the unique
        %indices and redefine the pairs
        if strcmp(param.GC.ElecPairEffect{i_effect},  'Pred_Pred') || ...
                strcmp(param.GC.ElecPairEffect{i_effect},  'PE_PE')
            ind_pairedelecs_from{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub},...
                'stable'), ...
                10}; %Take all elec indices from both columns
            ind_pairedelecs_to{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub},...
                'stable'), ...
                10}; %Take all elec indices from both columns
            
        elseif strcmp(param.GC.ElecPairEffect{i_effect},  'Pred_PE') || ...
                strcmp(param.GC.ElecPairEffect{i_effect},  'PE_Pred')
            ind_pairedelecs_from{i_effect} = ...
                SelElecs.(label_effect1)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1),...
                'stable'),...
                10};%Take all indices from first (source) column only
            ind_pairedelecs_to{i_effect} = ...
                SelElecs.(label_effect2)...
                {unique(...
                SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2),...
                'stable'),...
                10};%Take all indices from second (target) column only
        end
        
        ind_pairedelecs_both    = [ind_pairedelecs_from{i_effect}; ind_pairedelecs_to{i_effect}];
        label_pairedelecs_both  = label_allelecs(ind_pairedelecs_both);
        
        %Sanity Check if labels from pairing and preproc data match
        clear ind2_* unique_ind*
        if strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_Pred') || ...
                strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_PE')
            %Determine elec index as listed in elec pairs data (for label comparison)
            ind2_pairedelecs        = ...
                [unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1)); ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2))];
            unique_ind2_pairedelecs = unique(ind2_pairedelecs, 'stable'); %elec indices for SelElecs
            unique_ind_pairedelecs = unique(ind_pairedelecs_both, 'stable');% elec indices for preproc data
            for i_selelec = 1:length(unique_ind_pairedelecs)
                if strcmp(...
                        label_allelecs(unique_ind_pairedelecs(i_selelec)), ...
                        SelElecs.(label_effect1){unique_ind2_pairedelecs(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_both(i_selelec))])
                    pause
                end
            end
        elseif strcmp(param.GC.ElecPairEffect{i_effect}, 'Pred_PE') || ...
                strcmp(param.GC.ElecPairEffect{i_effect}, 'PE_Pred')
            %Determine elec index as listed in elec pairs data (for label comparison)
            ind2_pairedelecs1 = ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,1),...
                'stable');
            ind2_pairedelecs2 = ...
                unique(SelElecs.Pairs.(param.GC.ElecPairEffect{i_effect}).(...
                param.GC.ElecPairRegion).Ind_UniquePairings_persub{i_sub}(:,2),...
                'stable');
            for i_selelec = 1:length(ind_pairedelecs_from{i_effect})
                if strcmp(...
                        label_allelecs(ind_pairedelecs_from{i_effect}(i_selelec)), ...
                        SelElecs.(label_effect1){ind2_pairedelecs1(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_from{i_effect}(i_selelec))])
                    pause
                end
            end
            for i_selelec = 1:length(ind_pairedelecs_to{i_effect})
                if strcmp(...
                        label_allelecs(ind_pairedelecs_to{i_effect}(i_selelec)), ...
                        SelElecs.(label_effect2){ind2_pairedelecs2(i_selelec),1}) == false
                    disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_to{i_effect}(i_selelec))])
                    pause
                end
            end
        end
        
        %Provide elec labels and display current elec info
        label_sourceelecs = [];
        for i_elec = 1:length(label_allelecs(ind_pairedelecs_from{i_effect}))
            label_sourceelecs = [label_sourceelecs, char(label_allelecs(ind_pairedelecs_from{i_effect}(i_elec))) ', ' ];
        end
        label_targetelecs = [];
        for i_elec = 1:length(label_allelecs(ind_pairedelecs_to{i_effect}))
            label_targetelecs = [label_targetelecs, char(label_allelecs(ind_pairedelecs_to{i_effect}(i_elec))) ', ' ];
        end
        disp(['-- Computing GC for ' sub ' for effect ' param.GC.ElecPairEffect{i_effect} ' (' ...
            param.GC.ElecPairRegion ')'])
        disp(['-- from source elecs (' label_sourceelecs ')'])
        disp(['-- to target elecs (' label_targetelecs ')' ])
        
        
        %% 3) Set up empty structs and select elecs and samples for GC input data
        
        %Since trial-specific analyses start here, we differentiate by
        %different simple/complex PE trials.
        
        for i_TD = 1:length(ToneDur_text)
            for i_complexPEeffect = 1:length(unique(...        %Simple PE trial differentiation          
                    DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34Predp34))
                                
                %3.1 Bring GC input data into requested format (dim: channels, timepoints, trials;)
                trialData = [];
                i_trialselection_complexPE{i_complexPEeffect} = ...
                    find(DataClean_CleanTrialsElecs_perTD{i_TD}.behav.stim.Index_AbsPitchDiff_p34Predp34 ...
                    == i_complexPEeffect);
                
                for i_trial = 1:length(i_trialselection_complexPE{i_complexPEeffect})
                    trialData(:,:,i_trial) = ...
                        DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata{i_trialselection_complexPE{i_complexPEeffect}(i_trial)}...
                        (:,1:numSamplesperTrial(i_TD));
                end
                
                %Restrict input data to selected electrodes and time frame
                trialData = trialData(ind_pairedelecs_both, :, :);
                
                %Set structs for Output parameters (dim: source elec * target elec)
                validA{i_effect}{i_TD}{i_complexPEeffect} = ... %Indicates if VAR coefficients matrix is valid
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                validF{i_effect}{i_TD}{i_complexPEeffect} = ...%Indicates if GC computation is valid
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                
                temporalGC{i_effect}{i_TD}{i_complexPEeffect}.source2target = ...%pairwise-conditional time-domain multivariate Granger causalities, source->target
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.source2target = ...%pairwise-conditional frequency-domain/ spectral multivariate Granger causalites, source->target
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC{i_effect}{i_TD}{i_complexPEeffect}.target2source = ...%pairwise-conditional time-domain multivariate Granger causalities, target->source
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.target2source = ...%pairwise-conditional frequency-domain/ spectral multivariate Granger causalites, target->source
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    size(param.GC.tone_aggregation,2));
                
                frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.source2target = ...
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    (param.GC.fs)+1, ...
                    size(param.GC.tone_aggregation,2));
                frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.target2source = ...
                    nan(size(ind_pairedelecs_from{i_effect},1), ...
                    size(ind_pairedelecs_to{i_effect},1), ...
                    (param.GC.fs)+1, ...
                    size(param.GC.tone_aggregation,2));                
                
                %dim: source elec * target elec * tone aggregation
            
            %% 5) Compute GC (for each tone snippet in each selected electrode pairing)
            %5.1 Preprocess to increase stationarity:
            %From Hardstone et al. (2021) and Chao et al. (2018)
            clear trialData_detrend trialData_tempnorm trialData_ensemble*
                for i_elec = 1:size(trialData,1)
                    for i_trial = 1:size(trialData,3)
                        %1) Detrend data for each trial, electrode, TD
                        trialData_detrend(i_elec,:,i_trial) = ...
                            detrend(trialData(i_elec,:,i_trial));
                        %2) Temporal normalization per trial, electrode, TD
                        %(subtraction of the across-sample mean of each
                        %time series and division by the across-sample standard deviation)
                        %Ensures that all variables have equal weights across the trial.
                        trialData_tempnorm(i_elec,:,i_trial) = ...
                            trialData_detrend(i_elec,:,i_trial) - ...
                            nanmean(trialData_detrend(i_elec,:,i_trial));
                        trialData_tempnorm...
                            (i_elec,:,i_trial) = ...
                            trialData_tempnorm(i_elec,:,i_trial) ./ ...
                            nanstd(trialData_tempnorm(i_elec,:,i_trial));
                    end
                    %3)Ensemble normalization per trial, electrode, TD
                    %(pointwise subtraction of the ensemble mean, calculated
                    %by averaging the measured amplitude values at a given channel at each
                    %time point across trials) and division by the ensemble standard deviation,
                    for i_trial = 1:size(trialData,3)
                        trialData_ensemble(i_elec,:,i_trial) = ...
                            trialData_tempnorm(i_elec,:,i_trial) - ...
                            nanmean(squeeze(trialData_tempnorm(i_elec,:,:)),2)';
                    end
                    for i_trial = 1:size(trialData,3)
                        trialData_ensemble2(i_elec,:,i_trial) = ...
                            trialData_ensemble(i_elec,:,i_trial) ./ ...
                            nanstd(squeeze(trialData(i_elec,:,:)),0,2)';
                    end
                end
            
                %4.1) determine appropriate model order for VAR model (tsdata_to_infocrit.m)
                %(balance maximum lag p to achieve best model fit while avoiding overfitting)
                for i_sourceelec = 1:length(ind_pairedelecs_from{i_effect})
                    for i_targetelec = 1:length(ind_pairedelecs_to{i_effect})
                        
                        %Avoid double-dipping (identical electrode as source + target)
                        if ind_pairedelecs_from{i_effect}(i_sourceelec) ~= ind_pairedelecs_to{i_effect}(i_targetelec)
                            
                            [label_allelecs{ind_pairedelecs_from{i_effect}(i_sourceelec)} ...
                                ' -> ' ...
                                label_allelecs{ind_pairedelecs_to{i_effect}(i_targetelec)}]
                            
                            %Calculate standerdized information criteria up to a priori specified maximum model order.
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                ptic('\n*** Determine model order (tsdata_to_infocrit)\n');
                                [AIC, BIC, ...
                                    moAIC(i_sourceelec, i_targetelec, i_tonegroup), ...
                                    moBIC(i_sourceelec, i_targetelec, i_tonegroup)] = ...
                                    tsdata_to_infocrit(...
                                    trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), ...
                                    :), ...
                                    param.GC.maxmorder, param.GC.regmode);
                                ptoc('*** tsdata_to_infocrit took ');
                                % figure(1); clf;
                                % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/param.GC.fs);
                                % title('Model order estimation');
                            end
                        else
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                moAIC(i_sourceelec, i_targetelec, i_tonegroup) = NaN;
                                moBIC(i_sourceelec, i_targetelec, i_tonegroup) = NaN;
                            end
                        end
                    end
                end
                
                %Select median model order across elec pairs * tone snippets as final model order
                %(similar to Hardston et al., Nat Commun 2021)
                param.GC.ComplexPE.morder_pereffectTD(i_effect,i_TD,i_complexPEeffect) = floor(nanmedian(moBIC, 'all'));
                fprintf('\nbest model order (BIC) = %d\n',param.GC.ComplexPE.morder_pereffectTD(i_effect,i_TD,i_complexPEeffect));
                
                %     %Optional: Sample autocovariance estimation
                %     G = tsdata_to_autocov(trialData_ensemble2,param.GC.morder_pereffectTD(i_effect,i_TD))
                
                num_Pairs{i_effect}{i_TD}{i_complexPEeffect}           = 0;
                num_InvalidPairs{i_effect}{i_TD}{i_complexPEeffect}    = 0;
                
                for i_sourceelec = 1:length(ind_pairedelecs_from{i_effect})
                    for i_targetelec = 1:length(ind_pairedelecs_to{i_effect})
                        
                        %Avoid double-dipping (identical electrode as source + target)
                        if ind_pairedelecs_from{i_effect}(i_sourceelec) ~= ind_pairedelecs_to{i_effect}(i_targetelec)
                            disp(['-- Comuting Simple PE GC for '...
                                label_pairedelecs_both{i_sourceelec} ' -> ' ...
                                label_pairedelecs_both{length(ind_pairedelecs_from{i_effect}) + i_targetelec}]);
                            
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                AS = struct; %Struct for VAR model (contains: A = VAR coefficients matrix, SIG = residuals covariance matrix)
                                FS = struct; %Struct for pairwise-conditional time-domain multivariate Granger causalities
                                
                                %4.2 Fit VAR model to time series data (for up to selected model order)
                                %(obtain estimates of model parameters which maximize the likelihood
                                %function for both full (i.e., X has influence on Y) and reduced
                                %(i.e., X has no influence on Y) VAR models - i.e., model fitting)
                                ptic('\n*** estimating VAR model (tsdata_to_var)... ');
                                [AS.A, AS.SIG] = ...
                                    tsdata_to_var...
                                    (trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), :), ...
                                    param.GC.ComplexPE.morder_pereffectTD(i_effect,i_TD,i_complexPEeffect), param.GC.regmode);
                                
                                %                 assert(~isbad(AS.A),'VAR estimation failed - bailing out');
                                ptoc;
                                %A = VAR coefficients matrix (dim: elec * elec * lag/model order)
                                %SIG = residuals covariance matrix (dim: elec * elec)
                                
                                %Check if VAR model estimation worked out (check for rank
                                %defficient VAR output matrix due to colinearlity in input data
                                %and that VAR spectral radius < 1, indicating stationarity or
                                %long-term memory of input data)
                                info = var_info(AS.A, AS.SIG);
                                %                 assert(~info.error,'VAR error(s) found - bailing out');
                                validA{i_effect}{i_TD}{i_complexPEeffect}(i_sourceelec, i_targetelec, i_tonegroup) = ...
                                    ~isbad(AS.A); %Record potential errors
                                
                                %4.3 Directly compute GC and corresponding p-value from full
                                %(i.e., X has influence on Y) and reduced/null (i.e., X has no influence on Y)
                                %VAR models using novel State-space version (not degraded by
                                %moving average components in input data)
                                ptic('*** computing time domain GC estimates and testing against reduced model (var_to_pwcgc)... ');
                                [FS.F, FS.pval] = ...
                                    var_to_pwcgc(...
                                    AS.A, AS.SIG, ...
                                    trialData_ensemble2...
                                    ([i_sourceelec length(ind_pairedelecs_from{i_effect})+i_targetelec], ...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(1),1):...
                                    Sample_Tone_StartStop{i_TD}(param.GC.tone_aggregation{i_tonegroup}(end),2), :), ...
                                    param.GC.regmode, param.GC.tstat);
                                %F  = Time-domain multivariate Granger causality: for the pairwise-
                                %conditional case, an n x n matrix containing pairwise-conditional causalities;
                                %the first index then references target ("to") variables
                                %and the second source ("from") variables,
                                %with NaNs on the diagonal.
                                
                                %F = pairwise conditional GC estimates (dim: elec*elec)
                                %pval = p-value from comparison full model vs. reduced/null model
                                %(i.e., does the past of source elec contain information
                                %about the future of target elec, above and beyond information
                                %contained in the past of the target elec)
                                ptoc;
                                
                                % Check for failed GC calculation
                                % assert(~isbad(FS.F,false),'GC calculation failed - bailing out');
                                validF{i_effect}{i_TD}{i_complexPEeffect}...
                                    (i_sourceelec,i_targetelec, i_tonegroup) = ...
                                    ~isbad(FS.F,false);
                                
                                %4.4 Calculate spectral pairwise-conditional causalities at given frequency
                                % resolution by state-space method.
                                ptic('*** computing frequency domain GC estimates  (var_to_spwcgc)... ');
                                FS.f = var_to_spwcgc...
                                    (AS.A, AS.SIG, ...
                                    param.GC.newfs);
                                % assert(~isbad(FS.f,false),'spectral GC calculation failed - bailing out');
                                ptoc;
                                
                                %4.5 Copy results in output struct
                                if validF{i_effect}{i_TD}{i_complexPEeffect}(i_sourceelec, i_targetelec, i_tonegroup) == 1
                                    
                                    %Richard (05.10.) 
                                    %F(1,2) = chan2(here i_targetelec) -> chan1(here i_sourceelec)
                                    %F(2,1) = chan1(here i_sourceelec ) -> chan2(here i_targetelec)

                                    %Direction source -> target
                                        temporalGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Source elec, Target elec, Tone)
                                        = FS.F(2,1);                                    
                                    temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Source elec, Target elec, Tone)
                                        = FS.pval(2,1);
                                    frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...%(Source elec, Target elec, Frequency, Tone)
                                        = FS.f(2,1,:);
                                    
                                    %Direction target -> source
                                    temporalGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Target elec, Source elec, Tone)
                                        = FS.F(1,2);
                                    temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...%(Target elec, Source elec, Tone)
                                        = FS.pval(1,2);
                                    frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...%(Target elec, Source elec, Frequency, Tone)
                                        = FS.f(1,2,:);
                                    
                                    num_Pairs{i_effect}{i_TD}{i_complexPEeffect} = ...
                                        num_Pairs{i_effect}{i_TD}{i_complexPEeffect} + 1;
                                    
                                else
                                    temporalGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                        = nan;
                                    temporalGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,i_tonegroup) ...
                                        = nan;
                                    frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                        (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                        = nan;
                                    
                                    num_Pairs{i_effect}{i_TD}{i_complexPEeffect} = ...
                                        num_Pairs{i_effect}{i_TD}{i_complexPEeffect} + 1;
                                    num_InvalidPairs{i_effect}{i_TD}{i_complexPEeffect} = ...
                                        num_InvalidPairs{i_effect}{i_TD}{i_complexPEeffect} + 1;
                                end
                            end
                        else %in case of double dipping, NaN results but don't count as invalid
                            for i_tonegroup = 1:size(param.GC.tone_aggregation,2)
                                temporalGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.source2target...
                                    (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                    = nan;
                                temporalGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                temporalGC_pval{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,i_tonegroup) ...
                                    = nan;
                                frequencyGC{i_effect}{i_TD}{i_complexPEeffect}.target2source...
                                    (i_sourceelec,i_targetelec,:,i_tonegroup) ...
                                    = nan;
                            end
                        end
                    end
                end
                disp(['GC computation for ' ToneDur_text{i_TD} 's TD finished - ' ...
                    num2str(num_Pairs{i_effect}{i_TD}{i_complexPEeffect} - ...
                    num_InvalidPairs{i_effect}{i_TD}{i_complexPEeffect}) '/' ...
                    num2str(num_Pairs{i_effect}{i_TD}{i_complexPEeffect}) ...
                    ' valid (GC computation did not bail out)']);
            end
        end
    end
    
    disp(['-- Computing GC from source elecs to target elecs for ' ...
        sub ' and effect ' ...
        param.GC.ElecPairEffect{i_effect} ' (' ...
        param.GC.ElecPairRegion ') took ' ...
        num2str(round(toc(tstart),0)) ' secs --'])
    
end

%% 5) Copy results in common output file
GCdata.ComplexPE.dim                  = 'ElecSelectEffect_TD_ComplexPEEffect';
GCdata.ComplexPE.temporalGC           = temporalGC;
GCdata.ComplexPE.pval_temporalGC      = temporalGC_pval;
GCdata.ComplexPE.spectralGC           = frequencyGC;
GCdata.ComplexPE.trialindices         = i_trialselection_complexPE;

GCdata.ComplexPE.info.validVAR        = validA;
GCdata.ComplexPE.info.validGC         = validF;
GCdata.ComplexPE.info.numPairs        = num_Pairs;
GCdata.ComplexPE.info.numInvalidPairs = num_InvalidPairs;
GCdata.ComplexPE.info.morder          = param.GC.ComplexPE.morder_pereffectTD;

end