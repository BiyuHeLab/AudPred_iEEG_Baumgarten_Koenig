function GCdata = ...
    NASTD_ECoG_Connectivity_CalculateGC_pertone ...
    (sub_list, i_sub, ...
    ToneDur_text, SelElecs, ...
    param, paths_NASTD_ECoG)
%Aim: Compute temporal and spectral GC for all selected electrode pairs in
%a specific subject. Return temporal ans spectral Granger values and p-values.

%Notes:
% Do computation subject-wise, then looped over each electrode pair
%1) Select electrodes from raw data
%2) Determine electrode pairs
%3) perform downsampling (factor 2) before
%4) GC computation
%4.1) determine model order (tsdata_to_infocrit.m)
%4.2) estimate VAR (vector auto regressive) model (tsdata_to_var.m)
%VAR estimates the (linear) relationship of several variables over time;
%Relates current observations of a variable with past observations
%of itself and past observations of other variables.
%Vectors are just used as a convenient form of data storage here
%which sums up the original matrix form of data storgare
%Contains n-equations (for n endogenous variables) and p-lags
%Lag selection is minimized via some model selection criteria
%4.3) check VAR model for problems
%4.4) Compute GC in time domain or frequency domain
%5) GC Statistics


%% 0.1) Set empty GC data struct and output paths
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/core/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/utils/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/stats/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/gc/');

path_fig = ([paths_NASTD_ECoG.ECoGdata_Connectivity 'GC/Figs/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) Load and prepare single-subject data

%Load in preproc single-trial data
sub = sub_list{i_sub};
disp(['Computing CG for subject: ' sub])

NASTD_ECoG_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;

loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
tic
disp(['Loading preprocessed data set for sub: ' sub])
load(loadfile_ECoGpreprocdata, 'DataClean_AllTrials');
disp(['done loading in ' num2str(toc) ' sec'])

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

%1.4 Determine TP/samples for each tone start+end
for i_TD = 1:length(ToneDur_text)
    curr_toneDur_text = ToneDur_text{i_TD};
    
    TP_Tone_StartStop{i_TD} = NaN(36,2);
    TP_Tone_StartStop{i_TD}(1,1) = 0; %p1 set as t = 0 in trial definition
    for i_tone = 2:36
        Dist = abs(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} - ((str2num(curr_toneDur_text)*i_tone) - str2num(curr_toneDur_text)));
        minDist = min(Dist);
        i_minDist = find(Dist == minDist);
        TP_Tone_StartStop{i_TD}(i_tone,1) = DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(i_minDist);
    end
    for i_tone = 1:35
        i_LastSampleTone = find(DataClean_CleanTrialsElecs_perTD{i_TD}.time{1} == TP_Tone_StartStop{i_TD}(i_tone + 1,1));
        TP_Tone_StartStop{i_TD}(i_tone,2) = DataClean_CleanTrialsElecs_perTD{i_TD}.time{1}(i_LastSampleTone);
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
end

%Cleanup
clear DataClean_CleanTrialsElecs DataClean_CleanTrials DataClean_AllTrials

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

%1.6 Remove trials where some electrodes have only NaN samples (NY798)
for i_TD = 1:length(ToneDur_text)
    DataClean_CleanTrialsElecs_perTD{i_TD} = ...
        NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
        (sub, ToneDur_text{i_TD}, DataClean_CleanTrialsElecs_perTD{i_TD});
end

%1.7 Apply zero-phase-shift filter to get input data
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

%1.8 Read out input data parameters
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


%% 2) Select and pair electrodes
%Determine preproc data indices of selected paired electrodes
if strcmp(param.GC.label_ElecPairSel{1}, 'Pred_Pred')
    label_effect1 = 'PredEffect';
    label_effect2 = 'PredEffect';  
    title_label1   = 'Pred2Pred';
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_PE')
    label_effect1 = 'PEEffect';
    label_effect2 = 'PEEffect';  
    title_label1   = 'PE2PE';    
elseif strcmp(param.GC.label_ElecPairSel{1}, 'Pred_PE')
    label_effect1 = 'PredEffect';
    label_effect2 = 'PEEffect';  
    title_label1   = 'Pred2PE';   
elseif strcmp(param.GC.label_ElecPairSel{1}, 'PE_Pred')
    label_effect1 = 'PEEffect';
    label_effect2 = 'PredEffect';      
    title_label1   = 'PE2Pred';   
end
        
if strcmp(param.GC.label_ElecPairSel{2}, 'AllRegions')
    title_label2   = 'AllReg';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Frontal_Temporal')
    title_label2   = 'Front2Temp';
elseif strcmp(param.GC.label_ElecPairSel{2}, 'Temporal_Frontal')
    title_label2   = 'Temp2Front';
end

%Determine elec index as listed in preproc data (for data selection)
elec_pairings(:,1) = SelElecs.(label_effect1){SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1),10};
elec_pairings(:,2) = SelElecs.(label_effect2){SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2),10};
ind_pairedelecs_from    = unique(elec_pairings(:,1), 'stable');
ind_pairedelecs_to      = unique(elec_pairings(:,2), 'stable');
ind_pairedelecs_double  = intersect(ind_pairedelecs_from,ind_pairedelecs_to);
ind_pairedelecs         = [ind_pairedelecs_from; ind_pairedelecs_to];
label_pairedelecs       = label_allelecs(ind_pairedelecs);

%Sanity Check if labels from pairing and preproc data match
if strcmp(param.GC.label_ElecPairSel{1}, 'Pred_Pred') || strcmp(param.GC.label_ElecPairSel{1}, 'PE_PE')
    %Determine elec index as listed in elec pairs data (for label comparison)
    ind2_pairedelecs        = ...
        [unique(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1)); ...
        unique(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2))];
    for i_selelec = 1:length(ind_pairedelecs)
        if strcmp(...
                label_allelecs(ind_pairedelecs(i_selelec)), ...
                SelElecs.(label_effect1){ind2_pairedelecs(i_selelec),1}) == false
            disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs(i_selelec))])
            pause
        end
    end
elseif strcmp(param.GC.label_ElecPairSel{1}, 'Pred_PE') || strcmp(param.GC.label_ElecPairSel{1}, 'PE_Pred')
    %Determine elec index as listed in elec pairs data (for label comparison)
    ind2_pairedelecs1 = unique(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1)); 
    ind2_pairedelecs2 = unique(SelElecs.Pairs.(param.GC.label_ElecPairSel{1}).(param.GC.label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2));    
    for i_selelec = 1:length(ind_pairedelecs_from)
        if strcmp(...
                label_allelecs(ind_pairedelecs_from(i_selelec)), ...
                SelElecs.(label_effect1){ind2_pairedelecs1(i_selelec),1}) == false
            disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_from(i_selelec))])
            pause
        end
    end    
    for i_selelec = 1:length(ind_pairedelecs_to)
        if strcmp(...
                label_allelecs(ind_pairedelecs_to(i_selelec)), ...
                SelElecs.(label_effect2){ind2_pairedelecs2(i_selelec),1}) == false
            disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs_to(i_selelec))])
            pause
        end
    end
end
   
%Provide elec info
label_sourceelecs = [];
for i_elec = 1:length(label_allelecs(ind_pairedelecs_from))
    label_sourceelecs = [label_sourceelecs, char(label_allelecs(ind_pairedelecs_from(i_elec))) ', ' ];
end
label_targetelecs = [];
for i_elec = 1:length(label_allelecs(ind_pairedelecs_to))
    label_targetelecs = [label_targetelecs, char(label_allelecs(ind_pairedelecs_to(i_elec))) ', ' ];
end
disp(['-- Computing GC from source elecs (' label_sourceelecs ') to target elecs (' label_targetelecs ')'])

clear ind2_pairedelecs* elec_pairings


%% 3) Set up empty structs and select elcs and samples for input data 
clear trialData validF validA temporalGC temporalGC_pval frequencyGC

for i_TD = 1:length(ToneDur_text)
    
    %3.1 Bring GC input data into requested format (dim: channels, timepoints, trials;)
    trialData{i_TD} = [];
    for i_trial = 1:numTrials(i_TD)
        trialData{i_TD}(:,:,i_trial) = ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata{i_trial}...
            (:,1:numSamplesperTrial(i_TD));
    end

    %Restrict input data to selected electrodes and time frame
    temp_trialData{i_TD} = trialData{i_TD}(ind_pairedelecs, :, :);
    
%     %Restrict trial data to sequence presentation samples
%     temp2_trialData{i_TD} = ...
%         temp_trialData{i_TD}(:,Sample_Tone_StartStop{i_TD}(1,1):Sample_Tone_StartStop{i_TD}(34,2),:);    

%     %Optional: downsample
%     if param.GC.downsample > 0
%         for i_TD = 1:length(ToneDur_text)
%             trialData{i_TD} = trialData{i_TD}(:,1:param.GC.downsample:end,:);
%         end
%     end

    trialData{i_TD}       = temp_trialData{i_TD};    
    
    %Set structs for Output parameters
    validA{i_TD}          = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1)); %Indicates if VAR coefficients matrix is valid
    validF{i_TD}          = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1)); %Indicates if GC computation is valid
   
    temporalGC{i_TD}.source2target      = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1), size(Sample_Tone_StartStop{i_TD},1)); %pairwise-conditional time-domain multivariate Granger causalities
    temporalGC_pval{i_TD}.source2target = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1), size(Sample_Tone_StartStop{i_TD},1)); %pairwise-conditional frequency-domain/ spectral multivariate Granger causalites
    temporalGC{i_TD}.target2source      = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1), size(Sample_Tone_StartStop{i_TD},1)); %pairwise-conditional time-domain multivariate Granger causalities
    temporalGC_pval{i_TD}.target2source = nan(size(ind_pairedelecs_from,1), size(ind_pairedelecs_to,1), size(Sample_Tone_StartStop{i_TD},1)); %pairwise-conditional frequency-domain/ spectral multivariate Granger causalites
    %dim: source elec * target elec * tone
    
    if param.GC.downsample > 0
        frequencyGC{i_TD}.source2target = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (param.GC.fs/param.GC.downsample)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));
        frequencyGC{i_TD}.target2source = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (param.GC.fs/param.GC.downsample)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));    
    else
        frequencyGC{i_TD}.source2target = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (param.GC.fs)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));        
        frequencyGC{i_TD}.target2source = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (param.GC.fs)+1, ...
            size(Sample_Tone_StartStop{i_TD},1)); 
    end
    %dim: source elec * target elec * frequencies * tone
end

clear temp_trialData

%% 4) Check input data for stationarity
%Optional: check with reverse arrangement test if input data is stationary
% for i_TD = 1:length(ToneDur_text)
%     h = waitbar(0, 'Testing stationarity across elecs, tone, and trials...');
%     
%     output_RAtest{i_TD} = ... %dim: elec pairs, tones, trials
%         nan(size(trialData{i_TD},1), ...
%         size(Sample_Tone_StartStop{i_TD},1), ...
%         size(trialData{i_TD},3));
%     
%     for i_elec = 1:size(trialData{i_TD},1)
%         for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)
%             for i_trial = 1:size(trialData{i_TD},3)
%                 waitbar(i_elec / size(trialData{i_TD},1));
%                 
%                 %Tone-wise data epochs not detrended prior to analysis
%                 %                 output_RAtest{i_TD}(i_elec, i_tone, i_trial) = ...
%                 %                     RA_test(...
%                 %                     trialData{i_TD}(...
%                 %                     i_elec,...
%                 %                     Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),...
%                 %                     i_trial),1);
%                 
%                 %Tone-wise data epochs detrended prior to analysis
%                 temp_inputdata = [];
%                 temp_inputdata = detrend(trialData{i_TD}...
%                     (i_elec,Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),i_trial));
%                 output_RAtest{i_TD}(i_elec, i_tone, i_trial) = ...
%                     RA_test(temp_inputdata,1);
%             end
%         end
%     end
%     close(h)
% 
%     temp_numhits = length(find(output_RAtest{i_TD}));
%     temp_numtotal = size(output_RAtest{i_TD},1)*size(output_RAtest{i_TD},2)*size(output_RAtest{i_TD},3);
%     disp([num2str(temp_numhits) '/' num2str(temp_numtotal) ...
%         ' (' num2str(round((temp_numhits/temp_numtotal)*100)) ...
%         '%) of elec*tone*trial combinations stationary (TD' (num2str(i_TD)) ')']);
% end


%% 5) Compute GC for each tone snippet in each selected electrode pairing
%5.1 Preprocess to increase stationarity: 
for i_TD = 1:length(ToneDur_text)
    for i_elec = 1:size(trialData{i_TD},1)
        for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)
            for i_trial = 1:size(trialData{i_TD},3)
                %1) Detrend tone snippets for each electrode
                trialData{i_TD}...
                    (i_elec,Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),i_trial) = ...
                    detrend(...
                    trialData{i_TD}...
                    (i_elec,Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),i_trial)); 
                
                %2) Subtract mean activity across electrodes per sample                
            end
        end
    end
end


for i_TD = 1:length(ToneDur_text)
    for i_sourceelec = 1:length(ind_pairedelecs_from)
        for i_targetelec = 1:length(ind_pairedelecs_to)            
            %Avoid double-dipping (identical electrode as source + target)
            if ind_pairedelecs_from(i_sourceelec) ~= ind_pairedelecs_to(i_targetelec)
                %4.1) determine appropriate model order for VAR model (tsdata_to_infocrit.m)
                %(balance maximum lag p to achieve best model fit while
                %avoiding overfitting)
                % Calculate standerdized information criteria up to a priori specified maximum model order.
                for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)
                    ptic('\n*** Determine model order (tsdata_to_infocrit)\n');
                    [AIC, BIC, ...
                        moAIC(i_sourceelec, i_targetelec, i_tone), ...
                        moBIC(i_sourceelec, i_targetelec, i_tone)] = ...
                        tsdata_to_infocrit(...
                        trialData{i_TD}...
                        ([i_sourceelec length(ind_pairedelecs_from)+i_targetelec], ...
                        Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2), ...
                        :), ...
                        param.GC.maxmorder, param.GC.regmode);
                    ptoc('*** tsdata_to_infocrit took ');
                    % figure(1); clf;
                    % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/param.GC.fs);
                    % title('Model order estimation');
                end
            else
                for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)
                    moAIC(i_sourceelec, i_targetelec, i_tone) = NaN;
                    moBIC(i_sourceelec, i_targetelec, i_tone) = NaN;
                end                
            end
        end
    end
    %Select median model order across elec pairs * tone snippets as final model order
    %(similar to Hardston et al., Nat Commun 2021)
    param.GC.morder = floor(nanmedian(moBIC, 'all'));
    fprintf('\nbest model order (BIC) = %d\n',param.GC.morder);
    
    %     %Optional: Sample autocovariance estimation
    %     G = tsdata_to_autocov(trialData{i_TD},param.GC.morder)
    
    num_Pairs           = 0;
    num_InvalidPairs    = 0;  
    
    for i_sourceelec = 1:length(ind_pairedelecs_from)
        for i_targetelec = 1:length(ind_pairedelecs_to)
            %Avoid double-dipping (identical electrode as source + target)
            if ind_pairedelecs_from(i_sourceelec) ~= ind_pairedelecs_to(i_targetelec)
                disp(['-- Comuting GC for '...
                    label_pairedelecs{i_sourceelec} ' -> ' ...
                    label_pairedelecs{length(ind_pairedelecs_from) + i_targetelec}]);
                
                for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)                    
                    AS = struct; %Struct for VAR model (contains: A = VAR coefficients matrix, SIG = residuals covariance matrix)
                    FS = struct; %Struct for pairwise-conditional time-domain multivariate Granger causalities
                    
                    %4.2 Fit VAR model to time series data (for up to selected model order)
                    %(obtain estimates of model parameters which maximize the likelihood
                    %function for both full (i.e., X has influence on Y) and reduced
                    %(i.e., X has no influence on Y) VAR models - i.e., model fitting)
                    ptic('\n*** estimating VAR model (tsdata_to_var)... ');
                    [AS.A, AS.SIG] = ...
                        tsdata_to_var...
                        (trialData{i_TD}([i_sourceelec (i_targetelec+length(ind_pairedelecs_from))], ...
                        Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2), :), ...
                        param.GC.morder, param.GC.regmode);
                    
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
                    validA{i_TD}(i_sourceelec, i_targetelec) = ...
                        ~isbad(AS.A); %Record potential errors
                    
                    %4.3 Directly compute GC and corresponding p-value from full
                    %(i.e., X has influence on Y) and reduced/null (i.e., X has no influence on Y)
                    %VAR models using novel State-space version (not degraded by
                    %moving average components in input data)
                    ptic('*** computing time domain GC estimates and testing against reduced model (var_to_pwcgc)... ');
                    [FS.F, FS.pval] = ...
                        var_to_pwcgc(...
                        AS.A, AS.SIG, ...
                        trialData{i_TD}([i_sourceelec (i_targetelec+length(ind_pairedelecs_from))], ...
                        Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2), :), ...
                        param.GC.regmode, param.GC.tstat);
                    %                 F  = Time-domain multivariate Granger causality: for the pairwise-
                    %                 conditional case, an n x n matrix containing pairwise-
                    %                 conditional causalities; the first index then references target
                    %                 ("to") variables and the second source ("from") variables, with
                    %                 NaNs on the diagonal.
                    ptoc;
                    %F = pairwise conditional GC estimates (dim: elec*elec)
                    %pval = p-value from comparison full model vs. reduced/null
                    %model (i.e., does the past of source elec contain information
                    %about the future of target elec, above and beyond information
                    %conatained in the past of the target elec)
                    % Check for failed GC calculation
                    assert(~isbad(FS.F,false),'GC calculation failed - bailing out');
                    validF{i_TD}(i_sourceelec,i_targetelec) = ~isbad(FS.F,false);
                    
                    %4.4 Calculate spectral pairwise-conditional causalities at given frequency
                    % resolution by state-space method.
                    ptic('*** computing frequency domain GC estimates  (var_to_spwcgc)... ');
                    FS.f = var_to_spwcgc...
                        (AS.A, AS.SIG, ...
                        param.GC.newfs);
                    assert(~isbad(FS.f,false),'spectral GC calculation failed - bailing out');
                    ptoc;
                    
                    %4.5 Copy results in output struct
                    if validF{i_TD}(i_sourceelec,i_targetelec) == 1
                        %Direction source - target
                        temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)        = FS.F(1,2);%(Source elec, Target elec, Tone)
                        temporalGC_pval{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)   = FS.pval(1,2);
                        frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,i_tone)     = FS.f(1,2,:);%(Source elec, Target elec, Frequency, Tone)
                        %Direction target - source
                        temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)        = FS.F(2,1);%(Target elec, Source elec, Tone)
                        temporalGC_pval{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)   = FS.pval(2,1);
                        frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,i_tone)     = FS.f(2,1,:);%(Target elec, Source elec, Frequency, Tone)                        
                        num_Pairs = num_Pairs + 1;
                    else
                        temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)         = nan;
                        temporalGC_pval{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)    = nan;
                        frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,i_tone)      = nan;
                        temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)         = nan;
                        temporalGC_pval{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)    = nan;
                        frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,i_tone)      = nan;                        
                        num_Pairs = num_Pairs + 1;
                        num_InvalidPairs = num_InvalidPairs + 1;
                    end
                    
                    %To do: % Check that spectral causalities average (integrate) to time-domain
                    % causalities, as they should according to theory.
                end
            else %in case of double dipping, NaN results but don't count as invalid
                for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)   
                    temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)         = nan;
                    temporalGC_pval{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)    = nan;
                    frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,i_tone)      = nan;
                    temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)         = nan;
                    temporalGC_pval{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)    = nan;
                    frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,i_tone)      = nan;
                end                
            end            
        end
    end    
    disp(['GC computation for ' ToneDur_text{i_TD} 's TD finished - ' ...
        num2str(num_Pairs - num_InvalidPairs) '/' num2str(num_Pairs) ' valid (GC computation bailed out)']);    
end

%Copy results in common output file    
GCdata  = struct;
GCdata.temporalGC           = temporalGC;
GCdata.pval_temporalGC      = temporalGC_pval;
GCdata.spectralGC           = frequencyGC;

GCdata.info.validVAR        = validA;
GCdata.info.validGC         = validF;
GCdata.info.numPairs        = num_Pairs;
GCdata.info.numInvalidPairs = num_InvalidPairs;

GCdata.info.sub             = sub;
GCdata.info.InputDataType   = param.GC.InputDataType;
GCdata.info.ElecPairSel     = param.GC.label_ElecPairSel;
GCdata.info.param           = param.GC;

GCdata.label_allelecs       = label_allelecs;
GCdata.ind_pairedelecs_from = ind_pairedelecs_from;
GCdata.ind_pairedelecs_to   = ind_pairedelecs_to;


% %% 5) plot GC estimate for time and frequency domain
% tone_grouping = [1 11; 12 22; 23 33; 34 34];
% colormap parula
% colormap_parula = colormap;
% colormap gray
% colormap_gray = colormap;
% close
% %5.1 Plot Source to Target GC as function of (aggregated) tones over the course of the sequence
% for i_TD = 1:length(ToneDur_text)   
% 
%     %Summary GC and p-values across pairs per tone group
%     f1 = figure;
%     set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     DimSubplot = [length(tone_grouping) 2];
%     sgtitle(['Time-domain pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
%     subplot_num = 0;
%     
%     %Determine Clim
%     for i_tonegroup = 1:length(tone_grouping)
%         maxC_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(temporalGC{i_TD}.source2target...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3)));
%     end
%     maxC_tempGC{i_TD} = max(maxC_tempGC{i_TD});
%     
%     for i_tonegroup = 1:length(tone_grouping)
%         %GC values
%         subplot_num = subplot_num + 1;
%         
%         f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
%         imagesc(mean(temporalGC{i_TD}.source2target...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
%         f_sub.XTick =  [1:length(ind_pairedelecs_to)];
%         f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
%         f_sub.XLabel.String = {'Target Elecs'};
%         f_sub.YTick =  [1:length(ind_pairedelecs_from)];
%         f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
%         f_sub.YLabel.String = {'Source Elecs'};
%         f_sub.CLim = [0 maxC_tempGC{i_TD}];
%         f_sub.Colormap = colormap_parula;
%         colorbar
%         title(['GC estimates - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
%         %p-values
%         subplot_num = subplot_num + 1;
%         
%         f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
%         imagesc(mean(temporalGC_pval{i_TD}.source2target...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
%         f_sub.XTick =  [1:length(ind_pairedelecs_to)];
%         f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
%         f_sub.XLabel.String = {'Target Elecs'};
%         f_sub.YTick =  [1:length(ind_pairedelecs_from)];
%         f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
%         f_sub.YLabel.String = {'Source Elecs'};
%         f_sub.CLim = [0 0.05];
%         f_sub.Colormap = colormap_gray;
%         colorbar
%         title(['GC p-values (vs. null model, uncorrected) - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
%     end
%     
%     
%     %Temporal GC per connection per tone group
%     DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
%     i_subplot = 0;
%     f1 = figure;
%     set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     %Determine common y-axis limits
%     for i_tonegroup = 1:length(tone_grouping)
%        maxY_tempGC{i_TD}(i_tonegroup,1) = max(max(mean(temporalGC{i_TD}.source2target...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3)));
%     end
%     maxY_tempGC{i_TD} = max(maxY_tempGC{i_TD});
% 
%     for i_sourceelec = 1:length(ind_pairedelecs_from)
%         for i_targetelec = 1:length(ind_pairedelecs_to)
%             i_subplot = i_subplot + 1;
%             f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
%             bar([mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
%                 mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
%                 mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
%                 temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)])
%             hold on
%             errorbar([mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
%                 mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
%                 mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
%                 temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)],...
%                 [std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11)), ...
%                 std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22)), ...
%                 std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33)), 0], ...
%                 'LineStyle', 'none', 'Color', 'k');
% 
%             f_sub.XTickLabel = {'1-11', '12-22', '23-33', '34'};
%             f_sub.XLabel.String = {'Tones'};
%             f_sub.YLabel.String = {'GC estimate'};
%             f_sub.YLim = [0 maxY_tempGC{i_TD}*1.25];
% 
%             title(['From ' label_allelecs{ind_pairedelecs_from(i_sourceelec)} ' to ' label_allelecs{ind_pairedelecs_to(i_targetelec)}])
%         end
%     end
%     sgtitle(['Time-domain pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
%     
%     
%     %Spectral GC per tone group
%     DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
%     i_subplot = 0;
%     f2 = figure;
%     set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     colormap cool
%     colormap_cool = colormap;
%     %Determine plot settings based on input freq
%     if strcmp(param.GC.InputDataType{1}, 'Broadband')
%         temp_xlim = [1 150];
%         temp_xtick = [0:25:150];
%         temp_xticklabel = {'0', '25', '50', '75', '100', '125', '150'};
%     elseif strcmp(param.GC.InputDataType{1}, 'HP05toLP30Hz')
%         temp_xlim = [0 30];
%         temp_xtick = [0:5:30];
%         temp_xticklabel = {'0', '5', '10', '15', '20', '25', '30'};
%     elseif strcmp(param.GC.InputDataType{1}, 'HighGamma_LogAmp')
%         temp_xlim = [70 150];
%         temp_xtick = [70:10:150];
%         temp_xticklabel = {'70', '80', '90', '100', '110', '120', '130', '140', '150'};
%     end
%     
%     %Determine common y-axis limits 
%     for i_tonegroup = 1:length(tone_grouping)
%         maxY_spectGC{i_TD}(i_tonegroup,1) = ...
%             max(max(max(mean(...
%             frequencyGC{i_TD}.source2target...
%             (:,:,temp_xlim(1):temp_xlim(2),tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),4))));
%     end    
%     maxY_spectGC{i_TD} = max(maxY_spectGC{i_TD});
%     
%     for i_sourceelec = 1:length(ind_pairedelecs_from)
%         for i_targetelec = 1:length(ind_pairedelecs_to)
%             i_subplot = i_subplot + 1;
%             f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
%             plot(squeeze(mean(frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,1:11),4)), 'Color', colormap_cool(50,:), 'LineWidth', 1);
%             hold on;
%             plot(squeeze(mean(frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,12:22),4)), 'Color', colormap_cool(100,:), 'LineWidth', 1);
%             plot(squeeze(mean(frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,23:33),4)), 'Color', colormap_cool(150,:), 'LineWidth', 1);
%             plot(squeeze(mean(frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,34),4)), 'Color', colormap_cool(200,:), 'LineWidth', 1);
%             
%             f_sub.XLim = temp_xlim;
%             f_sub.XTick = temp_xtick;
%             f_sub.XTickLabel = temp_xticklabel;
%             
%             f_sub.XLabel.String = {'Frequency [Hz]'};
%             f_sub.YLabel.String = {'GC'};
%             f_sub.YLim = [0 maxY_spectGC{i_TD}];
% 
%             if i_sourceelec ==1 && i_targetelec == 1
%                legend({'Tone 1-11', '12-22', '23-33', '34'}) 
%             end
% 
%             title(['From ' label_allelecs{ind_pairedelecs_from(i_sourceelec)} ' to ' label_allelecs{ind_pairedelecs_to(i_targetelec)}])
%         end
%     end
%     sgtitle(['Spectral pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
%      
%     
% end
% 
% 
% %5.2 Plot reversed (Target to Source) GC as function of (aggregated) tones over the course of the sequence
% for i_TD = 1:length(ToneDur_text)   
% 
%     %Summary GC and p-values across pairs per tone group
%     f1 = figure;
%     set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     DimSubplot = [length(tone_grouping) 2];
%     sgtitle(['Reversed (target-to-source) Time-domain pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - ' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
%     subplot_num = 0;
%     
%     for i_tonegroup = 1:length(tone_grouping)
%         %GC values
%         subplot_num = subplot_num + 1;
%         
%         f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
%         imagesc(mean(temporalGC{i_TD}.target2source...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
%         f_sub.XTick =  [1:length(ind_pairedelecs_to)];
%         f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
%         f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
%         f_sub.YTick =  [1:length(ind_pairedelecs_from)];
%         f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
%         f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
%         f_sub.CLim = [0 maxC_tempGC{i_TD}];
%         f_sub.Colormap = colormap_parula;
%         colorbar
%         title(['GC estimates - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
%         
%         %p-values
%         subplot_num = subplot_num + 1;
%         
%         f_sub = subplot(DimSubplot(1),DimSubplot(2),subplot_num);
%         imagesc(mean(temporalGC_pval{i_TD}.target2source...
%             (:,:,tone_grouping(i_tonegroup,1):tone_grouping(i_tonegroup,2)),3))
%         f_sub.XTick =  [1:length(ind_pairedelecs_to)];
%         f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
%         f_sub.XLabel.String = {'Source (revsered Target) Elecs'};
%         f_sub.YTick =  [1:length(ind_pairedelecs_from)];
%         f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
%         f_sub.YLabel.String = {'Target (reversed Source) Elecs'};
%         f_sub.CLim = [0 0.05];
%         f_sub.Colormap = colormap_gray;
%         colorbar
%         title(['GC p-values (vs. null model, uncorrected) - Tone ' num2str(tone_grouping(i_tonegroup,1)) ' to ' num2str(tone_grouping(i_tonegroup,2))]);
%     end
%      
%     
%     %Temporal GC per connection per tone group
%     DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
%     i_subplot = 0;
%     f1 = figure;
%     set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
% 
%     for i_sourceelec = 1:length(ind_pairedelecs_from)
%         for i_targetelec = 1:length(ind_pairedelecs_to)
%             i_subplot = i_subplot + 1;
%             f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
%             bar([mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11),3), ...
%                 mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22),3),  ...
%                 mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33),3), ...
%                 temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,34)])
%             hold on
%             errorbar([mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11),3), ...
%                 mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22),3),  ...
%                 mean(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33),3), ...
%                 temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,34)],...
%                 [std(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,1:11)), ...
%                 std(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,12:22)), ...
%                 std(temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,23:33)), 0], ...
%                 'LineStyle', 'none', 'Color', 'k');
% 
%             f_sub.XTickLabel = {'1-11', '12-22', '23-33', '34'};
%             f_sub.XLabel.String = {'Tones'};
%             f_sub.YLabel.String = {'GC estimate'};
%             f_sub.YLim = [0 maxY_tempGC{i_TD}];
% 
%             title(['From ' label_allelecs{ind_pairedelecs_to(i_targetelec)} ' to '  label_allelecs{ind_pairedelecs_from(i_sourceelec)}])
%         end
%     end
%     sgtitle(['Reversed Time-domain pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
% 
%     
%     %Spectral GC per tone group
%     DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
%     i_subplot = 0;
%     f2 = figure;
%     set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
%     colormap cool
%     colormap_cool = colormap;
%     for i_sourceelec = 1:length(ind_pairedelecs_from)
%         for i_targetelec = 1:length(ind_pairedelecs_to)
%             i_subplot = i_subplot + 1;
%             f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
%             plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,1:11),4)), 'Color', colormap_cool(50,:), 'LineWidth', 1);
%             hold on;
%             plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,12:22),4)), 'Color', colormap_cool(100,:), 'LineWidth', 1);
%             plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,23:33),4)), 'Color', colormap_cool(150,:), 'LineWidth', 1);
%             plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,34),4)), 'Color', colormap_cool(200,:), 'LineWidth', 1);
% 
%             if strcmp(param.GC.InputDataType{1}, 'Broadband')
%                 f_sub.XLim = [1 150];
%                 f_sub.XTick = [0:25:150];
%                 f_sub.XTickLabel = {'0', '25', '50', '75', '100', '125', '150'};
%             elseif strcmp(param.GC.InputDataType{1}, 'HP05toLP30Hz')
%                 f_sub.XLim = [0 30];
%                 f_sub.XTick = [0:5:30];
%                 f_sub.XTickLabel = {'0', '5', '10', '15', '20', '25', '30'};                
%             elseif strcmp(param.GC.InputDataType{1}, 'HighGamma_LogAmp')
%                 f_sub.XLim = [70 150];
%                 f_sub.XTick = [70:10:150];
%                 f_sub.XTickLabel = {'70', '80', '90', '100', '110', '120', '130', '140', '150'};     
%             end
%             f_sub.XLabel.String = {'Frequency [Hz]'};
%             f_sub.YLabel.String = {'GC'};
%             f_sub.YLim = [0 maxY_spectGC{i_TD}];
% 
%             if i_sourceelec ==1 && i_targetelec == 1
%                legend({'Tone 1-11', '12-22', '23-33', '34'}) 
%             end
% 
%             title(['From ' label_allelecs{ind_pairedelecs_to(i_targetelec)} ' to '  label_allelecs{ind_pairedelecs_from(i_sourceelec)}])
%         end
%     end
%     sgtitle(['Reversed spectral pw GC estimates per tone group - ' ...
%         sub ' - ' ToneDur_text{i_TD} 's TD - Reversed' param.GC.label_ElecPairSel{1} ', ' param.GC.label_ElecPairSel{2} ' - ' param.GC.InputDataType{1}], 'Interpreter', 'none')
%        
% end

end