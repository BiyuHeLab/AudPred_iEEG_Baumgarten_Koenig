function [temporalGC, temporalGC_pval, frequencyGC, validA, validF] = ...
    NASTD_ECoG_Connectivity_CalculateGC ...
    (sub_list, i_sub, ...
    ToneDur_text, InputDataType, SelElecs, label_ElecPairSel, ...
    GCparam)
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


%% 0.1) Set empty GC data struct
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/core/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/utils/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/stats/');
addpath('/isilon/LFMI/VMdrive/Thomas/toolboxes/mvgc_v1.3/gc/');


%% 1) Load and prepare single-subject data

%Load in preproc single-trial data
sub = sub_list{i_sub};
disp(['Computing CG for subject: ' sub])

NASTD_ECoG_subjectinfo %load subject info file (var: si)
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
if strcmp(InputDataType, 'Broadband') %1) broadband data
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP150Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP150Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    end
elseif strcmp(InputDataType, 'HP05toLP30Hz') %2) slow freqs (0.5-30 Hz)
    for i_TD = 1:length(ToneDur_text)
        [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.HP05toLP30Hz] = ...
            NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
            (sub, ToneDur_text{i_TD}, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, ...
            plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);
    end    
elseif strcmp(InputDataType, 'HighGamma_LogAmp') %3) high gamma-band amplitude (70-150Hz) 
    FrequencyBands          = [70, 150]; %Hz [HP,LP]
    FrequencyBand_Labels    = {'HighGamma'};    
    for i_TD = 1:length(ToneDur_text)
        for i_freqbands = 1:length(FrequencyBands)            
            FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];            
            [DataClean_CleanTrialsElecs_perTD{i_TD}.filteredtrialdata, ...
                DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, ToneDur_text{i_TD}, ...
                FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel,...
                DataClean_CleanTrialsElecs_perTD{i_TD}, ...
                plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG);            
        end
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
if strcmp(label_ElecPairSel{1}, 'Pred_Pred')
    label_effect1 = 'PredEffect';
    label_effect2 = 'PredEffect';    
elseif strcmp(label_ElecPairSel{1}, 'PE_PE')
    label_effect1 = 'PEEffect';
    label_effect2 = 'PEEffect';    
elseif strcmp(label_ElecPairSel{1}, 'Pred_PE')
    label_effect1 = 'PredEffect';
    label_effect2 = 'PEEffect';  
elseif strcmp(label_ElecPairSel{1}, 'PE_Pred')
    label_effect1 = 'PEEffect';
    label_effect2 = 'PredEffect';      
end
        
%Determine elec index as listed in preproc data (for data selection)
elec_pairings(:,1) = SelElecs.(label_effect1){SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1),10};
elec_pairings(:,2) = SelElecs.(label_effect2){SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2),10};
ind_pairedelecs_from    = unique(elec_pairings(:,1), 'stable');
ind_pairedelecs_to      = unique(elec_pairings(:,2), 'stable');
ind_pairedelecs         = [ind_pairedelecs_from; ind_pairedelecs_to];
label_pairedelecs       = label_allelecs(ind_pairedelecs);

%Sanity Check if labels from pairing and preproc data match
if strcmp(label_ElecPairSel{1}, 'Pred_Pred') || strcmp(label_ElecPairSel{1}, 'PE_PE')
    %Determine elec index as listed in elec pairs data (for label comparison)
    ind2_pairedelecs        = ...
        [unique(SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1)); ...
        unique(SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2))];
    for i_selelec = 1:length(ind_pairedelecs)
        if strcmp(...
                label_allelecs(ind_pairedelecs(i_selelec)), ...
                SelElecs.(label_effect1){ind2_pairedelecs(i_selelec),1}) == false
            disp(['Electrode label mismatch for elec #' num2str(ind_pairedelecs(i_selelec))])
            pause
        end
    end
elseif strcmp(label_ElecPairSel{1}, 'Pred_PE') || strcmp(label_ElecPairSel{1}, 'PE_Pred')
    %Determine elec index as listed in elec pairs data (for label comparison)
    ind2_pairedelecs1 = unique(SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,1)); 
    ind2_pairedelecs2 = unique(SelElecs.Pairs.(label_ElecPairSel{1}).(label_ElecPairSel{2}).Ind_UniquePairings_persub{i_sub}(:,2));    
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
%     if GCparam.downsample > 0
%         for i_TD = 1:length(ToneDur_text)
%             trialData{i_TD} = trialData{i_TD}(:,1:GCparam.downsample:end,:);
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
    
    if GCparam.downsample > 0
        frequencyGC{i_TD}.source2target = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (GCparam.fs/GCparam.downsample)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));
        frequencyGC{i_TD}.target2source = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (GCparam.fs/GCparam.downsample)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));    
    else
        frequencyGC{i_TD}.source2target = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (GCparam.fs)+1, ...
            size(Sample_Tone_StartStop{i_TD},1));        
        frequencyGC{i_TD}.target2source = ...
            nan(size(ind_pairedelecs_from,1), ...
            size(ind_pairedelecs_to,1), ...
            (GCparam.fs)+1, ...
            size(Sample_Tone_StartStop{i_TD},1)); 
    end
    %dim: source elec * target elec * frequencies * tone
end

clear temp_trialData

%% 4) Check input data for stationarity
%Optional: check with reverse arrangement test if input data is stationary
for i_TD = 1:length(ToneDur_text)
    h = waitbar(0, 'Testing stationarity across elecs, tone, and trials...');
    
    output_RAtest{i_TD} = ... %dim: elec pairs, tones, trials
        nan(size(trialData{i_TD},1), ...
        size(Sample_Tone_StartStop{i_TD},1), ...
        size(trialData{i_TD},3));
    
    for i_elec = 1:size(trialData{i_TD},1)
        for i_tone = 1:size(Sample_Tone_StartStop{i_TD},1)
            for i_trial = 1:size(trialData{i_TD},3)
                waitbar(i_elec / size(trialData{i_TD},1));
                
                %Tone-wise data epochs not detrended prior to analysis
                %                 output_RAtest{i_TD}(i_elec, i_tone, i_trial) = ...
                %                     RA_test(...
                %                     trialData{i_TD}(...
                %                     i_elec,...
                %                     Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),...
                %                     i_trial),1);
                
                %Tone-wise data epochs detrended prior to analysis
                temp_inputdata = [];
                temp_inputdata = detrend(trialData{i_TD}...
                    (i_elec,Sample_Tone_StartStop{i_TD}(i_tone,1):Sample_Tone_StartStop{i_TD}(i_tone,2),i_trial));
                output_RAtest{i_TD}(i_elec, i_tone, i_trial) = ...
                    RA_test(temp_inputdata,1);
            end
        end
    end
    close(h)

    temp_numhits = length(find(output_RAtest{i_TD}));
    temp_numtotal = size(output_RAtest{i_TD},1)*size(output_RAtest{i_TD},2)*size(output_RAtest{i_TD},3);
    disp([num2str(temp_numhits) '/' num2str(temp_numtotal) ...
        ' (' num2str(round((temp_numhits/temp_numtotal)*100)) ...
        '%) of elec*tone*trial combinations stationary (TD' (num2str(i_TD)) ')']);
end


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
                    GCparam.maxmorder, GCparam.regmode);
                ptoc('*** tsdata_to_infocrit took ');
                % figure(1); clf;
                % plot_tsdata([AIC BIC]',{'AIC','BIC'},1/GCparam.fs);
                % title('Model order estimation');
            end
        end
    end
    %Select median model order across elec pairs * tone snippets as final model order
    %(similar to Hardston et al., Nat Commun 2021)
    GCparam.morder = floor(median(moBIC, 'all'));
    fprintf('\nbest model order (BIC) = %d\n',GCparam.morder);
    
    %     %Optional: Sample autocovariance estimation
    %     G = tsdata_to_autocov(trialData{i_TD},GCparam.morder)
    
    for i_sourceelec = 1:length(ind_pairedelecs_from)
        for i_targetelec = 1:length(ind_pairedelecs_to)
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
                    GCparam.morder, GCparam.regmode);
                assert(~isbad(AS.A),'VAR estimation failed - bailing out');
                ptoc;
                %A = VAR coefficients matrix (dim: elec * elec * lag/model order)
                %SIG = residuals covariance matrix (dim: elec * elec)
                
                %Check if VAR model estimation worked out (check for rank
                %defficient VAR output matrix due to colinearlity in input data
                %and that VAR spectral radius < 1, indicating stationarity or
                %long-term memory of input data)
                info = var_info(AS.A, AS.SIG);
                assert(~info.error,'VAR error(s) found - bailing out');
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
                    GCparam.regmode, GCparam.tstat);
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
                else
                    temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)         = nan;
                    temporalGC_pval{i_TD}.source2target(i_sourceelec,i_targetelec,i_tone)    = nan;
                    frequencyGC{i_TD}.source2target(i_sourceelec,i_targetelec,:,i_tone)      = nan;
                    temporalGC{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)         = nan;
                    temporalGC_pval{i_TD}.target2source(i_sourceelec,i_targetelec,i_tone)    = nan;
                    frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,i_tone)      = nan;
                end
                
                %To do: % Check that spectral causalities average (integrate) to time-domain
                % causalities, as they should according to theory.
                
                
                %%Classic stepwise version from mvgc1.0
                %             %4.3 Compute autocovariance sequence for VAR model
                %             [GS(i_sourceelec,i_targetelec).G, GS(i_sourceelec,i_targetelec).info] = ...
                %                 var_to_autocov(AS(i_sourceelec, i_targetelec).A, AS(i_sourceelec, i_targetelec).SIG, 1000);
                %
                %             %4.4 Compute pairwise-conditional time-domain multivariate Granger
                %             %causalities from autocovariance sequence
                %             FS(i_sourceelec,i_targetelec).F = autocov_to_pwcgc(GS(i_sourceelec,i_targetelec).G);
                %             %Check if mvgc computation worked out
                %             validF{i_TD}(i_sourceelec,i_targetelec) = ~isbad(FS(i_sourceelec,i_targetelec).F,false);
                %
                %             F = FS(i_sourceelec,i_targetelec).F;
                %
                %             %If MVGC computation for electrode pair is valid, compute p-value
                %             %and store output data
                %             if validF{i_TD}(i_sourceelec,i_targetelec) == 1
                %                 %Test MVGC against null distribution
                %                 pval = mvgc_pval...
                %                     (F, GCparam.morder, numSamplesperTrial(i_TD), numTrials(i_TD), ...
                %                     1, 1, GCparam.nvars-2, GCparam.tstat); % take careful note of arguments!
                %
                %                 %Store p-values and GC-values for current electrode pair
                %                 temporalGC_pval{i_TD}(i_sourceelec,i_targetelec) = pval(1,2);     %(Source -> Target)
                %                 temporalGC_pval{i_TD}(i_targetelec,i_sourceelec) = pval(2,1);     %(Target -> Source)
                %                 temporalGC{i_TD}(i_sourceelec,i_targetelec) = F(1,2);      %(Source -> Target)
                %                 temporalGC{i_TD}(i_targetelec,i_sourceelec) = F(2,1);      %(Target -> Source)
                %
                %                 %4.5 Calculate pairwise-conditional frequency-domain MVGCs from autocovariance sequence
                %                 f = autocov_to_spwcgc(GS(i_sourceelec,i_targetelec).G,GCparam.fs/GCparam.downsample);
                %                 frequencyGC{i_TD}(i_sourceelec,i_targetelec,:) = f(1,2,:);   %(Source -> Target, Frequency)
                %                 frequencyGC{i_TD}(i_targetelec,i_sourceelec,:) = f(2,1,:);   %(Target -> Source, Frequency)
                %             else
                %                 temporalGC_pval{i_TD}(i_sourceelec,i_targetelec) = nan;
                %                 temporalGC_pval{i_TD}(i_targetelec,i_sourceelec) = nan;
                %                 temporalGC{i_TD}(i_sourceelec,i_targetelec) = nan;
                %                 temporalGC{i_TD}(i_targetelec,i_sourceelec) = nan;
                %                 frequencyGC{i_TD}(i_sourceelec,i_targetelec,:) = nan;
                %                 frequencyGC{i_TD}(i_targetelec,i_sourceelec,:) = nan;
                %             end
            end
        end
    end
end

%% 5) plot GC estimate for time and frequency domain

%plot GC as function of (aggregated) tones over the course of the sequence
i_TD = 1;
DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
i_subplot = 0;
f1 = figure;
set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
for i_sourceelec = 1:length(ind_pairedelecs_from)
    for i_targetelec = 1:length(ind_pairedelecs_to)
        i_subplot = i_subplot + 1;
        f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
        bar([mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
            mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
            mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
            temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)])
        hold on
        errorbar([mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11),3), ...
            mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22),3),  ...
            mean(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33),3), ...
            temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,34)],...
            [std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,1:11)), ...
            std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,12:22)), ...
            std(temporalGC{i_TD}.source2target(i_sourceelec,i_targetelec,23:33)), 0], ...
            'LineStyle', 'none');
        
        f_sub.XTickLabel = {'1-11', '12-22', '23-33', '34'};
        f_sub.XLabel.String = {'Tones'};
        f_sub.YLabel.String = {'GC estimate'};
                
        title([label_allelecs{ind_pairedelecs_from(i_sourceelec)} ' -> ' label_allelecs{ind_pairedelecs_to(i_targetelec)}])
    end
end
sgtitle(['Time-domain GC estimates - ' ...
    sub ' - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2}], 'Interpreter', 'none')

%plot spectral GC as function of (aggregated) tones over the course of the sequence
i_TD = 1;
DimSubplot = [length(ind_pairedelecs_from) length(ind_pairedelecs_to)];
i_subplot = 0;
f2 = figure;
set(f2,'units','normalized','outerposition',[0 0 1 1]) %full screen
colormap cool
colormap_cool = colormap;
for i_sourceelec = 1:length(ind_pairedelecs_from)
    for i_targetelec = 1:length(ind_pairedelecs_to)
        i_subplot = i_subplot + 1;
        f_sub = subplot(DimSubplot(1),DimSubplot(2),i_subplot);
        plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,1:11),4)), 'Color', colormap_cool(50,:));
        hold on;
        plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,12:22),4)), 'Color', colormap_cool(100,:));
        plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,23:33),4)), 'Color', colormap_cool(150,:));
        plot(squeeze(mean(frequencyGC{i_TD}.target2source(i_sourceelec,i_targetelec,:,34),4)), 'Color', colormap_cool(200,:));
        
        f_sub.XLim = [1 150];
        f_sub.XTick = [0:25:150];
        f_sub.XTickLabel = {'0', '25', '50', '75', '100', '125', '150'};
        f_sub.XLabel.String = {'Frequency [Hz]'};
        f_sub.YLabel.String = {'GC'};
        
        if i_sourceelec ==1 && i_targetelec == 1
           legend({'1-11', '12-22', '23-33', '34'}) 
        end
                
        title([label_allelecs{ind_pairedelecs_from(i_sourceelec)} ' -> ' label_allelecs{ind_pairedelecs_to(i_targetelec)}])
    end
end
sgtitle(['Spectral GC estimates - ' ...
    sub ' - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2}], 'Interpreter', 'none')




f1 = figure;
set(f1,'units','normalized','outerposition',[0 0 1 1]) %full screen
DimSubplot = [2 3];
sgtitle(['Pairwise-conditional GC estimates - ' ...
    sub ' - ' label_ElecPairSel{1} ', ' label_ElecPairSel{2}], 'Interpreter', 'none')

colormap parula
colormap_parula = colormap;
colormap gray
colormap_gray = colormap;
  
%TD 1 
%time domain GC values
f_sub = subplot(DimSubplot(1),DimSubplot(2),1);
imagesc(mean(temporalGC{1}.source2target(:,:,1:11),3))
f_sub.XTick =  [1:length(ind_pairedelecs_to)];
f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
f_sub.XLabel.String = {'Target Elecs'};
f_sub.YTick =  [1:length(ind_pairedelecs_from)];
f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
f_sub.YLabel.String = {'Source Elecs'};
f_sub.CLim = [0 0.05];
f_sub.Colormap = colormap_parula;
colorbar
title(['Time domain pairwise-conditional GC (' label_effect1 '->' label_effect2 ') - TD 1']);  
%p-values
f_sub = subplot(DimSubplot(1),DimSubplot(2),2);
imagesc(mean(temporalGC_pval{1}.source2target(:,:,1:11),3))
f_sub.XTick =  [1:length(ind_pairedelecs_to)];
f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
f_sub.XLabel.String = {'Target Elecs'};
f_sub.YTick =  [1:length(ind_pairedelecs_from)];
f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
f_sub.YLabel.String = {'Source Elecs'};
f_sub.CLim = [0 0.05];
f_sub.Colormap = colormap_gray;
colorbar
title(['p-value (' label_effect1 '->' label_effect2 ') - TD 1']); 
%frequency domain GC values
%PROBLEM _ style doesn't work for multiple target electrodes
f_sub = subplot(DimSubplot(1),DimSubplot(2),3);
for i_elec = 1:length(ind_pairedelecs_from)
    hold on;
    plot(squeeze(mean(frequencyGC{1}.source2target(i_elec,:,:,1:11),4))')
end
f_sub.XLim = [1 150];
f_sub.XTick = [0:25:150];
f_sub.XTickLabel = {'0', '25', '50', '75', '100', '125', '150'};
f_sub.XLabel.String = {'Frequency [Hz]'};
f_sub.YLabel.String = {'GC'};
title(['Frequency domain pairwise-conditional GC (' label_effect1 '->' label_effect2 ') - TD 1']);
legend(label_allelecs(ind_pairedelecs_from))

%TD 2 GC values
%GC values
f_sub = subplot(DimSubplot(1),DimSubplot(2),4);
imagesc(temporalGC{2}(:,:,1))
f_sub.XTick =  [1:length(ind_pairedelecs_to)];
f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
f_sub.XLabel.String = {'Target Elecs'};
f_sub.YTick =  [1:length(ind_pairedelecs_from)];
f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
f_sub.YLabel.String = {'Source Elecs'};
f_sub.CLim = [0 0.001];
f_sub.Colormap = colormap_parula;
colorbar
title(['Time domain pairwise-conditional GC (' label_effect1 '->' label_effect2 ') - TD 2']);  
%p-values
f_sub = subplot(DimSubplot(1),DimSubplot(2),5);
imagesc(temporalGC_pval{2}(:,:,1))
f_sub.XTick =  [1:length(ind_pairedelecs_to)];
f_sub.XTickLabel = label_allelecs(ind_pairedelecs_to);
f_sub.XLabel.String = {'Target Elecs'};
f_sub.YTick =  [1:length(ind_pairedelecs_from)];
f_sub.YTickLabel = label_allelecs(ind_pairedelecs_from);
f_sub.YLabel.String = {'Source Elecs'};
f_sub.CLim = [0 0.05];
f_sub.Colormap = colormap_gray;
colorbar
title(['p-value (' label_effect1 '->' label_effect2 ') - TD 2']); 
%frequency domain GC values
f_sub = subplot(DimSubplot(1),DimSubplot(2),6);
for i_elec = 1:length(ind_pairedelecs_from)
    hold on;
    plot(squeeze(frequencyGC{2}(i_elec,:,:,1)))
end
f_sub.XLim = [1 150];
f_sub.XTick = [0:25:150];
f_sub.XTickLabel = {'0', '25', '50', '75', '100', '125', '150'};
f_sub.XLabel.String = {'Frequency [Hz]'};
f_sub.YLabel.String = {'GC'};
title(['Frequency domain pairwise-conditional GC (' label_effect1 '->' label_effect2 ') - TD 2']);
legend(label_allelecs(ind_pairedelecs_from))

end