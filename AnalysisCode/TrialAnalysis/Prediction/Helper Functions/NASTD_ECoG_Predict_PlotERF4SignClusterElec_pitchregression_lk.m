function NASTD_ECoG_Predict_PlotERF4SignClusterElec_pitchregression_lk...
    (subs, FuncInput_EffectType, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    save_poststepFigs, paths_NASTD_ECoG)

%Aim: Plot a summmary figure showing sign. electrodes on surface brain for
%all subjects. For each of these electrodes, show also the traces per p*34
%conditionduring p1-p33 and during p33. Group ERFs by anatomical location.

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add base dir and own script dir
% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); %path to freesurfer where read_surf function is

if param.FDRcorrect == 1
    FDR_label = 'FDRcorr';
else
    FDR_label = 'uncorr';
end

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction ...
    '/PredEffects/Allsub_n' num2str(length(subs)) ...
    '/Figs/PredEffects_ERF/' FuncInput_EffectType ...
    '/ClusterCorr/' param.ElecSelect '/' FDR_label '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end


%% 1) For all subjects, prepare data underlying prediction effect analysis
% Read out preprocessed data for current effect, signal component, TD
tic
for i_sub = 1:length(subs)
    
    %1.1 Load in preprocessed neural and behavioral data
    sub = subs{i_sub};
    disp(['-- Preparing pre-processed data for subject: ' sub ' --'])
    
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    %Load in preprocessed neural and behavioral data
    loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
    tic
    load(loadfile_ECoGpreprocdata);
    
    for i_TD = 1:length(FuncInput_ToneDur_text)
        %1.2 Select relevant electrodes and trials
        temp_DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);
        temp_Index_Valid_Elecs = setdiff(temp_DataClean_CleanTrials.cfg.info_elec.selected.index4EDF, ...
            subs_PreProcSettings.(sub).rejectedChan_index); %via index
        temp_Label_Valid_Elecs = setdiff(temp_DataClean_CleanTrials.cfg.info_elec.selected.Label, ...
            subs_PreProcSettings.(sub).rejectedChan_label); %via label
        if ~isempty(find((strcmp(...
                sort(temp_DataClean_CleanTrials.label(temp_Index_Valid_Elecs)), ...
                sort(temp_Label_Valid_Elecs))) == 0))
            disp(['Mismatch in electrode selection']);
        end
        
        cfg         = [];
        cfg.channel = (temp_Label_Valid_Elecs);
        temp_DataClean_CleanTrialsElecs = ft_selectdata(cfg, temp_DataClean_CleanTrials);
        
        %Load prediction data and select current effect
        temp_DataClean_CleanTrialsElecs = ...
            NASTD_ECoG_Preproc_SelTrialsperTD...
            (FuncInput_ToneDur_text{i_TD}, temp_DataClean_CleanTrialsElecs);
        
        disp(['Total trial number for TD ' FuncInput_ToneDur_text{i_TD} 's: ' ...
            num2str(length(temp_DataClean_CleanTrialsElecs.trial))])
        
        %Adjust info subfield
        temp_DataClean_CleanTrialsElecs.info.elec = ...
            temp_DataClean_CleanTrialsElecs.cfg.previous.info_elec;
        temp_DataClean_CleanTrialsElecs.info.trigger = ...
            temp_DataClean_CleanTrialsElecs.cfg.previous.info_trigger;
        temp_DataClean_CleanTrialsElecs.info.ref = ...
            temp_DataClean_CleanTrialsElecs.cfg.previous.info_ref;
        
        temp_DataClean_CleanTrialsElecs.cfg.previous = ...
            rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_elec');
        temp_DataClean_CleanTrialsElecs.cfg.previous = ...
            rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_trigger');
        temp_DataClean_CleanTrialsElecs.cfg.previous = ...
            rmfield(temp_DataClean_CleanTrialsElecs.cfg.previous, 'info_ref');
        
        %1.3 Baseline Correction: Normalize whole-trial MEG activity relative to prestimulus (i.e., tone sequence start) baseline
        temp_BL_win = [-0.5 0];
        for i_trial = 1:length(temp_DataClean_CleanTrialsElecs.trial)
            for i_elec = 1:length(temp_DataClean_CleanTrialsElecs.label)
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,:) = ...
                    temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,:) - ...
                    nanmean(temp_DataClean_CleanTrialsElecs.trial{i_trial}(i_elec,...
                    find(temp_DataClean_CleanTrialsElecs.time{1} == temp_BL_win(1)): ...
                    find(temp_DataClean_CleanTrialsElecs.time{1} == temp_BL_win(2))));
            end
        end
        
        %1.4 Remove trials where some electrodes have only NaN samples (NY798)
        temp_DataClean_CleanTrialsElecs = ...
            NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials...
            (sub, FuncInput_ToneDur_text{i_TD}, temp_DataClean_CleanTrialsElecs);
        
        %1.5 Read out current signal component
        if strcmp(FuncInput_DataType, 'LP35Hz')
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.LP35Hz] = ...
                NASTD_ECoG_FiltNaNinterp_LP35Hz...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'HP01toLP30Hz')
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP01toLP30Hz] = ...
                NASTD_ECoG_FiltNaNinterp_HP01toLP30Hz...
                (sub, FuncInput_ToneDur_text, ...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
        
        elseif strcmp(FuncInput_DataType, 'HP05toLP30Hz')
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP05toLP30Hz] = ...
                NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
                (sub, FuncInput_ToneDur_text, ...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'HP1toLP30Hz')
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP1toLP30Hz] = ...
                NASTD_ECoG_FiltNaNinterp_HP1toLP30Hz...
                (sub, FuncInput_ToneDur_text, ...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'HP2toLP30Hz')
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.HP2toLP30Hz] = ...
                NASTD_ECoG_FiltNaNinterp_HP2toLP30Hz...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'Alpha_LogAmp')
            temp_FrequencyBands = [8 12];
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.Alpha_LogAmp] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'Beta_LogAmp')
            temp_FrequencyBands = [15 30];
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.Beta_LogAmp] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'Gamma_LogAmp')
            temp_FrequencyBands = [30 70];
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.Gamma_LogAmp] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        elseif strcmp(FuncInput_DataType, 'HighGamma_LogAmp')
            temp_FrequencyBands = [70 150];
            [temp_DataClean_CleanTrialsElecs.trial, ...
                temp_DataClean_CleanTrialsElecs.info.FiltInterp.HighGamma_LogAmp] = ...
                NASTD_ECoG_FiltNaNinterp_AmpEnvel...
                (sub, FuncInput_ToneDur_text{i_TD}, ...
                temp_FrequencyBands(1), temp_FrequencyBands(2), FuncInput_DataType,...
                temp_DataClean_CleanTrialsElecs, ...
                0, 0, paths_NASTD_ECoG);
            
        end
        
        %1.6 Extract p1-33, p33, and p34 for each trial
        %Define parameters depending on TD
        temp_nTrials             = length(temp_DataClean_CleanTrialsElecs.trial);
        temp_nSensors            = size(temp_DataClean_CleanTrialsElecs.trial{1},1);
        temp_ToneDur_Sec         = str2num(FuncInput_ToneDur_text{i_TD});
        
        %Determine TP/samples for each tone start+end
        StimTiming.TP_Tone_StartStop{i_sub, i_TD} = NaN(36,2);
        StimTiming.TP_Tone_StartStop{i_sub, i_TD}(1,1) = 0; %p1 set as t = 0 in trial definition
        for i_tone = 2:36
            temp_Dist = abs(temp_DataClean_CleanTrialsElecs.time{1} - ...
                ((str2num(FuncInput_ToneDur_text{i_TD})*i_tone) - str2num(FuncInput_ToneDur_text{i_TD})));
            temp_minDist = min(temp_Dist);
            i_minDist = find(temp_Dist == temp_minDist);
            StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,1) = temp_DataClean_CleanTrialsElecs.time{1}(i_minDist);
        end
        for i_tone = 1:35
            i_LastSampleTone = find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone + 1,1));
            StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,2) = temp_DataClean_CleanTrialsElecs.time{1}(i_LastSampleTone);
        end
        StimTiming.TP_Tone_StartStop{i_sub, i_TD} = StimTiming.TP_Tone_StartStop{i_sub, i_TD}(1:34,:);
        
        %Check if all tones are of equal length
        %(if not then choose min length by deleting the last sample of longer trials)
        temp_minSeqLength_sec = min(StimTiming.TP_Tone_StartStop{i_sub, i_TD}(:,2)-StimTiming.TP_Tone_StartStop{i_sub, i_TD}(:,1));
        for i_tone = 1:34
            if (StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,2) - StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,1)) > temp_minSeqLength_sec
                StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,2) = ...
                    temp_DataClean_CleanTrialsElecs.time{1}(...
                    find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,2))-1);
            end
        end
        
        %Determine samples corresponding to TP
        StimTiming.Sample_Tone_StartStop{i_sub, i_TD} = NaN(34,2);
        for i_tone = 1:34
            StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(i_tone,1) = ...
                find(temp_DataClean_CleanTrialsElecs.time{1} == StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,1));
            StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(i_tone,2) = ...
                find(StimTiming.TP_Tone_StartStop{i_sub, i_TD}(i_tone,2) == temp_DataClean_CleanTrialsElecs.time{1});
        end
        
        %Initialize and fill 3D data arrays (dim: nSens, nSamples per tone,
        %nTrials) and copy electrode lebels so that we find electrode indices
        ERFdata.p33{i_sub, i_TD} = zeros(temp_nSensors, ...
            length(StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(1,1):StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(1,2)), ...
            temp_nTrials);
        ERFdata.p34{i_sub, i_TD} = ERFdata.p33{i_sub, i_TD};
        ERFdata.p32{i_sub, i_TD} = ERFdata.p33{i_sub, i_TD};
        ERFdata.p31{i_sub, i_TD} = ERFdata.p33{i_sub, i_TD};
        ERFdata.elec_labels{i_sub, i_TD} = temp_DataClean_CleanTrialsElecs.label;
        
        %copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
        for i_trial = 1:temp_nTrials
            ERFdata.p33{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex,2));
            ERFdata.p34{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex+1,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex+1,2));
            ERFdata.p32{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-1,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-1,2));
            ERFdata.p31{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-2,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-2,2));
            ERFdata.p1p33{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(1,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex,2));
        end
        
        %1.7 Separate trials by Predp34
        label_Predp34 = [-1 0 1];% low, medium, high;
        for i_Predp34 = 1:length(label_Predp34)
            temp_TrialFilt_Predp34 = temp_DataClean_CleanTrialsElecs.behav.stim.predID == label_Predp34(i_Predp34);
            ERFdata.Index_Predp34{i_sub, i_TD}{i_Predp34} = find(temp_TrialFilt_Predp34 == 1);
            %Determine trial number in case we want to restrict analysis to equal
            %amount of trials across predp34conditions
            ERFdata.NumTrialsperpredp34(i_sub, i_TD, i_Predp34) = length(ERFdata.Index_Predp34{i_sub, i_TD}{i_Predp34});
        end
        ERFdata.MinNumTrialsperpredp34(i_sub, i_TD) = ...
            min(ERFdata.NumTrialsperpredp34(i_sub, :));
        
        %% 1.8 Separate trials by p32 (frequency at tone 32)
        % Read out p32 frequency for each trial
        temp_p32_freq = nan(temp_nTrials,1);
        
        for i_trial = 1:temp_nTrials
            temp_p32_freq(i_trial) = ...
                temp_DataClean_CleanTrialsElecs.behav.stim.series_f{i_trial}(32);
        end
        
        % Identify unique frequencies (conditions)
        [p32_unique_freqs, ~, freq_idx] = unique(temp_p32_freq);
        
        % Store labels (useful later for plotting)
        ERFdata.p32_freq_labels{i_sub, i_TD} = p32_unique_freqs;
        
        % Initialize indexing structure
        for i_freq = 1:length(p32_unique_freqs)
            
            % Find trials belonging to this frequency
            ERFdata.Index_p32{i_sub, i_TD}{i_freq} = find(freq_idx == i_freq);
            
            % Store number of trials per frequency
            ERFdata.NumTrialsper_p32(i_sub, i_TD, i_freq) = ...
                length(ERFdata.Index_p32{i_sub, i_TD}{i_freq});
        end
        
        %% 1.9 Make p32-evoked
        nElec = size(ERFdata.p32{i_sub,i_TD},1);
        nSamples = size(ERFdata.p32{i_sub,i_TD},2);
        nTrials = size(ERFdata.p32{i_sub,i_TD},3);

        % Initialize p32evoked
        ERFdata.p32evoked{i_sub,i_TD} = zeros(nElec, nSamples, nTrials);
        baseline_ms = 50;
        
        % Loop over trials
        SampleFreq      = DataClean_AllTrials.fsample;
        for i_trial = 1:nTrials
            % Loop over electrodes
            for i_elec = 1:nElec
                % Extract p32 data
                p32_data = squeeze(ERFdata.p32{i_sub,i_TD}(i_elec,:,i_trial));

                % Extract p31 data (tone before p32)
                p31_data = squeeze(ERFdata.p31{i_sub,i_TD}(i_elec,:,i_trial));

                % Determine number of samples corresponding to last 50 ms
                nSamplesBL = round(baseline_ms/1000 * SampleFreq);
                baseline_vals = p31_data(end-nSamplesBL+1:end);

                % Compute mean baseline
                baseline_mean = mean(baseline_vals);

                % Subtract baseline from p32
                ERFdata.p32evoked{i_sub,i_TD}(i_elec,:,i_trial) = p32_data - baseline_mean;
            end
        end
        
        %% 1.10 Sensory tracking (pitch → neural regression in 50 ms windows)
        % --- Parameters ---
%         win_ms = 50;
%         win_samples = round(win_ms/1000 * SampleFreq);
%         
%         data = ERFdata.p1p33{i_sub,i_TD}; % (elec × time × trial)
%         
%         nElec   = size(data,1);
%         nTrials = size(data,3);
%         
%         nTones = 32;
%         
%         % --- Extract pitch for all tones ---
%         pitch_alltones = nan(nTrials, nTones);
%         
%         for i_trial = 1:nTrials
%             seq = temp_DataClean_CleanTrialsElecs.behav.stim.series_f{i_trial};
%             pitch_alltones(i_trial,:) = seq(1:nTones);
%         end
%         
%         % --- Initialize ---
%         if i_TD == 1
%             maxWins = 4;
%         else 
%             maxWins = 8;
%         end
%         
%         ERFdata.pitch_beta_pooled{i_sub,i_TD} = nan(nElec, maxWins);
% 
%         for i_elec = 1:nElec
%             for i_win = 1:maxWins
%                 neural_all = [];
%                 pitch_all  = [];
%                 
%                 for i_tone = 1:nTones
%                     % --- Absolute indices ---
%                     s1_abs = StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(i_tone,1);
%                     s2_abs = StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(i_tone,2);
%                     
%                     % --- Convert to relative indices ---
%                     s1 = s1_abs - StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(1,1) + 1;
%                     s2 = s2_abs - StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(1,1) + 1;
%                     
%                     % --- Define windows within this tone ---
%                     win_edges = round(linspace(s1, s2+1, maxWins+1));
%                     
%                     idx_start = win_edges(i_win);
%                     idx_end   = win_edges(i_win+1) - 1;
%                     
%                     % --- Neural data (trials × 1) ---
%                     temp   = squeeze(data(i_elec, idx_start:idx_end, :)); % samples × trials
%                     neural = squeeze(mean(temp,1)); % trials
%                     
%                     % --- Pitch (trials × 1) ---
%                     pitch = pitch_alltones(:,i_tone);
%                     
%                     % --- Stack across tones ---
%                     neural_all = [neural_all; neural(:)];
%                     pitch_all  = [pitch_all; pitch(:)];
%                 end
%                 
%                 valid_idx = ~isnan(neural_all) & ~isnan(pitch_all);
%                 
%                 if sum(valid_idx) < 20
%                     continue
%                 end
%                 
%                 x = pitch_all(valid_idx);
%                 y = neural_all(valid_idx);
%                 
%                 % --- z-score pitch ---
%                 if std(x) == 0
%                     continue
%                 end
%                 x = (x - mean(x)) ./ std(x);
%                 
%                 % --- Regression ---
%                 X = [ones(length(x),1) x];
%                 b = X \ y;
%                 
%                 % --- Store slope ---
%                 ERFdata.pitch_beta_pooled{i_sub,i_TD}(i_elec,i_win) = b(2);
%                 
%             end
%         end
%         
%         clear temp*
%     
%         %% Permutation to assess significance of regression betas
%         nPerm = 1000;
% 
%         beta_real = ERFdata.pitch_beta_pooled{i_sub,i_TD}; % elec × win
% 
%         [nElec, nWins] = size(beta_real);
% 
%         ERFdata.pitch_beta_perm{i_sub,i_TD} = nan(nElec, nWins, nPerm);
%         ERFdata.pitch_pval_perm{i_sub,i_TD} = nan(nElec, nWins);
% 
%         for i_elec = 1:nElec
%             for i_win = 1:nWins
%                 % =========================
%                 % BUILD REAL DATA (reuse logic)
%                 % =========================
%                 neural_all = [];
%                 pitch_all  = [];
% 
%                 for i_tone = 1:nTones
% 
%                     s1_abs = StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(i_tone,1);
%                     s2_abs = StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(i_tone,2);
% 
%                     s1 = s1_abs - StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(1,1) + 1;
%                     s2 = s2_abs - StimTiming.Sample_Tone_StartStop{i_sub,i_TD}(1,1) + 1;
% 
%                     win_edges = round(linspace(s1, s2+1, nWins+1));
% 
%                     idx_start = win_edges(i_win);
%                     idx_end   = win_edges(i_win+1)-1;
% 
%                     temp   = squeeze(data(i_elec, idx_start:idx_end, :));
%                     neural = squeeze(mean(temp,1)); % trials
% 
%                     pitch = pitch_alltones(:,i_tone);
% 
%                     neural_all = [neural_all; neural(:)];
%                     pitch_all  = [pitch_all; pitch(:)];
%                 end
% 
%                 valid_idx = ~isnan(neural_all) & ~isnan(pitch_all);
% 
%                 if sum(valid_idx) < 20
%                     continue
%                 end
% 
%                 y = neural_all(valid_idx);
%                 x = pitch_all(valid_idx);
% 
%                 % z-score
%                 if std(x) == 0
%                     continue
%                 end
%                 x = (x - mean(x)) ./ std(x);
% 
%                 X = [ones(length(x),1) x];
%                 b_real = X \ y;
%                 b_real = b_real(2);
% 
%                 % =========================
%                 % PERMUTATIONS
%                 % =========================
%                 beta_perm = nan(nPerm,1);
% 
%                 for perm = 1:nPerm
%                     pitch_perm_all = [];
%                     % shuffle WITHIN each trial
%                     for i_trial = 1:nTrials
%                         pitch_perm_all = [pitch_perm_all; ...
%                             pitch_alltones(i_trial, randperm(nTones))'];
%                     end
% 
%                     % match neural vector length
%                     x_perm = pitch_perm_all(valid_idx);
% 
%                     % z-score
%                     if std(x_perm) == 0
%                         continue
%                     end
%                     x_perm = (x_perm - mean(x_perm)) ./ std(x_perm);
% 
%                     % regression
%                     Xp = [ones(length(x_perm),1) x_perm];
%                     b_perm = Xp \ y;
% 
%                     beta_perm(perm) = b_perm(2);
%                 end
% 
%                 % =========================
%                 % P-VALUE
%                 % =========================
%                 pval = mean(abs(beta_perm) >= abs(b_real));
% 
%                 ERFdata.pitch_pval_perm{i_sub,i_TD}(i_elec,i_win) = pval;
%                 ERFdata.pitch_beta_perm{i_sub,i_TD}(i_elec,i_win,:) = beta_perm;
% 
%             end
%         end
%     
%         %% FDR correction across electrodes within time windows
%     
%         % skip empty entries
%         if isempty(ERFdata.pitch_pval_perm{i_sub,i_TD})
%             continue
%         end
%         
%         pval_mat = ERFdata.pitch_pval_perm{i_sub,i_TD}; % elec × win
%         [nElec, nWins] = size(pval_mat);
%         
%         sig_fdr = nan(nElec, nWins);
%         
%         % =========================
%         % LOOP WINDOWS
%         % =========================
%         for i_win = 1:nWins
%             
%             pvals = pval_mat(:,i_win);
%             
%             % valid electrodes
%             valid_idx = ~isnan(pvals);
%             
%             if sum(valid_idx) < 5
%                 continue
%             end
%             
%             % --- FDR across electrodes ---
%             [h, ~, qvals] = fdr_bh(pvals(valid_idx), 0.05);
%             
%             % reconstruct full vector
%             temp_fdr = nan(nElec,1);
%             temp_fdr(valid_idx) = h;
%             
%             sig_fdr(:,i_win) = temp_fdr;
%         end
%         
%         % =========================
%         % STORE
%         % =========================
%         ERFdata.pitch_sig_fdr{i_sub,i_TD} = sig_fdr;
        

       %% =========================================
       % 1.11 P1 REGRESSION (sample-wise, FAST)
       % =========================================
       
       Data_p1 = ERFdata.p1p33{i_sub,i_TD}; % elec × time × trial
       
       [nSensors, nSamples, nTrials] = size(Data_p1);
       
       % --- Extract pitch of p1 ---
       pitch_p1 = nan(nTrials,1);
       for i_trial = 1:nTrials
           seq = temp_DataClean_CleanTrialsElecs.behav.stim.series_f{i_trial};
           pitch_p1(i_trial) = seq(1);
       end
       
       % --- z-score pitch ---
       pitch_p1 = (pitch_p1 - nanmean(pitch_p1)) ./ nanstd(pitch_p1);
       
       % --- Design matrix ---
       X = [ones(nTrials,1), pitch_p1];   % trials × 2
       XtX_inv = inv(X' * X);             % 2 × 2
       
       % --- Parameters ---
       param.clusteralpha = 0.05;
       param.numreps      = 1000;
       
       % --- Preallocate ---
       ERFdata.P1Regression{i_sub,i_TD}.tval  = cell(nSensors,1);
       ERFdata.P1Regression{i_sub,i_TD}.pval  = cell(nSensors,1);
       ERFdata.P1Regression{i_sub,i_TD}.beta  = cell(nSensors,1);
       ERFdata.P1Regression{i_sub,i_TD}.clusterstat = cell(nSensors,1);
       
       ERFdata.P1Regression{i_sub,i_TD}.ClusterMaxStat_shuff = cell(nSensors,1);
       
       %% =========================================
       % REAL DATA
       % =========================================
       
       for i_elec = 1:nSensors
           
           % --- Get all samples at once ---
           Y = squeeze(Data_p1(i_elec,:,:)); % samples × trials
           Y = Y'; % trials × samples
           
           % --- Regression (vectorized across samples) ---
           B = XtX_inv * (X' * Y); % 2 × samples
           
           Y_hat = X * B;
           res   = Y - Y_hat;
           
           s2 = sum(res.^2) / (nTrials - 2); % 1 × samples
           se = sqrt(s2 * XtX_inv(2,2));
           
           tvals = B(2,:) ./ se;
           
           % --- p-values (two-tailed) ---
           pvals = 2 * (1 - tcdf(abs(tvals), nTrials - 2));
           
           % --- Store ---
           ERFdata.P1Regression{i_sub,i_TD}.tval{i_elec} = tvals;
           ERFdata.P1Regression{i_sub,i_TD}.pval{i_elec} = pvals;
           ERFdata.P1Regression{i_sub,i_TD}.beta{i_elec} = B(2,:);
           
           % --- CLUSTERING ---
           ERFdata.P1Regression{i_sub,i_TD}.clusterstat{i_elec} = ...
               find_temporal_clusters(tvals, pvals, param.clusteralpha);
           
           % --- Preallocate shuffle max stat ---
           ERFdata.P1Regression{i_sub,i_TD}.ClusterMaxStat_shuff{i_elec} = nan(param.numreps,1);
       end
       
       %% =========================================
       % PERMUTATION TEST (FAST)
       % =========================================
       tic
       for i_rep = 1:param.numreps
           
           ind_shuff = randperm(nTrials);
           
           % --- Shuffle once ---
           Data_shuff = Data_p1(:,:,ind_shuff);
           
           for i_elec = 1:nSensors
               
               % --- All samples at once ---
               Y = squeeze(Data_shuff(i_elec,:,:)); % samples × trials
               Y = Y'; % trials × samples
               
               % --- Regression ---
               B = XtX_inv * (X' * Y);
               
               Y_hat = X * B;
               res   = Y - Y_hat;
               
               s2 = sum(res.^2) / (nTrials - 2);
               se = sqrt(s2 * XtX_inv(2,2));
               
               tvals_shuff = B(2,:) ./ se;
               
               % --- p-values ---
               pvals_shuff = 2 * (1 - tcdf(abs(tvals_shuff), nTrials - 2));
               
               % --- CLUSTERING ---
               cluster_shuff = find_temporal_clusters( ...
                   tvals_shuff, pvals_shuff, param.clusteralpha);
               
               % --- Store only what you NEED ---
               ERFdata.P1Regression{i_sub,i_TD}.ClusterMaxStat_shuff{i_elec}(i_rep) = ...
                   cluster_shuff.maxStatSumAbs;
           end
           
           
       end
       toc
       disp(num2str(toc))
       %% =========================================
       % CLUSTER-LEVEL P-VALUES (per electrode)
       % =========================================
       
       nElec = nSensors;
       
       cluster_pvals = nan(nElec,1);
       
       for i_elec = 1:nElec
           
           real_cluster = ERFdata.P1Regression{i_sub,i_TD}.clusterstat{i_elec};
           
           % --- If no clusters → p = 1 ---
           if isempty(real_cluster) || real_cluster.nClusters == 0
               cluster_pvals(i_elec) = 1;
               continue
           end
           
           real_stat = real_cluster.maxStatSumAbs;
           
           null_dist = ERFdata.P1Regression{i_sub,i_TD}.ClusterMaxStat_shuff{i_elec};
           null_dist = null_dist(~isnan(null_dist));
           
           if isempty(null_dist)
               cluster_pvals(i_elec) = NaN;
               continue
           end
           
           % --- cluster permutation p-value ---
           cluster_pvals(i_elec) = ...
               (sum(null_dist >= real_stat) + 1) / (length(null_dist) + 1);
       end
       
       % --- Store ---
       ERFdata.P1Regression{i_sub,i_TD}.cluster_pval = cluster_pvals;
       
       %% =========================================
       % FDR CORRECTION ACROSS ELECTRODES
       % =========================================
       
       pvals = ERFdata.P1Regression{i_sub,i_TD}.cluster_pval;
       
       valid_idx = ~isnan(pvals);
       
       if sum(valid_idx) > 0
           
           [h, ~, qvals] = fdr_bh(pvals(valid_idx), 0.05);
           
           % --- reconstruct full vectors ---
           sig_fdr = false(nElec,1);
           q_full  = nan(nElec,1);
           
           sig_fdr(valid_idx) = h;
           q_full(valid_idx)  = qvals;
           
       else
           sig_fdr = nan(nElec,1);
           q_full  = nan(nElec,1);
       end
       
       % --- Store ---
       ERFdata.P1Regression{i_sub,i_TD}.sig_fdr = sig_fdr;
       ERFdata.P1Regression{i_sub,i_TD}.qvals   = q_full;
    end

    %1.8 Cleanup
    clear DataClean_AllTrials loadfile_ECoGpreprocdata
    
end
disp(['-- Preparing pre-processed data for all subjects finished after: ' ...
    num2str(round(toc/60,2)) 'min --']) %About 20 min for n = 9 per TD

%save(['ERFdata_' FuncInput_DataType '_pitchregression.mat'], 'ERFdata', '-v7.3');
save(['ERFdata_' FuncInput_DataType '_pitchregression_P1.mat'], 'ERFdata', '-v7.3');

%ERFdata=load(['ERFdata_' FuncInput_DataType '_pitchregression.mat']);
ERFdata=load(['ERFdata_' FuncInput_DataType '_pitchregression_P1.mat']);
ERFdata = ERFdata.ERFdata;

%% Identify significant electrodes for pitch regression (with at least 1 significant window)
% Loop over subjects and tone durations
% for i_sub = 1:size(ERFdata.pitch_sig_fdr,1)
%     for i_TD = 1:size(ERFdata.pitch_sig_fdr,2)
%         
%         sig_mat = ERFdata.pitch_sig_fdr{i_sub,i_TD}; % electrodes × windows
%         [nElec, nWins] = size(sig_mat);
%         
%         % Initialize
%         electrodes_with_sig = [];
%         sig_windows_per_elec = cell(nElec,1);
%         
%         for i_elec = 1:nElec
%             sig_wins = find(sig_mat(i_elec,:) > 0); % find windows with significance
%             if ~isempty(sig_wins)
%                 electrodes_with_sig(end+1) = i_elec; %#ok<AGROW>
%                 sig_windows_per_elec{i_elec} = sig_wins;
%             end
%         end
%         
%         % Display results
%         fprintf('Subject %d, TD %d:\n', i_sub, i_TD);
%         if isempty(electrodes_with_sig)
%             fprintf('  No electrodes with significant windows.\n');
%         else
%             for i = 1:length(electrodes_with_sig)
%                 elec = electrodes_with_sig(i);
%                 fprintf('  Electrode %d: significant windows = %s\n', elec, mat2str(sig_windows_per_elec{elec}));
%             end
%         end
%         fprintf('\n');
%     end
% end

%% FDR correction across time windows
% for i_sub = 1:length(subs)
%     for i_TD = 1:length(FuncInput_ToneDur_text)
% 
%         pval_mat = ERFdata.pitch_pval_perm{i_sub,i_TD}; % elec × win
%         
%         if isempty(pval_mat)
%             continue
%         end
%         
%         [nElec, nWins] = size(pval_mat);
%         
%         % =========================
%         % FLATTEN
%         % =========================
%         pvals_vec = pval_mat(:);   % (nElec*nWins) × 1
%         
%         valid_idx = ~isnan(pvals_vec);
%         
%         % =========================
%         % FDR across ALL tests
%         % =========================
%         [h, ~, qvals] = fdr_bh(pvals_vec(valid_idx), 0.05);
%         
%         % reconstruct full vectors
%         sig_vec = nan(size(pvals_vec));
%         qval_vec = nan(size(pvals_vec));
%         
%         sig_vec(valid_idx)  = h;
%         qval_vec(valid_idx) = qvals;
%         
%         % =========================
%         % RESHAPE BACK
%         % =========================
%         sig_fdr_all = reshape(sig_vec, [nElec, nWins]);
%         qvals_all   = reshape(qval_vec, [nElec, nWins]);
%         
%         % =========================
%         % STORE
%         % =========================
%         ERFdata.pitch_sig_fdr_all{i_sub,i_TD} = sig_fdr_all;
%         ERFdata.pitch_qvals_all{i_sub,i_TD}   = qvals_all;
% 
%     end
% end

%% Aggregate pvals across subjects for full pitch sequence regression
% PitchRegResults = struct('sig_fdr', [], 'pval', [], 'pval_fdr', [], 'pval_fdr_derivative', []);
% for i_sub = 1:length(subs)
%     for i_TD = 1:length(FuncInput_ToneDur_text)
%         sig_mat = ERFdata.pitch_sig_fdr{i_sub,i_TD}; % elec × win
%         pvals = ERFdata.pitch_pval_perm{i_sub,i_TD};
%         pvals_masked = pvals;
%         pvals_masked(sig_mat == 0) = NaN;
%         pvals_deriv = -(log10(pvals_masked));
%         
%         nElec = size(sig_mat,1);
%         
%         % --- Initialize padded matrix ---
%         sig_mat_8 = nan(nElec, 8);
%         pvals_8 = nan(nElec, 8);
%         pvals_masked_8 = nan(nElec, 8);
%         pvals_deriv_8 = nan(nElec, 8);
%         
%         if i_TD == 1
%             % TD1: only 4 windows → fill first 4
%             sig_mat_8(:,1:4) = sig_mat;
%             pvals_8(:,1:4) = pvals;
%             pvals_masked_8(:,1:4) = pvals_masked;
%             pvals_deriv_8(:, 1:4) = pvals_deriv; 
%         else
%             % TD2: already 8 windows
%             sig_mat_8 = sig_mat;
%             pvals_8 = pvals;
%             pvals_masked_8 = pvals_masked; 
%             pvals_deriv_8 = pvals_deriv;
%         end
%         
%         if i_sub == 1 && i_TD == 1
%             PitchRegResults.sig_fdr = sig_mat_8;
%             PitchRegResults.pval = pvals_8;
%             PitchRegResults.pval_fdr = pvals_masked_8;
%             PitchRegResults.pval_fdr_derivative = pvals_deriv_8;
%         else
%             PitchRegResults.sig_fdr = [PitchRegResults.sig_fdr; sig_mat_8];
%             PitchRegResults.pval = [PitchRegResults.pval; pvals_8];
%             PitchRegResults.pval_fdr = [PitchRegResults.pval_fdr; pvals_masked_8];
%             PitchRegResults.pval_fdr_derivative = [PitchRegResults.pval_fdr_derivative; pvals_deriv_8];
%         end
%     end
% end

%% Aggregate pvals across subjects for full P1 pitch regression
PitchRegResults = struct();
PitchRegResults.TD1 = struct('sig_fdr', [], 'pval', [], 'pval_fdr', [], 'pval_fdr_derivative', []);
PitchRegResults.TD2 = struct('sig_fdr', [], 'pval', [], 'pval_fdr', [], 'pval_fdr_derivative', []);

for i_sub = 1:length(subs)
    for i_TD = 1:length(FuncInput_ToneDur_text)
        
        sig_mat = ERFdata.P1Regression{i_sub,i_TD}.sig_fdr;
        pvals   = ERFdata.P1Regression{i_sub,i_TD}.cluster_pval;
        
        pvals_masked = pvals;
        pvals_masked(sig_mat == 0) = NaN;
        pvals_deriv = -log10(pvals);
        
        % =========================
        % PAD TO DOUBLE LENGTH
        % =========================
        nElec = length(sig_mat);
        
        sig_mat         = [sig_mat; nan(nElec,1)];
        pvals           = [pvals; nan(nElec,1)];
        pvals_masked    = [pvals_masked; nan(nElec,1)];
        pvals_deriv     = [pvals_deriv; nan(nElec,1)];
        
        % --- choose TD field ---
        if i_TD == 1
            TDfield = 'TD1';
        else
            TDfield = 'TD2';
        end
        
        % --- initialize OR append ---
        if isempty(PitchRegResults.(TDfield).sig_fdr)
            PitchRegResults.(TDfield).sig_fdr = sig_mat;
            PitchRegResults.(TDfield).pval    = pvals;
            PitchRegResults.(TDfield).pval_fdr = pvals_masked;
            PitchRegResults.(TDfield).pval_fdr_derivative = pvals_deriv;
        else
            PitchRegResults.(TDfield).sig_fdr = [PitchRegResults.(TDfield).sig_fdr; sig_mat];
            PitchRegResults.(TDfield).pval = [PitchRegResults.(TDfield).pval; pvals];
            PitchRegResults.(TDfield).pval_fdr = [PitchRegResults.(TDfield).pval_fdr; pvals_masked];
            PitchRegResults.(TDfield).pval_fdr_derivative = [PitchRegResults.(TDfield).pval_fdr_derivative; pvals_deriv];
        end
    end
end


%% 2) Load current effect data and aggregate data across subjects
clear Param2plot
coords_allsub                       = [];
labels_allsub                       = [];
anatlabels_allsub                   = [];
Param2plot.all_subs.sub_index       = [];
Param2plot.all_subs.filt_StimCorr   = [];
Param2plot.all_subs.label_AnatCat   = [];
Param2plot.all_subs.index_AnatCat   = [];
Param2plot.all_subs.TD_index        = [];
usedElecs_chanposIndex              = [];

for i_sub = 1:length(subs)
    
    tic
    sub = subs{i_sub};
    disp(['Loading data for sub: ' sub])
    NASTD_ECoG_subjectinfo %load subject info file (var: si)
    subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    for i_TD = 1:length(FuncInput_ToneDur_text)
        %Load prediction data and select current effect
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
            sub '/Data/Samplewise/'];
        load([path_inputdata ...
            sub '_PredEffectsCluster_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            FuncInput_EffectType , 'labels_loadedData');
        CurrentEffect = eval(FuncInput_EffectType);
        clear PredEffect SimplePredErrEffect ComplexPredErrEffect
        
        %Also load ECoG preproc data for channel labels and position
        loadfile_ECoGpreprocdata = [si.path_preprocdata_sub];
        load(loadfile_ECoGpreprocdata);
        
        %Determine basic parameters
        SampleFreq      = DataClean_AllTrials.fsample;
        nSensors_all    = size(CurrentEffect.stats,2);
        nSamples(i_TD)  = size(CurrentEffect.clusterstat{1}.cluster_timecourse,2);
        
        %Load stimulus correlation data and select current effect
        if strcmp(param.ElecSelect, 'StimCorr')
            path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
            load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' ...
                FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
                'corr_ttest', 'SensorLabels');
            filt_signelecs_StimCorr = corr_ttest.p < 0.05;%only sign. stim corr elec
        else
            filt_signelecs_StimCorr  = true(length(labels_loadedData),1);  %all elecs
        end
        ind_signelecs_StimCorr  = find(filt_signelecs_StimCorr);
        nSensors_sel            = sum(filt_signelecs_StimCorr);
        
        %aggregate electrode selection across subjects
        Param2plot.all_subs.filt_StimCorr = ...
            [Param2plot.all_subs.filt_StimCorr; filt_signelecs_StimCorr];
        
        %Read electrode labels, coordinates, and anatomical labels for all
        %analyzed electrodes and aggregate them across subjects
        for i_elec = 1:length(ind_signelecs_StimCorr)
            usedElecs_chanposIndex(i_elec,1) = ...
                find(strcmp(...
                labels_loadedData{ind_signelecs_StimCorr(i_elec)}, ...
                DataClean_AllTrials.elec.label));
        end
        coords_sub{i_sub, i_TD}   = ...
            DataClean_AllTrials.elec.chanpos...
            (usedElecs_chanposIndex,:);
        coords_allsub       = [coords_allsub; coords_sub{i_sub, i_TD}];
        anatlabels_allsub   = ...
            [anatlabels_allsub; ...
            DataClean_AllTrials.elec.T1AnatLabel(usedElecs_chanposIndex)];
        
        for i_elec = 1:length(ind_signelecs_StimCorr)
            labels_allsub{end+1,1} = ...
                [labels_loadedData{ind_signelecs_StimCorr(i_elec)} ' ' sub];
        end
        
        %Determine max. number of sign. clusters (uncorrected) and restrict matrix to this range
        maxnum_cluster = 10; %Across subject max cluster number estimate
        %     maxnum_cluster = 0; %Individual determination doesn't work since different across subjects
        %     for i_elec = 1:nSensors_sel
        %         temp_signcluster = [];
        %         temp_signcluster = ...
        %             sum(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval ...
        %             < param.pval_plotting);
        %         if temp_signcluster > maxnum_cluster
        %            maxnum_cluster = temp_signcluster;
        %         end
        %     end
        
        %Store cluster information from all
        Param2plot.per_sub.clusterstat{i_sub, i_TD}           = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.pval{i_sub, i_TD}                  = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.pval_derivative{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_onset{i_sub, i_TD}          = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_offset{i_sub, i_TD}         = nan(nSensors_sel, maxnum_cluster);
        Param2plot.per_sub.effect_duration{i_sub, i_TD}       = nan(nSensors_sel, maxnum_cluster);
        
        for i_elec = 1:nSensors_sel
            %Determine clusterorder based on minimal p-value
            clusterorder_minp = [];
            [~, clusterorder_minp] = ...
                sort(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval);
            %For each cluster, read out corresponding information and order it according to minpval
            index_cluster_placement = 0;
            for i_cluster = clusterorder_minp
                index_cluster_placement = index_cluster_placement + 1;
                Param2plot.per_sub.clusterstat{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %Clusterstat
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_statSum(i_cluster);
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec, index_cluster_placement) = ... %pval
                    CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster);
                %If there is a valid cluster,
                if ~isnan(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_pval(i_cluster))
                    %read out cluster timing info
                    Param2plot.per_sub.effect_onset{i_sub, i_TD}(i_elec,index_cluster_placement ) = ... %First sample of cluster to compute relative onset
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(1);
                    Param2plot.per_sub.effect_offset{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Last sample of cluster to compute relative offset
                        CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster}(end);
                    Param2plot.per_sub.effect_duration{i_sub, i_TD}(i_elec, index_cluster_placement) = ...%Relative cluster duration in samples
                        length(CurrentEffect.clusterstat{ind_signelecs_StimCorr(i_elec)}.cluster_samples{i_cluster});
                end
            end
        end
        
        %Restrict cluster data to current electrode selection (StimCorr or All)
        for i_elec = 1:nSensors_sel
            if filt_signelecs_StimCorr(ind_signelecs_StimCorr(i_elec)) == 0 %Restrict p-values to selected elecs (StimCorr or All)
                Param2plot.per_sub.pval{i_sub, i_TD}(i_elec, :) = NaN;
            end
        end
        
        %Compute derivative of p-value to estimate strength of effect
        Param2plot.per_sub.pval_derivative{i_sub, i_TD} = ...
            -(log10(Param2plot.per_sub.pval{i_sub, i_TD}));
        
        %Perform FDR_correction on cluster-p-values
        %NOTE: FDR-correction is performed across electrodes for each cluster
        %seperately, since it is much too strict if we treat different clusters
        %as different electrodes/measurements
        if param.FDRcorrect == 1
            FDR_label = 'FDRcorr';
            cluster_pvalFDR_perelec = [];
            for i_cluster = 1:size(Param2plot.per_sub.pval{i_sub, i_TD},2)
                [~, cluster_critpFDR,~ , cluster_pvalFDR_perelec(:,i_cluster)] = ...
                    fdr_bh(Param2plot.per_sub.pval{i_sub, i_TD}(:,i_cluster), param.pval_FDR, 'pdep','no');
            end
            %Replace uncorrected p-value matrix with FDR-corrected p-values
            Param2plot.per_sub.pval{i_sub, i_TD} = cluster_pvalFDR_perelec;
        else
            FDR_label = 'uncorr';
        end
        
        %Read out cluster timecourse for sign. clusters only (to later easily
        %use it for plotting sign. clusters)
        for i_elec = 1:length(Param2plot.per_sub.pval{i_sub, i_TD})
            %Set empty time course as proxy
            Param2plot.per_sub.cluster_timecourse{i_sub, i_TD}(i_elec, :) = ...
                zeros(1,length(CurrentEffect.clusterstat...
                {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse));
            
            for i_cluster = 1:length(Param2plot.per_sub.pval{i_sub, i_TD}(i_elec,:))
                %If a cluster is sign, enter its timecourse in the blank proxy
                %nd denote it with the cluster index
                if Param2plot.per_sub.pval{i_sub, i_TD}(i_elec,i_cluster) < param.pval_plotting
                    Param2plot.per_sub.cluster_timecourse{i_sub, i_TD}...
                        (i_elec, find(CurrentEffect.clusterstat...
                        {ind_signelecs_StimCorr(i_elec)}.cluster_timecourse ...
                        == i_cluster)) = i_cluster;
                end
            end
        end
        
        %Data aggregation over subjects
        %Array to differentiate subject entries
        Param2plot.all_subs.sub_index = ...
            [Param2plot.all_subs.sub_index; ...
            ones(length(coords_sub{i_sub, i_TD}),1)*i_sub];
        %Array to differentiate TD entries for each electrode
        for i_elec = 1:nSensors_all
            Param2plot.all_subs.TD_index    = ...
                [Param2plot.all_subs.TD_index; i_TD];
        end
        
        %Clusterinfo aggregation over subjects (Dim: elec * cluster)
        if i_sub == 1 && i_TD == 1
            Param2plot.all_subs.clusterstat     = ...
                Param2plot.per_sub.clusterstat{i_sub, i_TD};
            Param2plot.all_subs.pval            = ...
                Param2plot.per_sub.pval{i_sub, i_TD};
            Param2plot.all_subs.pval_derivative = ...
                Param2plot.per_sub.pval_derivative{i_sub, i_TD};
            Param2plot.all_subs.effect_onset    = ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD};
            Param2plot.all_subs.effect_offset   = ...
                Param2plot.per_sub.effect_offset{i_sub, i_TD};
            Param2plot.all_subs.effect_duration = ...
                Param2plot.per_sub.effect_duration{i_sub, i_TD};
        else
            Param2plot.all_subs.clusterstat     = ...
                [Param2plot.all_subs.clusterstat; ...
                Param2plot.per_sub.clusterstat{i_sub, i_TD}];
            Param2plot.all_subs.pval            = ...
                [Param2plot.all_subs.pval; ...
                Param2plot.per_sub.pval{i_sub, i_TD}];
            Param2plot.all_subs.pval_derivative = ...
                [Param2plot.all_subs.pval_derivative; ...
                Param2plot.per_sub.pval_derivative{i_sub, i_TD}];
            Param2plot.all_subs.effect_onset    = ...
                [Param2plot.all_subs.effect_onset; ...
                Param2plot.per_sub.effect_onset{i_sub, i_TD}];
            Param2plot.all_subs.effect_offset   = ...
                [Param2plot.all_subs.effect_offset; ...
                Param2plot.per_sub.effect_offset{i_sub, i_TD}];
            Param2plot.all_subs.effect_duration = ...
                [Param2plot.all_subs.effect_duration; ...
                Param2plot.per_sub.effect_duration{i_sub, i_TD}];
        end
        
        
        %% 3) Categorize electrodes according to anatomical regions
        AnatReg_allSubs{i_sub, i_TD} = ...
            NASTD_ECoG_AssignAnatRegions(DataClean_AllTrials, labels_loadedData);
        
        %Restrict anatomical parcellation output to selected electrodes
        %(e.g., sequence tracking selection or all elecs)
        temp_filt_allelecs_currsub = ...
            intersect(...
            find(Param2plot.all_subs.sub_index == i_sub), ...
            find(Param2plot.all_subs.TD_index == i_TD));
        temp_filt_selelecs_currsub = ...
            logical(Param2plot.all_subs.filt_StimCorr(temp_filt_allelecs_currsub));
        
        AnatReg_allSubs{i_sub, i_TD}.CatIndex = ...
            AnatReg_allSubs{i_sub, i_TD}.CatIndex(temp_filt_selelecs_currsub,:);
        AnatReg_allSubs{i_sub, i_TD}.Info_perelec = ...
            AnatReg_allSubs{i_sub, i_TD}.Info_perelec(temp_filt_selelecs_currsub,:);
        
        %Aggregate category labels and indices across subjects
        Param2plot.all_subs.label_AnatCat  = ...
            [Param2plot.all_subs.label_AnatCat; ...
            AnatReg_allSubs{i_sub, i_TD}.Info_perelec(:,3)];
        Param2plot.all_subs.index_AnatCat  = ...
            [Param2plot.all_subs.index_AnatCat; ...
            AnatReg_allSubs{i_sub, i_TD}.CatIndex];
        
        %Cleanup
        usedElecs_chanposIndex = [];
        clear CurrentEffect labels_loadedData corr_ttest DataClean_AllTrials temp*
    end
    disp(['done loading in ' num2str(toc) ' sec'])
end

%% 4) Determine sign. electrodes
%SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
%SignElecs.array         = any(PitchRegResults.pval_fdr < param.pval_plotting,2);

for i_TD = 1:2
    if i_TD == 1
        TDfield = 'TD1';
    else
        TDfield = 'TD2';
    end
        
    SignElecs.array         = PitchRegResults.(TDfield).pval < param.pval_plotting;
    SignElecs.index         = find(SignElecs.array);%1D index array denoting sign. elecs (independent of number of clusters)
    SignElecs.num_elecs     = length(SignElecs.index);
    %SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

    %Determine anatomical parcellations for sign. electrodes
    SignElecs.label_AnatCat = Param2plot.all_subs.label_AnatCat(SignElecs.index);
    SignElecs.index_AnatCat = max(Param2plot.all_subs.index_AnatCat(SignElecs.index,:),[],2);

    %Determine labels of sign. elecs and add number to subtitle
    SignElecs.labels        = labels_allsub(SignElecs.index);
    SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
    SignElecs.fulllabels    = [];
    for i_elec = 1:length(SignElecs.labels)
        SignElecs.fulllabels{i_elec,1} = ...
            [SignElecs.labels{i_elec}, ' ', SignElecs.anatlabels{i_elec}];
    end

    if ~isempty(SignElecs.labels)
        sign_title = ...
            [num2str(SignElecs.num_elecs)...
            ' / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
    else
        sign_title = ...
            ['0 / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
    end

    %Create electrode labels for plotting
    for i_elec = 1:length(labels_allsub)%No electrode labels
        labels_plotting_empty{i_elec} = '';
    end
    counter_elecs           = 0;
    labels_plotting_number  = labels_plotting_empty;
    for i_elec = SignElecs.index'
        counter_elecs = counter_elecs +1;
        labels_plotting_number{i_elec} = num2str(counter_elecs);
    end

    %% 5) Plot figure showing surf with sign. elecs and p33 ERFs per p*34 condition
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = 2;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
        
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(PitchRegResults.(TDfield).pval_fdr_derivative)
        PlotInput(i_elec,1) = max(PitchRegResults.(TDfield).pval_fdr_derivative(i_elec,:));
    end
    
    %Colorlim
    clims               = [0 max(PlotInput)];
    Label_Colorbar      = '- log10 (cluster p-value) ';
    
    dat.dimord          = 'chan_time';
    dat.time            = 0;
    dat.label           = labels_allsub;
    dat.avg             = PlotInput;
    dat.sign_elecs      = SignElecs.array;
    dat.textcolor_rgb   = [1 0 0];
    SizeFactor          = 4;
    
    chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    cmap        = 'parula';
    
    SubplotPosition = [0 0 0 0];
    
    %Project all electrodes on one hemisphere,
    coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
    
    %Plot surface
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
        chanSize, clims, cmap, dat.textcolor_rgb, ...
        DimSubplot, [1 2 DimSubplot(1)+1 DimSubplot(1)+2], SubplotPosition, [], []);
    
    sp_handle_surf = sp_handle_surf_temp.L;
    pos = sp_handle_surf.Position;   % [left bottom width height]
    sp_handle_surf.Position = [pos(1), pos(2)-0.03, pos(3), pos(4)];
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    h_title = sgtitle({ ...
    ['Input data: ' FuncInput_DataType], ...
    [TDfield '; Significant pitch tracking (sample-wise) (' ...
     num2str(SignElecs.num_elecs) ' / ' num2str(length(labels_allsub)) ' electrodes)'] ...
    }, ...
    'FontSize', 16, 'FontWeight', 'bold');
    
    if save_poststepFigs == 1
        filename     = ['Surf1H_P1pitchregression_' TDfield '_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_AllSigwindowperElec_' ...
             FuncInput_DataType '_AllTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end

%% Plot electrodes per time window
for iWin = 1:8
    num_subplot = 2;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    SignElecs.array         = any(PitchRegResults.pval_fdr(:,iWin) < param.pval_plotting,2);
    SignElecs.index         = find(SignElecs.array);%1D index array denoting sign. elecs (independent of number of clusters)
    SignElecs.num_elecs     = length(SignElecs.index);
    SignElecs.num_cluster   = sum(sum(Param2plot.all_subs.pval < param.pval_plotting));

    %Determine anatomical parcellations for sign. electrodes
    SignElecs.label_AnatCat = Param2plot.all_subs.label_AnatCat(SignElecs.index);
    SignElecs.index_AnatCat = max(Param2plot.all_subs.index_AnatCat(SignElecs.index,:),[],2);

    %Determine labels of sign. elecs and add number to subtitle
    SignElecs.labels        = labels_allsub(SignElecs.index);
    SignElecs.anatlabels    = anatlabels_allsub(SignElecs.index);
    SignElecs.fulllabels    = [];
    for i_elec = 1:length(SignElecs.labels)
        SignElecs.fulllabels{i_elec,1} = ...
            [SignElecs.labels{i_elec}, ' ', SignElecs.anatlabels{i_elec}];
    end

    if ~isempty(SignElecs.labels)
        sign_title = ...
            [num2str(SignElecs.num_elecs)...
            ' / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
    else
        sign_title = ...
            ['0 / ' num2str(length(labels_allsub)) ' sign. elecs']; %1 line
    end

    %Create electrode labels for plotting
    for i_elec = 1:length(labels_allsub)%No electrode labels
        labels_plotting_empty{i_elec} = '';
    end
    counter_elecs           = 0;
    labels_plotting_number  = labels_plotting_empty;
    for i_elec = SignElecs.index'
        counter_elecs = counter_elecs +1;
        labels_plotting_number{i_elec} = num2str(counter_elecs);
    end
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(PitchRegResults.pval_fdr_derivative)
        PlotInput(i_elec,1) = PitchRegResults.pval_fdr_derivative(i_elec,iWin);
    end
    
    %Colorlim
    clims               = [0 max(PlotInput)];
    Label_Colorbar      = '- log10 (cluster p-value) ';
    
    dat.dimord          = 'chan_time';
    dat.time            = 0;
    dat.label           = labels_allsub;
    dat.avg             = PlotInput;
    dat.sign_elecs      = SignElecs.array;
    dat.textcolor_rgb   = [1 0 0];
    SizeFactor          = 4;
    
    chanSize = ones(1,length(dat.avg))*SizeFactor; %electrode size (arbitrary)
    cmap        = 'parula';
    
    SubplotPosition = [0 0 0 0];
    
    %Project all electrodes on one hemisphere,
    coords_allsub(:,1) = abs(coords_allsub(:,1)) * -1;
    
    %Plot surface
    sp_handle_surf_temp = NASTD_ECoG_Plot_SubplotSignElecsSurf_Label_LH...
        (coords_allsub, labels_plotting_number, dat.avg, SignElecs.array,...
        chanSize, clims, cmap, dat.textcolor_rgb, ...
        DimSubplot, [1 2 DimSubplot(1)+1 DimSubplot(1)+2], SubplotPosition, [], []);
    
    sp_handle_surf = sp_handle_surf_temp.L;
    pos = sp_handle_surf.Position;   % [left bottom width height]
    sp_handle_surf.Position = [pos(1), pos(2)-0.03, pos(3), pos(4)];
    
    title(['Window ' num2str(iWin) ', ' num2str(SignElecs.num_elecs) ' sig elecs'], ...
        'FontSize', 12);
    
    h_title = sgtitle({ ...
    ['Input data: ' FuncInput_DataType], ...
    ['Significant pitch tracking'] ...
    }, ...
    'FontSize', 16, 'FontWeight', 'bold');
    
    if save_poststepFigs == 1
        filename     = ['Surf1H_pitchregression_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_Window_n_' num2str(iWin) '_'...
             FuncInput_DataType '_AllTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end