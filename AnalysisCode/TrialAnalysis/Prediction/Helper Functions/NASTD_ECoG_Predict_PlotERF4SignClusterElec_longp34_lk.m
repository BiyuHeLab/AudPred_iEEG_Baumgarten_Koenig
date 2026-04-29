function NASTD_ECoG_Predict_PlotERF4SignClusterElec_longp34_lk...
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
        
        % --- Compute sampling rate and post-tone window ---
        fs = 1 / mean(diff(temp_DataClean_CleanTrialsElecs.time{1}));
        nSamples_post = round(1 * fs); % 1 second after tone

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
        
        % Base tone length (same as before)
        tone_len = length(StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(1,1): ...
            StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(1,2));
        
        % Extended length for p34 (+1s)
        p34_len = tone_len + nSamples_post;
        
        ERFdata.p33{i_sub, i_TD} = zeros(temp_nSensors, tone_len, temp_nTrials);
        ERFdata.p34{i_sub, i_TD} = nan(temp_nSensors, p34_len, temp_nTrials); % <-- now padded-ready
        ERFdata.p32{i_sub, i_TD} = ERFdata.p33{i_sub, i_TD};
        
        ERFdata.elec_labels{i_sub, i_TD} = temp_DataClean_CleanTrialsElecs.label;
        
        %copy trial wise info (i.e., all samples for selected tone for all channels and all trials into data arrays
        for i_trial = 1:temp_nTrials
            ERFdata.p33{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex,2));
            
            % --- p34 with +1s post-tone and NaN padding ---
            start_idx = StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex+1,1);
            end_idx_nominal = StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex+1,2);
            end_idx_extended = end_idx_nominal + nSamples_post;
            max_idx = size(temp_DataClean_CleanTrialsElecs.trial{i_trial}, 2);
            % Actual available end index
            end_idx_actual = min(end_idx_extended, max_idx);
            % Extract available data
            temp_segment = temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, start_idx:end_idx_actual);
            % Create NaN-padded container (preallocated size already matches)
            temp_padded = nan(temp_nSensors, p34_len);
            % Fill available portion
            temp_padded(:,1:size(temp_segment,2)) = temp_segment;
            % Store
            ERFdata.p34{i_sub, i_TD}(:,:,i_trial) = temp_padded;
            
            ERFdata.p32{i_sub, i_TD}(:, :, i_trial) = ...
                temp_DataClean_CleanTrialsElecs.trial{i_trial}(:, ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-1,1) : ...
                StimTiming.Sample_Tone_StartStop{i_sub, i_TD}(param.ToneIndex-1,2));
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

        
        clear temp*
    end
    
    %1.8 Cleanup
    clear DataClean_AllTrials loadfile_ECoGpreprocdata
    
end
disp(['-- Preparing pre-processed data for all subjects finished after: ' ...
    num2str(round(toc/60,2)) 'min --']) %About 20 min for n = 9 per TD

save(['ERFdata_' FuncInput_DataType '_longP34.mat'], 'ERFdata', '-v7.3');
ERFdata=load(['ERFdata_' FuncInput_DataType '_longP34.mat']);
ERFdata = ERFdata.ERFdata;

%% Plot p34 + 1s after tone end to look at sound-evoked activity
fs = 512; % sampling rate

nSubs = size(ERFdata.p34,1);
nTD   = size(ERFdata.p34,2);

TD_labels = {'0.2s TD','0.4s TD'};
TD_durations = [0.2 0.4]; % in seconds
colors = lines(nTD);

figure; hold on;

for i_TD = 1:nTD
    
    all_data_concat = []; % will become [time x pooled observations]
    
    for i_sub = 1:nSubs
        
        data = ERFdata.p34{i_sub, i_TD}; 
        % size: sensors x time x trials
        
        % reshape → collapse sensors + trials
        data_reshaped = reshape(data, size(data,1), size(data,2), size(data,3));
        
        % average over sensors first (optional but cleaner)
        data_sens_avg = squeeze(nanmean(data_reshaped,1)); 
        % now: time x trials
        
        % concatenate trials across subjects
        all_data_concat = [all_data_concat, data_sens_avg]; 
        % grows: time x (all trials across subs)
    end
    
    % average across all trials (and thus subjects implicitly)
    mean_timecourse = nanmean(all_data_concat,2); % time x 1
    
    % time axis (seconds)
    t = (0:length(mean_timecourse)-1) / fs;
    
    % plot
    h_plot(i_TD) = plot(t, mean_timecourse, ...
    'LineWidth', 2, 'Color', colors(i_TD,:));
    
    % --- TD-specific tone offset line ---
    xline(TD_durations(i_TD), '--', ...
        'Color', colors(i_TD,:), ...
        'LineWidth', 1.5);
    
end

tone_len1 = size(ERFdata.p33{1,1},2); % same across subs within TD
tone_len2 = size(ERFdata.p33{1,2},2); % same across subs within TD

xline(tone_len1/fs, '--', 'Tone offset', ...
    'Color', colors(1,:), 'LineWidth', 1.5, 'HandleVisibility','off');

xline(tone_len2/fs, '--', 'Tone offset', ...
    'Color', colors(2,:), 'LineWidth', 1.5, 'HandleVisibility','off');

xline(tone_len1/fs + 0.4, ':', 'Blank offset', ...
    'Color', colors(1,:), 'LineWidth', 1.5, 'HandleVisibility','off');

xline(tone_len2/fs + 0.4, ':', 'Blank offset', ...
    'Color', colors(2,:), 'LineWidth', 1.5, 'HandleVisibility','off');

xlabel('Time (s)');
ylabel('Amplitude');
yl = ylim;
ylim([yl(1), yl(2)+0.01]);
title(['p34 ERF (' FuncInput_DataType ') across subjects, sensors, trials'], ...
    'Interpreter', 'none');
legend(h_plot, TD_labels);
grid on;
filename     = ['extended_p34_ERF_' FuncInput_DataType '.png'];
figfile      = [path_fig filename];
saveas(gcf, figfile, 'png'); %save png version


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
SignElecs.array         = any(Param2plot.all_subs.pval < param.pval_plotting,2); %1D filter denoting sign. elecs (independent of number of clusters)
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

if ~isempty(SignElecs.index)
    %% 5) Plot figure showing surf with sign. elecs and p33 ERFs per p*34 condition
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
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
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder'
        subplot_counter = subplot_counter + 1;
        
        %Determine information for current electrode
        ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
        ind_curr_TD                 = Param2plot.all_subs.TD_index(SignElecs.index(i_signelec)); %TD of current elec
        ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
        ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
            find(strcmp(ERFdata.elec_labels{ind_curr_sub, ind_curr_TD}, ...
            strtok(SignElecs.labels{i_signelec},' ')));
        label_curr_elec = SignElecs.labels{i_signelec};
        
        %For current electrode, Compute AVG + STD across trials per Predp34 condition
        for i_Predp34 = 1:length(label_Predp34)
            temp_ERFdata_perp34 = [];
            for i_trial = 1:length(ERFdata.Index_Predp34{ind_curr_sub, ind_curr_TD}{i_Predp34})
                temp_ERFdata_perp34(i_trial,:) = ...
                    ERFdata.p33{ind_curr_sub, ind_curr_TD}...
                    (ind_curr_elec_currsub,:,ERFdata.Index_Predp34{ind_curr_sub, ind_curr_TD}{i_Predp34}(i_trial));
            end
            Avgp33{i_Predp34} = nanmean(temp_ERFdata_perp34,1);
            STDp33{i_Predp34} = nanstd(temp_ERFdata_perp34,1,1);
        end
        
        %Determine trace colors
        plot_param.color = ...
            {[0, 0.4470, 0.7410, 0.7], ...
            [0, 0.75, 0.75, 0.7],[0.8500, ...
            0.3250, 0.0980, 0.7]};
        
        %Plot ERF for sign  elecs
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
        
        plot(1:length(Avgp33{1}(:)),Avgp33{1}(:),...
            'color',plot_param.color{1},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{2}(:)),Avgp33{2}(:),...
            'color',plot_param.color{2},'LineWidth', 2)
        hold on;
        plot(1:length(Avgp33{3}(:)),Avgp33{3}(:),...
            'color',plot_param.color{3},'LineWidth', 2)
        axis tight
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('p33 [s]')
        
        %Highlight sign. samples/clusters
        grey  = [100 100 100]./255;
        maxY = max([max(Avgp33{1}(:)), max(Avgp33{2}(:)), max(Avgp33{3}(:))]);
        minY = min([min(Avgp33{1}(:)), min(Avgp33{2}(:)), min(Avgp33{3}(:))]);
        
        pad = abs(max([minY maxY]) * 0.05);
        ylims_per_elec{i_signelec} = [minY - pad, maxY + pad];
        ERFdata.ylims_per_elec = ylims_per_elec;
        
        sign_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,:) > 0;
        area(1:length(sign_samples), ...
            sign_samples * ...
            (max([minY maxY])+abs(max([minY maxY])*0.05)),...
            'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        if min([minY maxY]) < 0
            area(1:length(sign_samples), ...
                sign_samples * ...
                (min([minY maxY])-abs(min([minY maxY])*0.05)),...
                'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        end
        
        %Add cluster pval as text to shading
        for i_cluster = 1:sum(~isnan(...
                Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,:)))
            if Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,i_cluster) ...
                    < param.pval_plotting
                shading_startsample = ...
                    find(Param2plot.per_sub.cluster_timecourse{ind_curr_sub, ind_curr_TD}...
                    (ind_curr_elec_currsub,:) == i_cluster);
                text(shading_startsample(1),max([minY maxY] * 0.8), ...
                    [' p = ' num2str(round(...
                    Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}...
                    (ind_curr_elec_currsub, i_cluster),3))], ...
                    'FontSize',8);
            end
        end
        
        %Scale y-axis to show sign. text if present
        ylim([min([minY maxY])-+abs(max([minY maxY])*0.05) ...
            max([minY maxY])+abs(max([minY maxY])*0.05)]); %With space for sign. info
        
        %subplot title
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec ' ' FuncInput_ToneDur_text{ind_curr_TD} 's TD'],'FontSize',8,'Interpreter','none')
        
    end
    
    %Adjust figureheader/title
    switch FuncInput_EffectType
        case 'PredEffect'
            effect_text = ['Effect: Prediction (explain p33 activity by p*34);' ...
                ' p < ' num2str(param.pval_plotting) ' ' FDR_label];
        case 'SimplePredErrEffect'
            effect_text = ['Effect: Simple Prediction error' ...
                ' (explain p34 activity by absolute p33-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
        case 'ComplexPredErrEffect'
            effect_text = ['Effect: Complex Prediction error' ...
                ' (explain p34 activity by absolute p*34-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
    end
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - ' effect_text] ...
        ['Input data: ' FuncInput_DataType ', Pooled TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp33_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' ...
            FuncInput_EffectType '_' FuncInput_DataType ...
            '_AllTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    

    %% 8) Plot ERF of p33 per p32 condition
    % Collect all p32 frequencies across subjects and TDs
    all_p32_freqs = [];

    for i_sub = 1:length(subs)
        for i_TD = 1:length(FuncInput_ToneDur_text)
            all_p32_freqs = [all_p32_freqs; ...
                ERFdata.p32_freq_labels{i_sub, i_TD}(:)];
        end
    end

    % Get unique sorted frequencies
    global_p32_freqs = unique(all_p32_freqs);

    % Store
    ERFdata.global_p32_freqs = global_p32_freqs;

    n_global = length(global_p32_freqs);

    cmap_global = parula(n_global);   % or turbo, jet, etc.
    
    %% To plot all eight frequencies separately comment out below, otherwise group by 200, 400 and 800 Hz freqs
    % Define frequency groups (in Hz)
    freq_groups = [200 400 800];

    % Assign each global frequency to a group
    global_freq_group_idx = zeros(size(global_p32_freqs));

    for i = 1:length(global_p32_freqs)
        [~, idx] = min(abs(global_p32_freqs(i) - freq_groups));
        global_freq_group_idx(i) = idx; % 1=200, 2=400, 3=800
    end

    ERFdata.global_freq_group_idx = global_freq_group_idx;
    ERFdata.freq_groups = freq_groups;

    % Define fixed colors for the 3 groups (reuse your old colors)
    cmap_groups = [
        0.85, 0.1, 0.1;     % red
        0.0,  0.5, 0.0;     % dark green
        0.5,  0.0, 0.5      % purple
        ];
    
    
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
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
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder'
        subplot_counter = subplot_counter + 1;
        
        %Determine information for current electrode
        ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
        ind_curr_TD                 = Param2plot.all_subs.TD_index(SignElecs.index(i_signelec)); %TD of current elec
        ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
        ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
            find(strcmp(ERFdata.elec_labels{ind_curr_sub, ind_curr_TD}, ...
            strtok(SignElecs.labels{i_signelec},' ')));
        label_curr_elec = SignElecs.labels{i_signelec};
        
        %For current electrode, Compute AVG + STD across trials per P32 condition
        p32_freqs = ERFdata.p32_freq_labels{ind_curr_sub, ind_curr_TD};
        %n_freqs   = length(p32_freqs);
        
        n_freqs = length(freq_groups);
        
        Avgp33 = cell(n_freqs,1);
        STDp33 = cell(n_freqs,1);
        
        % Comment below if not using grouped frequencies
        trials_per_group = cell(n_freqs,1);
        for i_freq = 1:length(p32_freqs)
            % find this freq in global list
            idx_global = find(abs(ERFdata.global_p32_freqs - p32_freqs(i_freq)) < tol, 1);
            % get its group (1,2,3)
            group_idx = ERFdata.global_freq_group_idx(idx_global);
            % force column vector
            trials = ERFdata.Index_p32{ind_curr_sub, ind_curr_TD}{i_freq}(:);
            % concatenate vertically (safe)
            trials_per_group{group_idx} = [ ...
                trials_per_group{group_idx}; ...
                trials ...
                ];
        end
        for i_group = 1:n_freqs
            trials_per_group{i_group} = unique(trials_per_group{i_group});
        end
        
        for i_group = 1:n_freqs
            trials_idx = trials_per_group{i_group};
            temp = [];
            for i_trial = 1:length(trials_idx)
                temp(i_trial,:) = ERFdata.p33{ind_curr_sub, ind_curr_TD}(...
                    ind_curr_elec_currsub,:,trials_idx(i_trial));
            end
            Avgp33{i_group} = nanmean(temp,1);
            STDp33{i_group} = nanstd(temp,1,1);
        end
        
        plot_param.color = cell(n_freqs,1);
        for i_group = 1:n_freqs
            plot_param.color{i_group} = [cmap_groups(i_group,:) 0.9];
        end
        
%         for i_freq = 1:n_freqs
%             
%             trials_idx = ERFdata.Index_p32{ind_curr_sub, ind_curr_TD}{i_freq};
%             
%             temp = [];
%             for i_trial = 1:length(trials_idx)
%                 temp(i_trial,:) = ERFdata.p33{ind_curr_sub, ind_curr_TD}(...
%                     ind_curr_elec_currsub,:,trials_idx(i_trial));
%             end
%             
%             Avgp33{i_freq} = nanmean(temp,1);
%             STDp33{i_freq} = nanstd(temp,1,1);
%         end
        
        %Determine trace colors
%         plot_param.color = cell(n_freqs,1);
%         tol = 1e-6; % for float comparison
%         for i_freq = 1:n_freqs
%             % match this electrode's frequency to global list
%             idx_global = find(abs(ERFdata.global_p32_freqs - p32_freqs(i_freq)) < tol, 1);
%             if isempty(idx_global)
%                 error('p32 frequency not found in global list');
%             end
%             % assign consistent color
%             plot_param.color{i_freq} = [cmap_global(idx_global,:) 0.9];
%         end
        
        %Plot ERF for sign  elecs
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
       
%         hold on;
%         for i_freq = 1:n_freqs
%             plot(1:length(Avgp33{i_freq}), Avgp33{i_freq}, ...
%                 'color', plot_param.color{i_freq}, 'LineWidth', 2)
%         end

        hold on;
        for i_group = 1:n_freqs
            plot(1:length(Avgp33{i_group}), Avgp33{i_group}, ...
                'color', plot_param.color{i_group}, 'LineWidth', 2)
        end
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('p33 [s]')
        
        %Highlight sign. samples/clusters
        grey  = [100 100 100]./255;
        maxY = max(cellfun(@(x) max(x(:)), Avgp33));
        minY = min(cellfun(@(x) min(x(:)), Avgp33));
        
        sign_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,:) > 0;
        area(1:length(sign_samples), ...
            sign_samples * ...
            (max([minY maxY])+abs(max([minY maxY])*0.05)),...
            'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        if min([minY maxY]) < 0
            area(1:length(sign_samples), ...
                sign_samples * ...
                (min([minY maxY])-abs(min([minY maxY])*0.05)),...
                'basevalue',0,'FaceColor',grey,'FaceAlpha', 0.3, 'LineStyle','none');
        end
        
        %Add cluster pval as text to shading
        for i_cluster = 1:sum(~isnan(...
                Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,:)))
            if Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}(ind_curr_elec_currsub,i_cluster) ...
                    < param.pval_plotting
                shading_startsample = ...
                    find(Param2plot.per_sub.cluster_timecourse{ind_curr_sub, ind_curr_TD}...
                    (ind_curr_elec_currsub,:) == i_cluster);
                text(shading_startsample(1),max([minY maxY] * 0.8), ...
                    [' p = ' num2str(round(...
                    Param2plot.per_sub.pval{ind_curr_sub, ind_curr_TD}...
                    (ind_curr_elec_currsub, i_cluster),3))], ...
                    'FontSize',8);
            end
        end
        
        %Scale y-axis to show sign. text if present
        ylim(ERFdata.ylims_per_elec{i_signelec});
        
        %subplot title
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec ' ' FuncInput_ToneDur_text{ind_curr_TD} 's TD'],'FontSize',8,'Interpreter','none')
        
    end
    
%     % Create invisible axes spanning the whole figure
%     ax_legend = axes('Position',[0 0 1 1],'Visible','off');
%     hold(ax_legend, 'on');
%     
%     % Create dummy lines for legend
%     h_legend = gobjects(n_global,1);
%     for i = 1:n_global
%         h_legend(i) = plot(ax_legend, nan, nan, ...
%             'color', cmap_global(i,:), 'LineWidth', 2);
%     end
%     
%     % Create legend attached to this invisible axes
%     lgd = legend(ax_legend, h_legend, ...
%         arrayfun(@(x) [num2str(round(x)) ' Hz'], global_p32_freqs, 'UniformOutput', false), ...
%         'Location', 'eastoutside');
%     
%     lgd.Units = 'normalized';
%     pos = lgd.Position;
%     
%     pos(1) = pos(1) - 0.02;   % move left (increase if needed)
%     pos(3) = pos(3) * 0.9;    % optional: slightly narrower
%     lgd.Position = pos;

    % Create invisible axes spanning the whole figure
    ax_legend = axes('Position',[0 0 1 1],'Visible','off');
    hold(ax_legend, 'on');

    % Number of groups (should be 3)
    n_groups = length(freq_groups);

    % Create dummy lines for legend
    h_legend = gobjects(n_groups,1);
    for i = 1:n_groups
        h_legend(i) = plot(ax_legend, nan, nan, ...
            'color', cmap_groups(i,:), 'LineWidth', 2);
    end

    % Labels for groups
    legend_labels = {'Low p32', 'Medium p32', 'High p32'};

    % Create legend
    lgd = legend(ax_legend, h_legend, legend_labels, ...
        'Location', 'eastoutside');

    % Adjust position
    lgd.Units = 'normalized';
    pos = lgd.Position;

    pos(1) = pos(1) - 0.02;   % move left
    pos(3) = pos(3) * 0.9;    % slightly narrower
    lgd.Position = pos;
    
    %Adjust figureheader/title
    switch FuncInput_EffectType
        case 'PredEffect'
            effect_text = ['Effect: Explain p33 activity by p32;' ...
                ' p < ' num2str(param.pval_plotting) ' ' FDR_label];
        case 'SimplePredErrEffect'
            effect_text = ['Effect: Simple Prediction error' ...
                ' (explain p34 activity by absolute p33-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
        case 'ComplexPredErrEffect'
            effect_text = ['Effect: Complex Prediction error' ...
                ' (explain p34 activity by absolute p*34-p34 difference); p < ' ...
                num2str(param.pval_plotting) ' ' FDR_label];
    end
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - ' effect_text] ...
        ['Input data: ' FuncInput_DataType ', Pooled TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp33_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' ...
            FuncInput_EffectType '_' FuncInput_DataType ...
            'byp32_AllTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
    
       %% 9) Plot ERF of p34 for 1sec after offset of tone
    
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
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
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder'
        subplot_counter = subplot_counter + 1;
        
        %Determine information for current electrode
        ind_curr_sub                = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec)); %Sub of current elec
        ind_curr_TD                 = Param2plot.all_subs.TD_index(SignElecs.index(i_signelec)); %TD of current elec
        ind_curr_elec_allelecs      = SignElecs.index(i_signelec); %Index of current elec in across-sub elec list
        ind_curr_elec_currsub  = ... %Index of current elec in current-sub elec list
            find(strcmp(ERFdata.elec_labels{ind_curr_sub, ind_curr_TD}, ...
            strtok(SignElecs.labels{i_signelec},' ')));
        label_curr_elec = SignElecs.labels{i_signelec};
        
        %For current electrode, Compute AVG + STD across trials
       
        temp = ERFdata.p34{ind_curr_sub, ind_curr_TD}(...
                ind_curr_elec_currsub,:,:);
        Avgp34 = nanmean(temp,3);
        STDp34 = nanstd(temp,1,3);
        
        %Plot ERF for sign  elecs
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
       
        hold on;
        plot(1:length(Avgp34), Avgp34, 'LineWidth', 2)
        tone_len_sec = size(temp,2) / SampleFreq;
        if ind_curr_TD == 1
            tone_offset = 0.2;
        else
            tone_offset = 0.4;
        end
        tone_offset_x = tone_offset * SampleFreq;
        
        xline(tone_offset_x, '--', 'Tone offset', ...
            'Color', [0 0 0], ...
            'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'left');
        
        blank_offset_x = tone_offset_x + 0.4 * SampleFreq;
        
        xline(blank_offset_x, ':', 'Blank offset', ...
            'Color', [0 0 0], ...
            'LabelVerticalAlignment', 'top', ...
            'LabelHorizontalAlignment', 'left');
        
        %Change x-axis labeling to time [s]
        temp_xticks = xticks;
        temp_text = {};
        for i_xtick = 1:length(temp_xticks)
            temp_text = [temp_text num2str(round(temp_xticks(i_xtick)/SampleFreq,2))];
        end
        xticklabels(temp_text)
        xlabel('p34 [s]')
        
        maxY = max(Avgp34);
        minY = min(Avgp34);
        
        %subplot title
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec ' ' FuncInput_ToneDur_text{ind_curr_TD} 's TD'],'FontSize',8,'Interpreter','none')
        ylim([minY-0.05*abs(minY) - 0.05, maxY+0.05*abs(maxY) + 0.05]);
    end
    
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - extended p34 ERF'] ...
        ['Input data: ' FuncInput_DataType ', Pooled TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp34_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' FuncInput_DataType ...
            '_AllTD_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
    
        %% 10) Plot ERF of p34 for 1sec after offset of tone for both TDs on each plot
    
    %Determine number of subplots (surf + 1 subplot per sign. elec)
    num_subplot = SignElecs.num_elecs + 4;
    DimSubplot = [ceil(sqrt(num_subplot)), ceil(sqrt(num_subplot))];
    SubplotPos_Surf = [1 2 DimSubplot(1)+1 DimSubplot(1)+2]; %4 subplots in upper left quadrant
    SubplotPos_ERF = [3:DimSubplot(1) DimSubplot(1)+3:num_subplot];
    
    %Prepare figure
    h = figure('visible','on'); %ensures figure doesn't pop up during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %full screen
    set(gcf,'renderer','opengl');
    
    [unique_elec_names, ia] = unique( ...
    cellfun(@(x) strtok(x,' '), SignElecs.labels, 'UniformOutput', false), ...
    'stable');

    index_unique = ia;

    [~, index_anatorder] = sort(SignElecs.index_AnatCat(index_unique));
    index_anatorder_unique = index_unique(index_anatorder);

    color_TD1 = [0.8500 0.3250 0.0980]; % orange
    color_TD2 = [0 0.4470 0.7410];      % blue

    %Determine parameter to be plotted and set up plotting struct
    for i_elec = 1:length(Param2plot.all_subs.pval_derivative)
        PlotInput(i_elec,1) = max(Param2plot.all_subs.pval_derivative(i_elec,:));
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
    
    %Colorbar
    ColorbarPosition    = [sp_handle_surf.Position(1) sp_handle_surf.Position(2) 1 1];
    h = colorbar;
    h.Position(4)       = h.Position(4)*ColorbarPosition(4); %makes colorbar shorter
    h.Label.String      = Label_Colorbar;
    h.FontSize          = 12;
    caxis(clims)
    
    %Plot ERFs for each sign. electrode
    %Order sign. elecs based on anatomical parcellation
    
    legend_handles = [];
    legend_labels  = {};
    first_subplot = true;

    [~, index_anatorder] = sort(SignElecs.index_AnatCat);
    subplot_counter = 0;
    for i_signelec = index_anatorder_unique'
        
        subplot_counter = subplot_counter + 1;
        
        % --- Identify electrode ---
        ind_curr_sub = Param2plot.all_subs.sub_index(SignElecs.index(i_signelec));
        % Find all entries corresponding to this electrode
        elec_name = strtok(SignElecs.labels{i_signelec},' ');

        all_matches = find(strcmp( ...
            cellfun(@(x) strtok(x,' '), SignElecs.labels, 'UniformOutput', false), ...
            elec_name));

        sig_TDs = Param2plot.all_subs.TD_index(SignElecs.index(all_matches));
        
        elec_name = strtok(SignElecs.labels{i_signelec},' ');
        label_curr_elec = SignElecs.labels{i_signelec};
        
        subplot(DimSubplot(1),DimSubplot(2),SubplotPos_ERF(subplot_counter))
        hold on;
        
        all_vals = [];
        
        % =========================
        % LOOP OVER TDs
        % =========================
        for i_TD_plot = 1:2
            
            labels_curr = ERFdata.elec_labels{ind_curr_sub, i_TD_plot};
            ind_elec = find(strcmp(labels_curr, elec_name));
            
            if isempty(ind_elec)
                continue
            end
            
            % --- Extract data ---
            temp = ERFdata.p34{ind_curr_sub, i_TD_plot}(ind_elec,:,:);
            Avg  = nanmean(temp,3);
            
            all_vals = [all_vals; Avg(:)];
            
            % --- Color ---
            if i_TD_plot == 1
                col = color_TD1; % orange
            else
                col = color_TD2; % blue
            end
            
            % --- Line style ---
            is_significant = ismember(i_TD_plot, sig_TDs);

            if length(sig_TDs) == 2
                % Significant in BOTH TDs → both solid
                ls = '-';
                lw = 1.5;
            else
                % Only one TD significant
                if is_significant
                    ls = '-';  lw = 1.8;
                else
                    ls = '--'; lw = 1.5;
                end
            end
            
            % --- Plot ---
            h = plot(1:length(Avg), Avg, ...
                'Color', col, ...
                'LineStyle', ls, ...
                'LineWidth', lw);
            
            % --- Store legend handles ONLY ONCE ---
            if first_subplot
                legend_handles(i_TD_plot) = h;
            end
            
        end
        
        % =========================
        % TIMING LINES (BOTH TDs)
        % =========================
        
        % TD1 (orange)
        xline(0.2*SampleFreq, '--', ...
            'Color', color_TD1, ...
            'HandleVisibility','off');
        
        xline((0.2+0.4)*SampleFreq, ':', ...
            'Color', color_TD1, ...
            'HandleVisibility','off');
        
        % TD2 (blue)
        xline(0.4*SampleFreq, '--', ...
            'Color', color_TD2, ...
            'HandleVisibility','off');
        
        xline((0.4+0.4)*SampleFreq, ':', ...
            'Color', color_TD2, ...
            'HandleVisibility','off');
        
        % =========================
        % AXIS FORMATTING
        % =========================
        
        temp_xticks = xticks;
        xticklabels(arrayfun(@(x) num2str(round(x/SampleFreq,2)), ...
            temp_xticks, 'UniformOutput', false))
        
        xlabel('Time (s)')
        
        % --- Title ---
        title([labels_plotting_number{SignElecs.index(i_signelec)} ' ' ...
            AnatReg_allSubs{2}.CatLabels{SignElecs.index_AnatCat(i_signelec)} ' ' ...
            label_curr_elec ' ' FuncInput_ToneDur_text{ind_curr_TD} 's TD'], ...
            'FontSize',8,'Interpreter','none')
        
        % =========================
        % Y-LIMITS (BOTH TDs)
        % =========================
        
        if ~isempty(all_vals)
            minY = min(all_vals);
            maxY = max(all_vals);
            padding = 0.05 * (maxY - minY);
            ylim([minY - padding, maxY + padding]);
        end
        
        % =========================
        % SIGNIFICANT TIME SHADING
        % =========================

        grey = [0.6 0.6 0.6];
        
        yl = ylim; % current y-limits
        
        % --- Get electrode name ---
        elec_name = strtok(SignElecs.labels{i_signelec},' ');
        
        % --- Find all entries corresponding to this electrode ---
        all_matches = find(strcmp( ...
            cellfun(@(x) strtok(x,' '), SignElecs.labels, 'UniformOutput', false), ...
            elec_name));
        
        % --- Find which TDs this electrode is significant in ---
        sig_TDs = Param2plot.all_subs.TD_index(SignElecs.index(all_matches));
        
        % --- Loop over significant TDs individually ---
        for td_idx = 1:length(sig_TDs)
            i_sigTD = sig_TDs(td_idx);  % scalar, safe for indexing
            
            % Find this electrode in the current TD
            labels_curr = ERFdata.elec_labels{ind_curr_sub, i_sigTD};
            ind_elec = find(strcmp(labels_curr, elec_name));
            
            if isempty(ind_elec)
                continue; % skip if electrode not found
            end
            
            % Get significant samples for this electrode & TD
%             if i_signelec == 2 
%                 sig_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub, i_sigTD}(ind_elec,:) > 0; 
%                 sig_samples(1,113:180) = 1;
            sig_samples = Param2plot.per_sub.cluster_timecourse{ind_curr_sub, i_sigTD}(ind_elec,:) > 0;
            
            % --- Find contiguous segments of significance ---
            d = diff([0 sig_samples 0]);
            start_idx = find(d == 1);
            end_idx   = find(d == -1) - 1;
            
            % --- Draw shaded patches for each contiguous segment ---
            for i_seg = 1:length(start_idx)
                x1 = start_idx(i_seg);
                x2 = end_idx(i_seg);
                
                patch([x1 x2 x2 x1], ...
                    [yl(1) yl(1) yl(2) yl(2)], ...
                    grey, ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility','off');
            end
        end
    end
    
    h_td1 = plot(nan, nan, '-', 'Color', color_TD1, 'LineWidth', 1.5);
    h_td2 = plot(nan, nan, '-', 'Color', color_TD2, 'LineWidth', 1.5);

    % --- Significance legend (black solid / dashed) ---
    h_sig = plot(nan, nan, '-k',  'LineWidth', 1.5);
    h_ns  = plot(nan, nan, '--k', 'LineWidth', 1.5);

    % --- Create combined legend ---
    lgd = legend([h_td1, h_td2, h_sig, h_ns], ...
        {'TD1 (0.2 s)', 'TD2 (0.4 s)', 'sig', 'n.s.'});

    % --- Position + styling ---
    lgd.Units = 'normalized';
    lgd.Position = [0.92 0.4 0.07 0.2]; % right side of figure
    lgd.Box = 'on';
    lgd.FontSize = 12;
    lgd.LineWidth = 1;
    
    Fig_title = {['Group level (n = ' num2str(length(subs)) ') - extended p34 ERF'] ...
        ['Input data: ' FuncInput_DataType ', Both TD - Output: ' ...
        num2str(SignElecs.num_elecs) ' sign. elecs, ' num2str(SignElecs.num_cluster) ' sign. cluster']} ;
    sgtitle(Fig_title,'FontSize',18,'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['Surf1HERFp34_SignElec' param.ElecSelect ...
            '_Allsubn' num2str(length(subs)) '_' FuncInput_DataType ...
            '_BothTDplotted_Stat' FDR_label '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, figfile, 'png'); %save png version
        close all;
    end
end