function PowerSpectrum_Subject = process_subject_power_spectra_FFTnowindow(sub, paths_NASTD_ECoG, vars, sub_list, ToneDur_text, subs_PreProcSettings)
%PROCESS_SUBJECT_POWER_SPECTRA Compute power spectra of amplitude envelopes for a single subject
%
% Inputs:
%   sub                     - string, subject ID
%   paths_NASTD_ECoG        - struct with project paths
%   vars                    - struct with variables (including validSubjs)
%   sub_list                - cell array of all subject IDs
%   ToneDur_text            - cell array of tone durations, e.g., {'0.2', '0.4'}
%   subs_PreProcSettings     - struct with per-subject preproc info
%
% Output:
%   PowerSpectrum_Subject   - struct containing PSD per frequency band and TD

%% 0) Load subject info and preprocessed data
disp(['Processing subject: ' sub])

% Load subject-specific info
NASTD_ECoG_subjectinfo; % loads 'si'

% Load preprocessed ECoG data
loadfile_ECoGpreprocdata = si.path_preprocdata_sub;
S = load(loadfile_ECoGpreprocdata);
DataClean_AllTrials = S.DataClean_AllTrials;
clear S

SampleFreq = DataClean_AllTrials.fsample;

%% 1) Prepare input data

% 1.1) Select only clean trials
DataClean_CleanTrials = NASTD_ECoG_Preproc_SelCleanTrials(sub, DataClean_AllTrials);

% 1.2) Select valid electrodes
Index_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.index4EDF, ...
    subs_PreProcSettings.(sub).rejectedChan_index);
Label_Valid_Elecs = setdiff(DataClean_CleanTrials.cfg.info_elec.selected.Label, ...
    subs_PreProcSettings.(sub).rejectedChan_label);

% Check electrode selection
if ~isempty(find(strcmp(sort(DataClean_CleanTrials.label(Index_Valid_Elecs)), ...
        sort(Label_Valid_Elecs)) == 0))
    warning('Mismatch in electrode selection');
end

cfg = [];
cfg.channel = Label_Valid_Elecs;
DataClean_CleanTrialsElecs = ft_selectdata(cfg, DataClean_CleanTrials);

% 1.3) Select trials per Tone Duration (TD)
for i_TD = 1:length(ToneDur_text)
    DataClean_CleanTrialsElecs_perTD{i_TD} = ...
        NASTD_ECoG_Preproc_SelTrialsperTD(ToneDur_text{i_TD}, DataClean_CleanTrialsElecs);

    % Adjust info subfield
    DataClean_CleanTrialsElecs_perTD{i_TD}.info.elec    = DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_elec;
    DataClean_CleanTrialsElecs_perTD{i_TD}.info.trigger = DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_trigger;
    DataClean_CleanTrialsElecs_perTD{i_TD}.info.ref     = DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous.info_ref;

    DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous = ...
        rmfield(DataClean_CleanTrialsElecs_perTD{i_TD}.cfg.previous, {'info_elec','info_trigger','info_ref'});
end

clear DataClean_CleanTrialsElecs DataClean_CleanTrials DataClean_AllTrials

%% 2) Baseline correction and NaN trial removal
BL_win = [-0.5 0]; % prestimulus baseline

for i_TD = 1:length(ToneDur_text)
    nTrials = length(DataClean_CleanTrialsElecs_perTD{i_TD}.trial);
    nCh = length(DataClean_CleanTrialsElecs_perTD{i_TD}.label);

    for i_trial = 1:nTrials
        for i_elec = 1:nCh
            x = DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:);
            tvec = DataClean_CleanTrialsElecs_perTD{i_TD}.time{1};
            idx = find(tvec == BL_win(1)) : find(tvec == BL_win(2));
            DataClean_CleanTrialsElecs_perTD{i_TD}.trial{i_trial}(i_elec,:) = x - nanmean(x(idx));
        end
    end

    % Remove trials where some electrodes are all NaNs
    DataClean_CleanTrialsElecs_perTD{i_TD} = ...
        NASTD_ECoG_FiltNaNinterp_RemoveNaNtrials(sub, ToneDur_text{i_TD}, ...
        DataClean_CleanTrialsElecs_perTD{i_TD});
end

%% 3) Compute frequency-band amplitude envelopes
FrequencyBands = [8 12; 15 30; 30 70; 70 150]; 
FrequencyBand_Labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};

for i_TD = 1:length(ToneDur_text)
    for i_freqbands = 1:size(FrequencyBands,1)
        FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
        [DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel), ...
            DataClean_CleanTrialsElecs_perTD{i_TD}.info.FiltInterp.(FieldLabel)] = ...
            NASTD_ECoG_FiltNaNinterp_AmpEnvel(sub, ToneDur_text{i_TD}, ...
            FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, 0, 1, paths_NASTD_ECoG);
    end
end

%% 4) Compute power spectra of amplitude envelopes
fs = DataClean_CleanTrialsElecs_perTD{1}.fsample;
nBands = length(FrequencyBand_Labels);

PowerSpectrum_Subject = struct();

for i_TD = 1:length(ToneDur_text)
    for i_freqbands = 1:nBands
        
        FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
        filtData = DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel);
        
        nTrials = numel(filtData);
        nCh = size(filtData{1},1);
        
        % --- Determine minimum samples across trials ---
        minSamples = min(cellfun(@(x) size(x,2), filtData));
        
        % --- Preallocate ---
        Pxx_all = [];
        
        for ch = 1:nCh
            for t = 1:nTrials
                
                X = filtData{t}(ch, 1:minSamples);
                X = X - mean(X);   % remove DC
                
                % --- FFT-based PSD (full epoch) ---
                N = length(X);
                Xf = fft(X);
                Pxx = (abs(Xf).^2) / (fs * N);
                
                % One-sided spectrum
                Pxx = Pxx(1:floor(N/2)+1);
                f = (0:floor(N/2)) * fs / N;
                
                Pxx_all(:, end+1) = Pxx;
            end
        end
        
        % --- Average across trials and channels ---
        Pxx_mean = mean(Pxx_all, 2);
        
        PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).freq = f;
        PowerSpectrum_Subject.(sub).bands(i_freqbands).TD(i_TD).Pxx = Pxx_mean;
        
    end
end
end