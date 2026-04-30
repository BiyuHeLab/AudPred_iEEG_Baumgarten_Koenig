function PowerSpectrum_Subject = process_subject_power_spectra_1swindow(sub, paths_NASTD_ECoG, vars, sub_list, ToneDur_text, subs_PreProcSettings)
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
            NASTD_ECoG_FiltNaNinterp_AmpEnvel_No100msWindow(sub, ToneDur_text{i_TD}, ...
            FrequencyBands(i_freqbands,1), FrequencyBands(i_freqbands,2), FieldLabel, ...
            DataClean_CleanTrialsElecs_perTD{i_TD}, 0, 1, paths_NASTD_ECoG);
    end
end

%% 4) Compute power spectra of amplitude envelopes
fs = DataClean_CleanTrialsElecs_perTD{1}.fsample;
window = 1; % seconds
hann_flag = 1;
n = round(fs * window);
NFFT = 2^nextpow2(n);
f = fs/2 * linspace(0,1,NFFT/2+1);
nFreq = numel(f);
nBands = length(FrequencyBand_Labels);

PowerSpectrum_Subject = struct();

for i_TD = 1:length(ToneDur_text)
    for i_freqbands = 1:nBands
        FieldLabel = [FrequencyBand_Labels{i_freqbands} '_LogAmp'];
        filtData = DataClean_CleanTrialsElecs_perTD{i_TD}.(FieldLabel);

        nTrials = numel(filtData);
        nCh = size(filtData{1},1);

        minSamples = min(cellfun(@(x) size(x,2), filtData));
        nWindows = floor(minSamples / n);

        % Preallocate
        Pxx_ch = zeros(nFreq, nCh);
        count_ch = zeros(nCh,1);

        for ch = 1:nCh
            Pxx_local = zeros(nFreq,1);
            count_local = 0;

            for t = 1:nTrials
                x = detrend(filtData{t}(ch,1:minSamples));
                X = reshape(x(1:n*nWindows), n, nWindows);
                
                if hann_flag
                    X = X .* repmat(hann(n), 1, nWindows);
                end
                
                Xf = fft(X, NFFT, 1)/NFFT;
                Pxx = 2*abs(Xf(1:NFFT/2+1,:)).^2;

                Pxx_local = Pxx_local + sum(Pxx,2);
                count_local = count_local + nWindows;
            end

            Pxx_ch(:,ch) = Pxx_local;
            count_ch(ch) = count_local;
        end

        % Reduce across channels
        Pxx_mean = sum(Pxx_ch,2) / sum(count_ch);

        PowerSpectrum_Subject.bands(i_freqbands).TD(i_TD).freq = f;
        PowerSpectrum_Subject.bands(i_freqbands).TD(i_TD).Pxx  = Pxx_mean;
    end
end
end