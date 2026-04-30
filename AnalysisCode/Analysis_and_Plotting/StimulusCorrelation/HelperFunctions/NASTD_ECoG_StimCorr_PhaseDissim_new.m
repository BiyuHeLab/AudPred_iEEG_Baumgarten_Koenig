function NASTD_ECoG_StimCorr_PhaseDissim_new(subs, FuncInput_DataType, FuncInput_ToneDur_text, param, save_poststepFigs, paths_NASTD_ECoG)
% NASTD_ECoG_StimCorr_PhaseDissim
% Computes phase-dissimilarity (ITPC_identical - ITPC_similar) per electrode
% following the Luo & Poeppel / Cosyne-style method (sliding-window FFT,
% compute ITPC per f,t, average over time, sample similar trials to match identical).
%
% Usage: NASTD_ECoG_StimCorr_PhaseDissim(subs, 'HighGamma_LogAmp', {'0.3','0.6'}, param, 1, paths_NASTD_ECoG)

%% 0) Setup
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer));

path_fig = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figs/PhaseDissim/'];
if ~exist(path_fig,'dir'), mkdir(path_fig); end

% FFT window settings (published method)
%% Settings
hanning_win = 1; % apply Hann window
PhaseDissim_all = [];  % electrodes x freqs
ElecLabels_PD   = {};
Freqs_global    = [];  % will store frequency vector

%% Determine fixed FFT length across all subjects
all_lengths = [];
for i_sub = 1:length(subs)
    sub = subs{i_sub};
    NASTD_ECoG_subjectinfo; % loads si variable
    load(si.path_preprocdata_sub); % must load DataClean_AllTrials
    nTrials = length(DataClean_AllTrials.trial);
    for tr = 1:nTrials
        all_lengths(end+1) = length(DataClean_AllTrials.trial{tr}) - round(DataClean_AllTrials.fsample*0.3);
    end
end
Nfft = min(all_lengths);  % truncate all trials to this length
fprintf('Using FFT length Nfft = %d samples\n', Nfft);

%% Loop subjects
for i_sub = 1:length(subs)
    sub = subs{i_sub};
    disp(['-- Subject: ' sub ' --']);
    NASTD_ECoG_subjectinfo; 
    load(si.path_preprocdata_sub); 
    fs = DataClean_AllTrials.fsample;
    nTrials = length(DataClean_AllTrials.trial);

    % Determine significant electrodes
    t_allTDs = []; elec_labels_allTDs = {};
    for i_TD = 1:length(FuncInput_ToneDur_text)
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
             'corr_ttest','SensorLabels');
        t_allTDs = [t_allTDs; corr_ttest.t];
        elec_labels_allTDs = [elec_labels_allTDs; SensorLabels];
    end
    [unique_labels, ~, ic] = unique(elec_labels_allTDs);
    t_max = accumarray(ic, t_allTDs, [], @max);
    sig_idx = find(t_max >= 3);
    if isempty(sig_idx), continue; end

    subject_PD = []; subject_labels = {};
    
    %% Loop electrodes
    for e = 1:numel(sig_idx)
        elec_name = unique_labels{sig_idx(e)};
        all_labels = strtrim(DataClean_AllTrials.elec.label(:));
        match_idx = find(strcmp(all_labels, strtrim(elec_name)));
        if isempty(match_idx), warning('Electrode %s not found', elec_name); continue; end
        elec_idx = match_idx(1);

        % --- Precompute FFT for all trials ---
        phase_tfn_all = zeros(floor(Nfft/2)+1, nTrials);
        power_tfn_all = zeros(floor(Nfft/2)+1, nTrials);
        for tr = 1:nTrials
            series = DataClean_AllTrials.trial{tr}(elec_idx, 1:Nfft);
            series = series(:);
            if hanning_win
                series = series .* hann(Nfft);
            end
            S = fft(series, Nfft);
            S = S(1:floor(Nfft/2)+1);
            phase_tfn_all(:, tr) = angle(S);
            power_tfn_all(:, tr) = abs(S).^2;
            if tr == 1
                Freqs = (0:floor(Nfft/2)) * fs / Nfft;
                nFreqs = numel(Freqs);
            end
        end

        % --- Compute Cphase per predID / seqID ---
        C_ident_pairs = []; C_sim_pairs = [];
        predIDs_all = unique(DataClean_AllTrials.behav.stim.predID);
        for ip = 1:length(predIDs_all)
            predID = predIDs_all(ip);
            idx_pred = find(DataClean_AllTrials.behav.stim.predID == predID);
            if isempty(idx_pred), continue; end
            seqIDs = unique(DataClean_AllTrials.behav.stim.uSeqID(idx_pred));
            for iseq = 1:length(seqIDs)
                seqID = seqIDs(iseq);
                idx_ident = find(DataClean_AllTrials.behav.stim.predID==predID & ...
                                 DataClean_AllTrials.behav.stim.uSeqID==seqID);
                if numel(idx_ident)<2, continue; end
                idx_sim_all = find(DataClean_AllTrials.behav.stim.predID==predID & ...
                                   DataClean_AllTrials.behav.stim.uSeqID~=seqID);
                if isempty(idx_sim_all), continue; end

                % Compute Cphase for identical
                phase_ident = phase_tfn_all(:, idx_ident);
                Cphase_ident = (mean(cos(phase_ident),2)).^2 + (mean(sin(phase_ident),2)).^2;

                % Compute Cphase for similar
                phase_sim = phase_tfn_all(:, idx_sim_all);
                Cphase_sim = (mean(cos(phase_sim),2)).^2 + (mean(sin(phase_sim),2)).^2;

                % Average over time -> vector
                Cident_timeavg = mean(Cphase_ident, 2)'; 
                Csim_timeavg   = mean(Cphase_sim, 2)';

                C_ident_pairs = [C_ident_pairs; Cident_timeavg];
                C_sim_pairs   = [C_sim_pairs;   Csim_timeavg];
            end
        end

        % Electrode PD
        if ~isempty(C_ident_pairs) && ~isempty(C_sim_pairs)
            PD_e = mean(C_ident_pairs,1) - mean(C_sim_pairs,1);
        else
            PD_e = nan(1,nFreqs);
        end

        subject_PD = [subject_PD; PD_e];
        subject_labels = [subject_labels; elec_name];
        fprintf('Sub %s: elec %s done (%d/%d)\n', sub, elec_name, e, numel(sig_idx));
    end

    % Concatenate across subjects
    if isempty(Freqs_global)
        Freqs_global = Freqs;
    elseif ~isequal(Freqs_global, Freqs)
        error('Frequency vector mismatch across subjects');
    end
    PhaseDissim_all = [PhaseDissim_all; subject_PD];
    ElecLabels_PD = [ElecLabels_PD; subject_labels];
end

%% Plot: mean with SEM, ticks every 0.5 Hz, labels every 2 Hz, vertical lines
if isempty(PhaseDissim_all)
    warning('No phase dissimilarity data computed (PhaseDissim_all empty).');
    return;
end

figure('Color','w'); hold on;
meanPD = nanmean(PhaseDissim_all,1);
semPD  = nanstd(PhaseDissim_all,[],1) ./ sqrt(sum(~any(isnan(PhaseDissim_all),2)));

% ribbon
h = fill([Freqs_global fliplr(Freqs_global)], ...
         [meanPD+semPD fliplr(meanPD-semPD)], [0.3 0.3 1]);
set(h,'FaceAlpha',0.15,'EdgeColor','none');

% mean line
plot(Freqs_global, meanPD, 'b', 'LineWidth', 2);

% x-axis ticks every 0.1 Hz (restricted 0-15 Hz)
xticks(0:0.1:15);

% labels every 1 Hz
label_freqs = 0:1:15;
set(gca,'XTick',label_freqs, 'XTickLabel', arrayfun(@(x) num2str(x), label_freqs, 'UniformOutput', false));

% vertical lines at stimulus frequencies
xline(2.5,'--','Color',[0.7 0 0],'LineWidth',1.2);
xline(5,'--','Color',[0.7 0 0],'LineWidth',1.2);

xlabel('Frequency (Hz)','FontSize',16);
ylabel('Phase Dissimilarity','FontSize',16);
%title(sprintf('Phase dissimilarity (t >= 3), n = %d electrodes', size(PhaseDissim_all,1)), ...
 %     'FontSize',18,'FontWeight','bold');

set(gca,'FontSize',14,'LineWidth',1.2);  % increase tick label size and axis line width
xlim([0 15]);
grid on; box off;


% Save
if save_poststepFigs
    filename = sprintf('PhaseDissim_allsubs_%s_tge3.png', FuncInput_DataType);
    saveas(gcf, fullfile(path_fig, filename), 'png');
    close(gcf);
end

end