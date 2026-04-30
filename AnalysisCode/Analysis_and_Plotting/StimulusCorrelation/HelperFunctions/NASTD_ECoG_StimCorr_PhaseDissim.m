function NASTD_ECoG_StimCorr_PhaseDissim(subs, FuncInput_DataType, FuncInput_ToneDur_text, param, save_poststepFigs, paths_NASTD_ECoG)

%% 0) Setup paths
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
addpath(genpath(paths_NASTD_ECoG.Freesurfer)); 

path_fig = ['/isilon/LFMI/VMdrive/Lua/NASTD/Figs/'];
if ~exist(path_fig,'dir'); mkdir(path_fig); end

%% 1) Aggregate electrodes across TDs for each subject
PhaseDissim_all = [];
ElecLabels_PD   = {};
Freqs = 0:0.1:15;
nFreqs          = length(Freqs);

parfor i_sub = 1:length(subs)
    sub = subs{i_sub};
    disp(['-- Loading data for sub: ' sub ' --']);
    NASTD_ECoG_subjectinfo; % loads si
    
    % store t-values and labels for aggregation across TDs
    t_allTDs = [];
    elec_labels_allTDs = {};
    
    for i_TD = 1:length(FuncInput_ToneDur_text)
        path_inputdata = [paths_NASTD_ECoG.ECoGdata_StimCorr '/' sub '/Data/'];
        load([path_inputdata sub '_StimCorr_' FuncInput_DataType '_' FuncInput_ToneDur_text{i_TD} 'sTD.mat'], ...
            'corr_ttest','SensorLabels');
        
        t_allTDs = [t_allTDs; corr_ttest.t];
        elec_labels_allTDs = [elec_labels_allTDs; SensorLabels];
    end
    
    % Combine across TDs: take max t-value per electrode
    [unique_labels, ~, ic] = unique(elec_labels_allTDs);
    t_max = accumarray(ic, t_allTDs, [], @max);
    
    % Select electrodes with t >= 3
    sig_idx = find(t_max >= 3);
    if isempty(sig_idx), continue; end
    
    % Load preprocessed ECoG data
    load(si.path_preprocdata_sub); % DataClean_AllTrials

    % Prepare phase dissimilarity matrix
    nTrials = length(DataClean_AllTrials.trial);
    nElecSig = length(sig_idx);
    PD = zeros(nElecSig, nFreqs);

    % Get predID and uSeqID info
    predIDs = unique(DataClean_AllTrials.behav.stim.predID);
    
    for e = 1:nElecSig
        elec_name = unique_labels{sig_idx(e)};
        elec_idx = find(strcmp(DataClean_AllTrials.elec.label, elec_name));

        phase_trials = zeros(nTrials, nFreqs);

        % 1) Compute phase for all trials
        for f = 1:nFreqs
            freq = Freqs(f);
            if freq == 0, continue; end
            [b,a] = butter(2, [freq-0.5 freq+0.5]/(DataClean_AllTrials.fsample/2));
            
            for t = 1:nTrials
                sig = filtfilt(b,a, DataClean_AllTrials.trial{t}(elec_idx,:) );
                analytic = hilbert(sig);
                phase_trials(t,f) = angle(analytic(round(end/2)));
            end
        end

        % 2) Compute ITPC for identical vs similar trials
        dissim_ident_all = [];
        dissim_sim_all   = [];
        
        for i_predID = 1:length(predIDs)
            idx_pred = find(DataClean_AllTrials.behav.stim.predID == predIDs(i_predID));
            seqIDs   = unique(DataClean_AllTrials.behav.stim.uSeqID(idx_pred));
            
            for i_seqID = 1:length(seqIDs)
                % Identical trials: same predID & uSeqID
                idx_ident = find(DataClean_AllTrials.behav.stim.predID == predIDs(i_predID) & ...
                                 DataClean_AllTrials.behav.stim.uSeqID == seqIDs(i_seqID));
                % Similar trials: same predID, different uSeqID
                idx_sim   = find(DataClean_AllTrials.behav.stim.predID == predIDs(i_predID) & ...
                                 DataClean_AllTrials.behav.stim.uSeqID ~= seqIDs(i_seqID));
                
                % Pairwise dissimilarity: identical
                for i = 1:length(idx_ident)-1
                    for j = i+1:length(idx_ident)
                        for f = 1:nFreqs
                            phase_diff = phase_trials(idx_ident(i),f) - phase_trials(idx_ident(j),f);
                            dissim_ident_all(end+1,f) = 1 - abs(mean(exp(1i*phase_diff)));
                        end
                    end
                end
                
                % Pairwise dissimilarity: similar
                for i = 1:length(idx_ident)
                    for j = 1:length(idx_sim)
                        for f = 1:nFreqs
                            phase_diff = phase_trials(idx_ident(i),f) - phase_trials(idx_sim(j),f);
                            dissim_sim_all(end+1,f) = 1 - abs(mean(exp(1i*phase_diff)));
                        end
                    end
                end
            end
        end
        
        % Average across all pairs for this electrode
        PD(e,:) = mean(dissim_ident_all,1) - mean(dissim_sim_all,1);
        disp(['Electrode ' elec_name ' done (' num2str(e) '/' num2str(nElecSig) ')']);
    end

    PhaseDissim_all = [PhaseDissim_all; PD];
    ElecLabels_PD   = [ElecLabels_PD; unique_labels(sig_idx)];
end

%% 2) Plot phase dissimilarity across frequencies
figure; hold on;

meanPD = nanmean(PhaseDissim_all,1);
semPD  = nanstd(PhaseDissim_all,[],1) ./ sqrt(size(PhaseDissim_all,1)); % SEM

% --- Light SEM ribbon ---
h = fill([Freqs fliplr(Freqs)], ...
         [meanPD+semPD fliplr(meanPD-semPD)], ...
         [0.3 0.3 1]);
set(h,'FaceAlpha',0.15,'EdgeColor','none');

% --- Mean line ---
plot(Freqs, meanPD, 'b', 'LineWidth', 2);

% --- Ticks every 0.5 Hz, labels every 2 Hz ---
set(gca, 'XTick', Freqs);                         % ticks at 0.5 Hz resolution
label_freqs = 0:2:Freqs(end);                     % only whole 2-Hz labels
set(gca, 'XTickLabel', ...
         arrayfun(@(f) sprintf('%g',f), label_freqs, 'UniformOutput', false));
set(gca,'XTick',label_freqs);                     % only label these

% --- Axes labels ---
xlabel('Frequency (Hz)');
ylabel('Phase Dissimilarity (Identical - Similar)');
title(['Phase dissimilarity across ' num2str(size(PhaseDissim_all,1)) ' electrodes with t ≥ 3']);

% --- Dashed vertical lines ---
xline(2.5, '--', 'Color', 'r', 'LineWidth', 1.2);
xline(5,   '--', 'Color', 'r', 'LineWidth', 1.2);

xlim([Freqs(1) Freqs(end)]);
grid on;
box off;

%% 3) Save figure
if save_poststepFigs
    filename = ['PhaseDissim_' FuncInput_DataType '_AllTD.png'];
    saveas(gcf, fullfile(path_fig, filename), 'png');
    close all;
end

end