%% Project: NASTD_ECoG
%Preprocessing of raw iEEG/ECoG data:

%Load in raw ECoG data (edf-file with electrode recordings and text file 
%with electrode locations in MNI-/ T1-coordinates) and then preprocess 
%according to the following steps:

%1) Create FT-compatible struct with selected (MNI-EDF matched, strip+grid) electrodes
%2) Read out auditory triggers and define blocks and trials
%3) Cut continious recording into blocks
%4) Detrend, demean, BS-filter line noise, and HP-filter
%5) Remove pulse artifacts
%6) Compute common reference data using good electrodes
%7) Cut block data into single trials
%8) Ensure fit between iEEG and behav trial matrix
%9) Manually reject bad electrodes and trials
%10) Save clean data

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/')
%Add project base path

NASTD_ECoG_setVars
paths_NASTD_ECoG = NASTD_ECoG_paths;

% addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));
%Add project base and script dir

%% 0.2) Determine subject-specific parameters (whole-recording)
sub_list = vars.sub_list(vars.validSubjs);

i_sub = 1; %Select current subject
sub = sub_list{i_sub};
disp(['Currently preprocessing: ' sub]);

NASTD_ECoG_subjectinfo %Load subject info file (var: si)
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings; 
%Load in file with individual preproc infos

save_rawFTstruct = 0;

%Define directory for preprocessed data
preproc_dir = [paths_NASTD_ECoG.Preproc_ECoGdata '/' sub '/'];
if (~exist(preproc_dir, 'dir')); mkdir(preproc_dir); end
cd(preproc_dir);

%Plotting options
plot_poststepFigs = 1;
save_poststepFigs = 1;

%% Step 1: Create FT-compatible struct 
%with selected (MNI-EDF matched, strip+grid) electrodes
data_ECoGraw = NASTD_ECoG_Preproc_CreateFTstruct...
    (i_sub, save_rawFTstruct, vars, preproc_dir);

% 1.2) Visually inspect data
%Plot ECoG data in databrowser for coarse overview
cfg             = [];
cfg.channel     = data_ECoGraw.info_elec.selected.Label; %select Grid&Strip elecs
cfg.viewmode    = 'vertical';
ft_databrowser(cfg,data_ECoGraw);

%plot selected channel
sel_channel_index = 50;
sel_channel_label = data_ECoGraw.label(sel_channel_index);
figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
plot(data_ECoGraw.trial{1}(sel_channel_index,:)')

%ECG data
if ~isempty(data_ECoGraw.info_elec.all.index_ECGchan)
    figure
%     plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(1),:)...
%         - data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(2),:))
%     %Difference
    plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(1),:),'r-');
    hold on;
    plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_ECGchan(2),:),'b-');
    %Both
    title([sub ' - ECG data'])
end

%Trigger data
figure
plot(data_ECoGraw.trial{1}(data_ECoGraw.info_elec.all.index_triggerchan,:))
title([sub ' - Trigger data (Chan: ' data_ECoGraw.label{data_ECoGraw.info_elec.all.index_triggerchan} ')'])

%Plot electrodes on MNI volume to inspect visual fit
LH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_left.mat']);
RH = load([paths_NASTD_ECoG.FieldTrip '/template/anatomy/surface_pial_right.mat']);

figure
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
title([sub ' Electrode position on MNI standard brain' ])
ft_plot_mesh(LH.mesh);
ft_plot_mesh(RH.mesh);
% ft_plot_sens(ElecStruct_MNI, 'label','label'); %from MNI file
ft_plot_sens(data_ECoGraw.elec, 'label','label','facecolor', [0 0 1], ...
    'fontsize',8); %from combined FT-struct

view([-90 20]) %Left
%     view([-270 20]) %Right
material dull
lighting gouraud
camlight

%% Step 2: read out triggers from auditory trigger channel (whole-recording)
[data_ECoGraw.info_trigger] = ...
    NASTD_ECoG_Preproc_ReadTriggers...
    (sub, data_ECoGraw, plot_poststepFigs, save_poststepFigs);

%% Step 3: Cut continious recording into blocks (whole-recording to blocks)
%3.1) Create trialfunction in FT-format defining blocks. Block start and
%end is based on trigger definition done in previous step and stored in
%info_trigger subfield.
offsetBlock_inSec = 1; %enter if/how long additional time in sec should be cut out before/after block start/end
offsetBlock_inSamples = offsetBlock_inSec*data_ECoGraw.fsample;
block = NASTD_ECoG_Preproc_DefineBlocks(data_ECoGraw,offsetBlock_inSamples);

%3.2) cut continuous data into block
cfg = [];
cfg.trl = block;
data_ECoGraw_blocks = ft_redefinetrial(cfg, data_ECoGraw);

%Copy info_elec and info_trigger subfields into cfg.subfield (not done by FT-functions)
data_ECoGraw_blocks.cfg.info_elec = data_ECoGraw.info_elec;
data_ECoGraw_blocks.cfg.info_trigger = data_ECoGraw.info_trigger;

%3.3 Check if data is present for all sensors & timepoints in all blocks
for i_block = 1:length(data_ECoGraw_blocks.trial)
    if any(find(isnan(data_ECoGraw_blocks.trial{i_block})))
        missing_data = find(isnan(data_ECoGraw_blocks.trial{i_block}(1,:)));
        disp(['Data missing in block ' num2str(i_block) ...
            ' from timepoint ' num2str(data_ECoGraw_blocks.time{1}(missing_data(1))) ' sec'])
    end
end

clear data_ECoGraw

%3.4.1) Visually inspect trigger data for newly cut blocks
cfg             = [];
cfg.channel     = data_ECoGraw_blocks.label(subs_PreProcSettings.(sub).Trigger_chanindex); %trigger chan
cfg.ylim        = [min(data_ECoGraw_blocks.trial{1}(subs_PreProcSettings.(sub).Trigger_chanindex,:)) ...
    max(data_ECoGraw_blocks.trial{1}(subs_PreProcSettings.(sub).Trigger_chanindex,:))];
%     cfg.channel     = 68:71; %5 chan
%     cfg.ylim        = 'maxabs';
cfg.viewmode    = 'vertical';
ft_databrowser(cfg,data_ECoGraw_blocks);

%3.4.2) Plot power spectrum of specific channel for all blocks
if plot_poststepFigs == 1
    for i_chan = data_ECoGraw_blocks.cfg.info_elec.selected.index4EDF' %for all selected chans
        
        disp(['Raw power spectrum documentation - elec ' '<strong>' num2str(i_chan) '</strong> /' ...
            num2str(length(data_ECoGraw_blocks.cfg.info_elec.selected.index4EDF))]);
        
        
        figure('visible','off');
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
        Figtitle = strcat([sub ' Power-Spectrum for Channel: ' data_ECoGraw_blocks.label{i_chan} '(' num2str(i_chan) ')']);
        sgtitle(strcat(Figtitle))
        
        for i_block = 1:length(data_ECoGraw_blocks.sampleinfo)
            subplot(length(data_ECoGraw_blocks.sampleinfo)/4, ...
                length(data_ECoGraw_blocks.sampleinfo)/3,i_block);
            pwelch(data_ECoGraw_blocks.trial{i_block}(i_chan,:),[],[],[],512);
            xlim([0 150]);
            ylim([-60 60]);
            subFigtitle = strcat({['Block ' num2str(i_block)]});
            title(subFigtitle)
        end
        
        if save_poststepFigs == 1
            path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/PSD_raw/' ];
            if (~exist(path_fig, 'dir')); mkdir(path_fig); end
            
            filename     = strcat([sub '_PSD_raw_AllBlocks_Chan' ...
                data_ECoGraw_blocks.label{i_chan} '(' ...
                num2str(i_chan) ').png']);
            
            figfile      = [path_fig filename];
            saveas(gcf, [figfile], 'png'); %save png version
        end
        close
    end
end

%% Step 4: Detrend and filter (block-wise)
%4.1 Define parameters
filterOrder = 3;
bandstop_freqs =  [60 120 180 240]; %band-passs freqs (power line noise) 
highPassFreq = 0.01; %in Hz

% %4.2 Detrend/Demean, BS- & HP-filter
data_ECoGfilt_blocks = data_ECoGraw_blocks; %copy struct
for i_block = 1:length(data_ECoGraw_blocks.sampleinfo)
    data_ECoGfilt_blocks.trial{i_block} = ...
        NASTD_ECoG_Preproc_DetrendLineNoiseFilt...
        (i_block, sub, subs_PreProcSettings, ...
        filterOrder, bandstop_freqs, highPassFreq,...
        data_ECoGraw_blocks);
end

%4.3 Check if data is present for all sensors & timepoints in all blocks
for i_block = 1:length(data_ECoGfilt_blocks.trial)
    if any(find(isnan(data_ECoGfilt_blocks.trial{i_block}))) == 1
        missing_data = find(isnan(data_ECoGfilt_blocks.trial{i_block}(1,:)));
        disp(['Data missing in block ' num2str(i_block) ...
            ' from timepoint ' num2str(data_ECoGfilt_blocks.time{1}(missing_data(1))) ' sec'])
    end
end

%4.4 Plot manualy filtered data for selected electrodes
if plot_poststepFigs == 1
    tic
    for i_chan = data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF'
        
        disp(['Filtered power spectrum documentation - elec ' '<strong>' num2str(i_chan) '</strong> /' ...
            num2str(length(data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF))]);
        
        figure('visible','off');
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%for i_chan = 100 %length(data_ECoGraw_blocks.label)
        for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)
            subplot(length(data_ECoGfilt_blocks.sampleinfo)/4,length(data_ECoGfilt_blocks.sampleinfo)/3,i_block);
            pwelch(data_ECoGfilt_blocks.trial{i_block}(i_chan,:),[],[],[],512); %filtered data
            subFigtitle = strcat({['Block:' num2str(i_block)]});
            title(subFigtitle)
            xlim([0 150]);
            %         xlim([55 65]);
            ylim([-60 60]);
        end
        
        Figtitle = [sub ' Power-Spectrum after manual BS-Filter for Channel:' ...
            data_ECoGfilt_blocks.label(i_chan) '(' num2str(i_chan) ')'];
        sgtitle(strcat(Figtitle{:}))
        
        if save_poststepFigs == 1
            path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/PSD_afterFilt/' ];
            if (~exist(path_fig, 'dir')); mkdir(path_fig); end
            
            filename     = strcat([sub '_PSD_LNfilterDetrend_AllBlocks_Chan' ...
                data_ECoGfilt_blocks.label{i_chan} '(' ...
                num2str(i_chan) ').png']);
            figfile      = [path_fig filename];
            
            saveas(gcf, [figfile], 'png'); %save png version
            close
        end
    end
    toc
end

%% 5) Remove Pulse Artifact
if istrue(subs_PreProcSettings.(sub).ApplyPulseArtifactRejection)
    pulseart_dir = [preproc_dir 'Figs/Preproc/PulseArtfct/'];
    if (~exist(pulseart_dir, 'dir')); mkdir(pulseart_dir); end
    cd(pulseart_dir);
    
    % subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings;
    
    data_ECoGfiltpulse_blocks = NASTD_ECoG_Preproc_RemovePulseArtifact...
        (sub, data_ECoGfilt_blocks, subs_PreProcSettings);
    
    %visually compare post- vs. pre-corrected data
    if plot_poststepFigs == 1
        tic
        for i_chan = data_ECoGfiltpulse_blocks.cfg.info_elec.selected.index4EDF'
            figure('visible','off');
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            for i_block = 1:length(data_ECoGfiltpulse_blocks.sampleinfo)
                subplot(length(data_ECoGfiltpulse_blocks.sampleinfo)/...
                    4,length(data_ECoGfiltpulse_blocks.sampleinfo)/3,i_block);
                plot(data_ECoGfilt_blocks.trial{i_block}(i_chan,2000:12000),'r'); %uncorrected data
                hold on;
                plot(data_ECoGfiltpulse_blocks.trial{i_block}(i_chan,2000:12000),'k'); %corrected data
                subFigtitle = strcat({['Block:' num2str(i_block)]});
                title(subFigtitle)
            end
            
            Figtitle = [sub ' Pulse Artifact Correction for Channel:' ...
                data_ECoGfiltpulse_blocks.label(i_chan) '(' num2str(i_chan) ')'];
            sgtitle(strcat(Figtitle{:}))
            
            if save_poststepFigs == 1
                path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/PulseArtfct/Elec/' ];
                if (~exist(path_fig, 'dir')); mkdir(path_fig); end
                
                filename     = strcat([sub '_PulseArtfctCorr_AllBlocks_Chan' ...
                    data_ECoGfiltpulse_blocks.label{i_chan} '(' ...
                    num2str(i_chan) ').png']);
                figfile      = [path_fig filename];
                
                saveas(gcf, [figfile], 'png'); %save png version
                close
            end
        end
        toc
    end
    
    data_ECoGfilt_blocks = data_ECoGfiltpulse_blocks;
    clear data_ECoGfiltpulse_blocks
end
% Heartbeat Artifact Removal
% For each subject that had an artifact-free ECG signal recorded and a heartbeat-related artifact present in
% the ECOG data (N=10), an algorithm was applied to remove this heartbeat-aligned component without
% distorting the rest of the signal (Tal & Abeles, 2013) . First, heartbeats were detected as threshold crossings of the ECG
% signal. Then for each ECoG electrode, the signal was split into a set of heartbeat-aligned trials which had
% the duration equal to twice the median of the inter-heartbeat interval, and were centered on the time
% of the heartbeat. The trial-averaged heartbeat-evoked waveform was then low-pass filtered (zero-
% phase-shift 3 rd -order Butterworth filter at <5 Hz), with a tapered window applied (Tukey window, 10%
% cutoff). This provided a template of the artifact component that could then be removed from the ECoG
% signal, time-synched to each heartbeat, without a discontinuity arising between neighboring heartbeats.
% For those subjects without a clean ECG signal (N=4), electrodes with heartbeat-related artifacts were
% removed from analyses.
% Tal I, Abeles M. Cleaning MEG artifacts using external cues. Journal of neuroscience methods
% 217, 31-38 (2013).

%% 6) Compute common reference data (block-wise)
%Aim: Compute average signal (per sample) across all Strip&Grid ECoG electrodes
%within a block and then subtract this common average from each individual
%electrode. However, for that to work, we need to determine if all channels
%are clean (i.e., low variance of signal), in order not to mess up
%the common average.

%6.1 Set Threshold (as multiple of STD) determining which channels are excluded from CAR
numstd_thresh = subs_PreProcSettings.(sub).numSTDthresh_forCAR;

%6.2 Determine suprathresh channels, take them out, compute common average
%reference (CAR) per block based on remaining channels and then subtract that
%reference from every Strip&Grid channel per block.
%Plot images for threshold and cross-correlation pre vs post-referencing
data_ECoGfiltref_blocks = NASTD_ECoG_Preproc_CommonAvgRef...
    (sub, numstd_thresh, data_ECoGfilt_blocks, subs_PreProcSettings, ...
    paths_NASTD_ECoG, plot_poststepFigs, save_poststepFigs);

%6.3 Plot power spectrum of re-referenced data for selected electrodes
if plot_poststepFigs == 1
    tic
    for i_chan = data_ECoGfiltref_blocks.cfg.info_elec.selected.index4EDF'
        
        disp(['CAR power spectrum documentation - elec ' '<strong>' num2str(i_chan) '</strong> /' ...
            num2str(length(data_ECoGfiltref_blocks.cfg.info_elec.selected.index4EDF))]);
        
        figure('visible','off');
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
        for i_block = 1:length(data_ECoGfiltref_blocks.sampleinfo)
            
            subplot(length(data_ECoGfiltref_blocks.sampleinfo)/4,...
                length(data_ECoGfiltref_blocks.sampleinfo)/3,i_block);
            
            pwelch(data_ECoGfiltref_blocks.trial{i_block}(i_chan,:),[],[],[],512); %filtered data
            subFigtitle = strcat({['Block:' num2str(i_block)]});
            title(subFigtitle)
            xlim([0 150]);
            %         xlim([55 65]);
            ylim([-60 60]);
        end
        
        Figtitle = [sub ' Power-Spectrum after common avg referencing for Channel:' ...
            data_ECoGfiltref_blocks.label(i_chan) '(' num2str(i_chan) ')'];
        sgtitle(strcat(Figtitle{:}))
        
        if save_poststepFigs == 1
            path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/CommonAvgRef/' ];
            if (~exist(path_fig, 'dir')); mkdir(path_fig); end
            
            filename     = strcat([sub '_PSD_CAref_AllBlocks_Chan' ...
                data_ECoGfiltref_blocks.label{i_chan} '(' ...
                num2str(i_chan) ').png']);
            figfile      = [path_fig filename];
            
            saveas(gcf, [figfile], 'png'); %save png version
        end
        close
    end
    toc
end

%6.4 Check post-referencing power spectrum to identify noisy channels -
%Document noisy channels in subjectPreProcSettings.m under unsureChan_afterCAR_label

%% 7) Cut block data into single trials
%Trial definition based on above-defined auditory triggers
%7.1) Create trialfunction in FT-format defining trials. Trial start and
%end is based on trigger definition done in previous step and stored in
%info_trigger subfield.
offsetTrial_inSec       = 1; %enter how long additional time in sec should be cut out before/after block start/end
offsetTrial_inSamples   = offsetTrial_inSec*data_ECoGfiltref_blocks.fsample;
trials                  = ...
    NASTD_ECoG_Preproc_DefineTrials...
    (data_ECoGfiltref_blocks,offsetTrial_inSamples);

%7.2) Cut block data into trials
cfg                     = [];
cfg.trl                 = trials;
data_ECoGfiltref_trials = ft_redefinetrial(cfg, data_ECoGfiltref_blocks);

%Copy info_elec and info_trigger subfields into cfg.subfield (not done by FT-functions)
data_ECoGfiltref_trials.cfg.info_elec       = data_ECoGfiltref_blocks.cfg.info_elec;
data_ECoGfiltref_trials.cfg.info_trigger    = data_ECoGfiltref_blocks.cfg.info_trigger;
data_ECoGfiltref_trials.cfg.info_ref        = data_ECoGfiltref_blocks.cfg.info_ref;
%Clear old fields from cfg.previous
data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_elec');
data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_trigger');
data_ECoGfiltref_trials.cfg.previous = rmfield(data_ECoGfiltref_trials.cfg.previous,'info_ref');

%7.3 Check if data is present for all sensors & timepoints in all trials
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
    if any(find(isnan(data_ECoGfiltref_trials.trial{i_trial}))) == 1
        missing_data = find(isnan(data_ECoGfiltref_trials.trial{i_trial}(1,:)));
        disp(['Data missing in trial ' num2str(i_trial) ...
            ' from timepoint ' num2str(data_ECoGfiltref_trials.time{1}(missing_data(1))) ...
            ' sec'])
    end
end

%7.4 Plot trigger channel for all trials to visually confirm trial
%definition
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
hold on;
plot(data_ECoGfiltref_trials.trial{i_trial}...
    (subs_PreProcSettings.(sub).Trigger_chanindex,:));
title([sub ' - Trigger per trial'])
end

%% 8) Match and fuse ECoG data with Behavioral/Matlab data and remove missed trials
%Note: In ECoG recording, missed AND repeated trials are present,
%with repeated trials added at the end of the block. In Matlab/behav data,
%repeated trials are present in subfields .timing and .pahandle, but NOT in
%subfields stim and resp (there, repeated trials overwrite entries of the
%first run).
%Thus, we have to restructure the ECoG trial matrix so that it includes
%only the redone trials and not the missed trials, and that the redone
%trial-indices are in agreement with the Matlab/behav data

%8.0 Load in Matlab/behavioral data
loadfile_behav = [si.path_behavioraldata_sub];%path to indiv behavioral/stimulus data
load(loadfile_behav);

%8.1 Create binary matrix of all trials for all blocks with 1 indicating redone trials (behav index)
RedoneTrial_BehavIndex = zeros(10,12);
NumberRepeats  = zeros(10,12);
for i_redoneTrials = 1:length(data.redoneTrialNums)
    RedoneTrial_BehavIndex(data.redoneTrialNums{i_redoneTrials}) = 1;
    
    for i_numberRepeats = 1:length(data.redoneTrialNums{i_redoneTrials})
        NumberRepeats(data.redoneTrialNums{i_redoneTrials}(i_numberRepeats)) = ...
            sum(data.redoneTrialNums{i_redoneTrials} == ...
            data.redoneTrialNums{i_redoneTrials}(i_numberRepeats));
    end
end

%8.2 Create TrialIndex-per-Block matrix (ECoG data)
for i_block = unique(data_ECoGfiltref_trials.trialinfo(:,2))'
    max_trials_perblock(i_block) = ...
        size(find(data_ECoGfiltref_trials.trialinfo(:,2) == i_block),1);
end
TrialsperBlock_ECoGIndex = NaN(max(max_trials_perblock),12);
%For each iEEG trial, find out to which block and wich position within this block it belongs
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
    currentBlock = ...
        data_ECoGfiltref_trials.trialinfo(find(...
        data_ECoGfiltref_trials.trialinfo(:,1) == i_trial),2);
    TrialPos_incurrentBlock = ...
        find(find(data_ECoGfiltref_trials.trialinfo(:,2) ...
        == currentBlock) == i_trial);
    TrialsperBlock_ECoGIndex(TrialPos_incurrentBlock, currentBlock) ...
        = i_trial;
end

%8.3 For each block, replace those trials at within-block-position indicated
%by the binary matrix with those repeated at the end of the block. Note that if
%multiple trials where repeated multiple times, those trials are repeatedly
%presented in order at the end of the block.
for i_block = unique(data_ECoGfiltref_trials.trialinfo(:,2))'
    RedoneTrials = find(RedoneTrial_BehavIndex(:,i_block) ~= 0);
    
    %Mechanism: Recreate order of repeated trials for each round of repetitions. 
    %Then select the last (i.e., successfull) repetition which indexes the
    %MEG-trial corresponding to the respective behavioral trial.
    for i_redoneTrials = 1:length(RedoneTrials) %across all redone trials
        
        RepTrialStart = 11; %First possible repeated trial index per block
        RepTrialMatrix = nan(max(NumberRepeats(:,i_block)),length(RedoneTrials));
        %Empty matrix holding later indices
        for i_repcycle = 1:NumberRepeats(RedoneTrials(i_redoneTrials),i_block) 
            %for each time the current trial was repeated
            NumRepsinRepCycle(i_repcycle) = length(find(NumberRepeats(:,i_block) >= i_repcycle));
            %Check how many trials were repeated this often or more, since
            %they will be also repeated and added to this repetition cycle
            %in ascending order of trialnumber
            if i_repcycle == 1
                RepTrialMatrix(i_repcycle,:) = ...
                    [RepTrialStart: RepTrialStart + NumRepsinRepCycle(i_repcycle) - 1];
                %List repeated trial indices in ascending order of trialnumber
                RepTrialIndex = (RepTrialMatrix(i_repcycle,i_redoneTrials));
                %Select current trial based on its order in the repated
                %trials (works only for first round, since here all
                %repeated trials are still present)
            else
                LastTrialInd = RepTrialMatrix(i_repcycle-1,~isnan(RepTrialMatrix(i_repcycle-1,:)));
                %Determine with which index to continue
                RepTrialMatrix(i_repcycle,1:NumRepsinRepCycle(i_repcycle)) = ...
                    [LastTrialInd(end) + 1 : LastTrialInd(end) + NumRepsinRepCycle(i_repcycle)];
                %List repeated trial indices in ascending order of trialnumber
                CurrentIndexShift = NumRepsinRepCycle(1) - NumRepsinRepCycle(i_repcycle);
                %Difference of trial number in this compared t the first repetition cycle
                %(relevant for index placement in RepTrialMatrix)
                RepTrialIndex = (RepTrialMatrix(i_repcycle, i_redoneTrials - CurrentIndexShift));
                %Read out corresponding trial index 
            end
        end            
            
        TrialsperBlock_ECoGIndex(RedoneTrials(i_redoneTrials),i_block) = ...
            TrialsperBlock_ECoGIndex(RepTrialIndex,i_block); 
        
%         TrialsperBlock_ECoGIndex(RedoneTrials(i_redoneTrials),i_block) = ...
%             TrialsperBlock_ECoGIndex(10 + i_redoneTrials ...
%             * NumberRepeats(RedoneTrials(i_redoneTrials),i_block),i_block);        
%         TrialsperBlock_ECoGIndex(10 + i_redoneTrials,i_block) = NaN;
        
    end
end
redoneTrialIndices_perBlock = ...
    TrialsperBlock_ECoGIndex(1:10,1:12); %Remove empty rows
%trial_index = ECoG trial number, trial-in-block-position = Behav trial number

%8.4 Verify ECoG-data to Behav-data fit by comparing ToneDur for every trial
trialcounter = 0;
for i_block = 1:size(redoneTrialIndices_perBlock,2)
    for i_trial = 1:size(redoneTrialIndices_perBlock,1)
        
        ToneDur_perTrial{i_block}(i_trial,1) = ...
            length(data_ECoGfiltref_trials.trial...
            {redoneTrialIndices_perBlock(i_trial,i_block)}); %Trial length ECoG data
        
        trialcounter = trialcounter + 1;
        ToneDur_perTrial{i_block}(i_trial,2) = ...
            data.stim.toneDur(trialcounter); %ToneDur behav
        
        if ToneDur_perTrial{i_block}(i_trial,1) < 7000 ...
                && ToneDur_perTrial{i_block}(i_trial,2) > 0.2 %short ECoG TD, long Behav TD
            disp(['ToneDur Mismatch in Block: ' ...
                num2str(i_block) ', Trial: ' num2str(i_trial)])
        elseif ToneDur_perTrial{i_block}(i_trial,1) > 7000 ...
                && ToneDur_perTrial{i_block}(i_trial,2) < 0.4 %long ECoG TD, short Behav TD
            disp(['ToneDur Mismatch in Block: ' ...
                num2str(i_block) ', Trial: ' num2str(i_trial)])
        end
        
    end
end
redoneTrialIndices = [];
for i_block = 1:size(redoneTrialIndices_perBlock,2)
    redoneTrialIndices = ...
        [redoneTrialIndices; redoneTrialIndices_perBlock(:,i_block)];
end

%8.5 Re-order trial fields
%ECoG trials
NumTrials_ECoG = length(data_ECoGfiltref_trials.trial); %count ECoG trials (includes reps)
dat_fields  = fieldnames(data_ECoGfiltref_trials);
for i_field = 1:length(dat_fields)
    if eval(['length(data_ECoGfiltref_trials.' dat_fields{i_field} ') == ' num2str(NumTrials_ECoG)])
        %if missed trials are encoded in current field
        disp(['Adjusting subfield: ' dat_fields{i_field}])
        if size(eval(['data_ECoGfiltref_trials.' dat_fields{i_field}]),1) == 1
            eval(['data_ECoGfiltref_trials.' dat_fields{i_field} ...
                ' = data_ECoGfiltref_trials.' dat_fields{i_field} '(redoneTrialIndices);']);
        else
            eval(['data_ECoGfiltref_trials.' dat_fields{i_field} ...
                ' = data_ECoGfiltref_trials.' dat_fields{i_field} '(redoneTrialIndices,:);']);
        end
    end
    
end

%Behav timing files
dat_fields  = fieldnames(data.timing);
for i_field = 1:length(dat_fields)
    if eval(['length(data.timing.' dat_fields{i_field} ') == ' ...
            num2str(NumTrials_ECoG)]) %if missed trials are encoded in current field
        disp(['Adjusting subfield: ' dat_fields{i_field}])
        eval(['data.timing.' dat_fields{i_field} ' = data.timing.' ...
            dat_fields{i_field} '(redoneTrialIndices);']);
    end
end

%8.6 Append ECoG and Behav data sets
data_ECoGfiltref_trials.behav = data;

%% 9) IED detection & removal

%9.0) Set IED detection parameters
%based on Automatic spike detection from Janca et al. 2015
defaultSpikeDetectionSettings = SpikeDetectionSettings;
defaultSpikeDetectionSettings.detectionBand_k1 = 3.65; %k1 threshold value for obvious spike decision (3.65 optimium according to Janca 2015)
defaultSpikeDetectionSettings.outputAtOriginalSamplingRate = true;
% defaultSpikeDetectionSettings.combineSpikesSeconds = 0.3;

Num_IED = zeros(length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF), ...
    length(data_ECoGfiltref_trials.trial));
%Chan*Trial matrix counting the unmber of IEDs per channel-trial

%9.1) Apply automatic IED detection
tic
for i_trial = 1:length(data_ECoGfiltref_trials.sampleinfo)
    
    disp(['IED detection - trial ' '<strong>' num2str(i_trial) '</strong> /' ...
        num2str(length(data_ECoGfiltref_trials.sampleinfo))]);
    
    [spikeDetected_High{i_trial}, ... %time*channel matrix locating (index) and ordering (numerical entry) found spikes
        spikeDetected_localMax{i_trial}, ...%time (resampled)*channel matrix wit binary (0/1) entry
        spikeDetected_HighK{i_trial}, ...%1*num spikes matrix threshold values or each spike
        ~,~,~, MuscleArtifact_defaultSettings{i_trial}] = ...
        AutomatedSpikeDetection...
        (data_ECoGfiltref_trials.trial{i_trial}...
        (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF,:)', ... %Only selected iEEG channels
        data_ECoGfiltref_trials.fsample, ...
        defaultSpikeDetectionSettings);
    
    %Count number of IED (independent of duration)
    for i_chan = 1:size(spikeDetected_High{i_trial},2)
        Num_IED(i_chan,i_trial) = 0;
        for i_sample = 1:length(spikeDetected_High{i_trial})
            if i_sample == 1 && spikeDetected_High{i_trial}(1,i_chan) > 0
                Num_IED(i_chan,i_trial) = Num_IED(i_chan,i_trial) + 1;
            elseif spikeDetected_High{i_trial}(i_sample,i_chan) > 0
                if spikeDetected_High{i_trial}(i_sample,i_chan) > ...
                        spikeDetected_High{i_trial}(i_sample-1,i_chan)
                    Num_IED(i_chan,i_trial) = Num_IED(i_chan,i_trial) + 1;
                end
            end
        end
    end    
end

%9.2) Plot summary showig IED-cotamination across elecs & trials
if plot_poststepFigs == 1
    Num_samples = zeros(length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF), length(data_ECoGfiltref_trials.trial));
    for i_trial = 1:length(data_ECoGfiltref_trials.trial)
        Num_samples(:,i_trial) = length(data_ECoGfiltref_trials.trial{i_trial});
    end
    
    relNum_IEDsamples = Num_IED./Num_samples;
    relNum_IEDsamples(:,end+1) = NaN;
    relNum_IEDsamples(:,end+1) = nanmean(relNum_IEDsamples,2);
    relNum_IEDsamples(:,end+1) = relNum_IEDsamples(:,end);
    relNum_IEDsamples(:,end+1) = relNum_IEDsamples(:,end);
    relNum_IEDsamples(:,end+1) = relNum_IEDsamples(:,end);
    
    figure('visible','on'); %ensures figure doesn't pup during plotting
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
    
    h = imagesc(log10(relNum_IEDsamples));
    caxis([-5 -2.5])
    colorbar
    colormap 'parula'
    xlabel('Trial','FontSize',12)
    ylabel('Channel','FontSize',12)
    yticks(1:3:length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF))
    chan = data_ECoGfiltref_trials.label(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF);
    yticklabels(chan(1:3:end))
    title({[sub ' - Summary IED-contaminated trials (log10(# IED-samples/# all samples)'],...
        [num2str(length(find(log10(relNum_IEDsamples(:,end)) > -3))) '/' num2str(length(relNum_IEDsamples)) ...
        ' chan with samples > 1' char(8240) ' contaminated']})
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/IEDdetect/' ];
        if (~exist(path_fig, 'dir')); mkdir(path_fig); end
        filename     = strcat([sub '_Summary_IEDdetection_ChanTrials.png']);
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
    end
    close
end

% %9.3) Determine avg.  shape (amp + duration) of IEDs to know how many
% %samples to exclude - Note: Doesn't really work, IEDs too different, now
% going by conservative duration estimate
% %9.3.1 Plot all IEDs (per channel) on top of each other
% for i_chan = 1:size(spikeDetected_High{1},2)
%     All_IED_perChan = nan(sum(Num_IED(i_chan,:)),200);
%     counter_IED = 0;
%     for i_trial = 1:length(data_ECoGfiltref_trials.trial)
%         if any(unique(spikeDetected_High{i_trial}(:,i_chan)) ~= 0)
%             all_IEDs = unique(spikeDetected_High{i_trial}(:,i_chan)');
%             all_IEDs = all_IEDs(all_IEDs > 0);
%             for i_IED = all_IEDs
%                 counter_IED = counter_IED + 1;
%                 sample_IED = find(spikeDetected_High{i_trial}(:,i_chan)' == i_IED);
%                 if sample_IED(end) <= length(data_ECoGfiltref_trials.trial{i_trial})
%                     samplelength_IED = length(sample_IED);
%                     All_IED_perChan(counter_IED,1:samplelength_IED) = data_ECoGfiltref_trials.trial{i_trial}...
%                         (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan),sample_IED);
%                 else
%                     samplelength_IED = length(data_ECoGfiltref_trials.trial{i_trial}...
%                         (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan),sample_IED(1):...
%                         length(data_ECoGfiltref_trials.trial{i_trial})));
%                     All_IED_perChan(counter_IED,1:samplelength_IED) = data_ECoGfiltref_trials.trial{i_trial}...
%                         (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan),sample_IED(1):...
%                         length(data_ECoGfiltref_trials.trial{i_trial}));
%                 end
%             end
%         end
%     end
% end
% figure;
% plot(All_IED_perChan','LineWidth',1);
% hold on;
% plot(nanmean(All_IED_perChan),'k','LineWidth',3)


%9.4) Plot all IED-contaminated trials for each channel channel
if plot_poststepFigs == 1    
    Num_samples_total = 0;
    for i_trial = 1:length(data_ECoGfiltref_trials.trial)
        Num_samples_total = Num_samples_total + length(data_ECoGfiltref_trials.trial{i_trial});
    end

    for i_chan = 1:size(spikeDetected_High{1},2)
        
        Num_IED_currChan = sum(Num_IED(i_chan,:));
        Num_IEDtrials_currChan = length(find(Num_IED(i_chan,:)));
        subplot_counter = 0;
        Num_samples_IED = 0;
        
        if Num_IED_currChan > 0          
            
            disp(['IED documentation - elec ' '<strong>' num2str(i_chan) '</strong> /' ...
                num2str(length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF))]);
            
            figure('visible','off'); %ensures figure doesn't pup during plotting
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
            
            for i_trial = find(Num_IED(i_chan,:))
                subplot_counter = subplot_counter + 1;
                subplot(round(sqrt(Num_IEDtrials_currChan)),ceil(sqrt(Num_IEDtrials_currChan)),subplot_counter);
                
                plot(1:length(data_ECoGfiltref_trials.trial{i_trial}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan),:)), ...
                    data_ECoGfiltref_trials.trial{i_trial}(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan),:), 'b');
                title(['Trial ' num2str(i_trial) ' (' num2str(Num_IED(i_chan,i_trial)) ' IEDs)'])
                
                hold on; plot(zeros(1,length(data_ECoGfiltref_trials.trial{i_trial})),'k')
                
                clear i_IEDs IED*
                i_allIEDs = find(spikeDetected_High{i_trial}(:,i_chan));
                IED_startsamples = [];
                IED_estimDur_samples = 102; %200ms
                for i_IED = 1:length(i_allIEDs)
                    if i_IED == 1
                        IED_startsamples = [IED_startsamples i_allIEDs(i_IED)];
                    elseif i_allIEDs(i_IED) > i_allIEDs(i_IED-1) + 1
                        IED_startsamples = [IED_startsamples i_allIEDs(i_IED)];
                    end
                end
                
                if ~isempty(IED_startsamples)
                    for i_IED = IED_startsamples
                        IED_selchanpos = zeros(1,length(spikeDetected_High{i_trial}));
                        IED_selchan_pos(i_IED:i_IED + IED_estimDur_samples) ...
                            = max(data_ECoGfiltref_trials.trial{i_trial}(i_chan,:));
                        IED_selchan_neg(i_IED:i_IED + IED_estimDur_samples) ...
                            = max(data_ECoGfiltref_trials.trial{i_trial}(i_chan,:)) * (-1);                        
                    end
                    
                        hold on;
                        plot2 = area(1:length(IED_selchan_pos), IED_selchan_pos, 0, ...
                            'EdgeColor', [1 0 0], 'FaceColor', [1 0 0]);
                        plot2.EdgeAlpha = 0.3;
                        plot2.FaceAlpha = 0.3;
                        hold on;
                        plot3 = area(1:length(IED_selchan_neg), IED_selchan_neg, 0, ...
                            'EdgeColor', [1 0 0], 'FaceColor', [1 0 0]);
                        plot3.EdgeAlpha = 0.3;
                        plot3.FaceAlpha = 0.3;                        
                    
                end
%                 Num_samples_IED = Num_samples_IED + length(find(spikeDetected_High{i_trial}(:,i_chan))); %Only samples with detected IEDs
                Num_samples_IED = Num_samples_IED + length(find(IED_selchan_pos ~= 0)); %Contaminated samples + manually-determined post period
            end
            
%             sgtitle([ sub ' - Elec ' data_ECoGfiltref_trials.label{data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan)} ' - ' ...
%                 num2str(sum(Num_IED(i_chan,:))) ' IED-samples found in ' num2str(Num_IEDtrials_currChan) ...
%                 '/' num2str(length(data_ECoGfiltref_trials.trial)) '  trials (' ...
%                 num2str(round((Num_samples_IED/Num_samples_total) * 1000,2)) ...
%                 char(8240) ' of all samples)']); %Only samples with detected IEDs in per mill
            sgtitle([ sub ' - Elec ' data_ECoGfiltref_trials.label{data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan)} ' - ' ...
                num2str(sum(Num_IED(i_chan,:))) ' IED-samples found in ' num2str(Num_IEDtrials_currChan) ...
                '/' num2str(length(data_ECoGfiltref_trials.trial)) '  trials (' ...
                num2str(round((Num_samples_IED/Num_samples_total) * 100,2)) ...
                '% of all samples)']); %Contaminated samples + manually-determined post period in percent
            
            if save_poststepFigs == 1
                path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/IEDdetect/' ];
                if (~exist(path_fig, 'dir')); mkdir(path_fig); end                
                filename     = strcat([sub '_IEDdetection_' num2str(Num_IED_currChan) 'IED_Chan' ...
                    data_ECoGfiltref_trials.label{data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF(i_chan)} '.png']);
                figfile      = [path_fig filename];                
                saveas(gcf, [figfile], 'png'); %save png version
            end
            close
        end
    end
end

%9.4) Remove IED-contaminated samples from data
%Clean: ensure temporal match (visual check shows only partial match), duration and
%then NaN of sesctions; take out channels with too much data loss
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
    for i_chan = 1:length(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF)
        
        i_allIEDs = find(spikeDetected_High{i_trial}(:,i_chan));
        IED_startsamples = [];
        IED_estimDur_samples = 102; %200ms
        for i_IED = 1:length(i_allIEDs)
            if i_IED == 1
                IED_startsamples = [IED_startsamples i_allIEDs(i_IED)];
            elseif i_allIEDs(i_IED) > i_allIEDs(i_IED-1) + 1
                IED_startsamples = [IED_startsamples i_allIEDs(i_IED)];
            end
        end
        
        if ~isempty(IED_startsamples)
            for i_IED = IED_startsamples
                IED_samplespost = zeros(1,length(data_ECoGfiltref_trials.trial{i_trial}));
                if i_IED + IED_estimDur_samples <= length(data_ECoGfiltref_trials.trial{i_trial})                    
                    IED_samplespost(i_IED:i_IED+IED_estimDur_samples) = 1;
                else                    
                    IED_samplespost(i_IED:length(data_ECoGfiltref_trials.trial{i_trial})) = 1;
                end
            end
            data_ECoGfiltref_trials.trial{i_trial}(i_chan,logical(IED_samplespost)) = NaN;       
        end
    end
end

disp(['Finished IED detection, documentation & cleaning after: ' num2str(round(toc/60,2)) ' min'])

%% 10) Manual Artifact Rejection
subs_PreProcSettings = NASTD_ECoG_Preproc_SubPreprocSettings; %load in file with individual preproc infos

cfg             = [];
cfg.ylim        = 'maxabs';
cfg.viewmode    = 'vertical';

% %1) Check trigger signals - does trigger placement make sense (i.e.,begin+end incl. offset)?
% cfg.channel     = data_ECoGfiltref_trials.label(subs_PreProcSettings.(sub).Trigger_chanindex); %trigger chan
% ft_databrowser(cfg,data_ECoGfiltref_trials)

%2) Check all Grid&Strip ECoG channels (do suspected noisy chans stick out?)
cfg.channel     = data_ECoGfiltref_trials.label...
    (data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF); 
    %selected channels
ft_databrowser(cfg,data_ECoGfiltref_trials)

%3) Check only suspected noisy channels (are they really noisy?) and 2 random channels for comparison
cfg.channel     = [subs_PreProcSettings.(sub).unsureChan_afterCAR_index ...
    round((subs_PreProcSettings.(sub).number_ECoGchan-1).*rand(2,1)+1)'];
ft_databrowser(cfg,data_ECoGfiltref_trials)
%Enter choice of rejected channels into: NASTD_ECoG_SubPreprocSettings under rejectedChan_

%4) Check only non-noisy channels - determine bad trials and save output
cfg.channel     = setdiff(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF, ...
    subs_PreProcSettings.(sub).rejectedChan_index); %selected channels - chans prev determined as bad
cfg = ft_databrowser(cfg,data_ECoGfiltref_trials);
%Enter choice of rejected trials into: NASTD_ECoG_SubPreprocSettings under rejectedTrials

%% 11) Manual IED-cleaning (necessary for heavily contaminated channels-of-interest (i.e., NY798)
%11.1) Plot only IED-contaminated channels in databrowser, mark IED-contaminated
%samples
IEDcontam_chanind = [33 34 36 41:46 49 50 51 52 53 54 55 57 58 62 63 ...
    112 116:119 124:126 128 130:132 134:136 142];%Specify contaminated channels

cfg = [];
cfg.channel = data_ECoGfiltref_trials.label(IEDcontam_chanind);
artfct_samples = ft_databrowser(cfg, data_ECoGfiltref_trials); %manually mark IEDs

%Store info
data_ECoGfiltref_trials.cfg.info_IEDartfct.excluded_samples = ...
    artfct_samples.artfctdef.visual.artifact;
data_ECoGfiltref_trials.cfg.info_IEDartfct.corrupted_chanind = ...
    IEDcontam_chanind;

%11.2) Relace IED-contaminated samples in selected channels with NaNs
cfg = [];
cfg.channel = data_ECoGfiltref_trials.label(IEDcontam_chanind);
data_IEDchannels = ft_preprocessing(cfg, data_ECoGfiltref_trials);

cfg = [];
cfg.artfctdef.visual.artifact = artfct_samples.artfctdef.visual.artifact; %Use manually determined samples
cfg.reject = 'nan'; %replace samples with NaNs
data_cleanedIEDchannels = ft_rejectartifact(cfg,data_IEDchannels);

%Visual Check: Temporarily interpolate NaNs and plot result for visual check
cfg = [];
cfg.method = 'linear'; %'nearest','linear','spline','pchip','cubic','v5cubic'
cfg.prewindow = 0.1; %100ms
cfg.postwindow = 0.1;
data_interpolIEDchannels = ft_interpolatenan(cfg, data_cleanedIEDchannels);
figure %Plot all outputs (original, NaNed, interploated)
    i_trial = randi(length(data_ECoGfiltref_trials.trial),1);
    i_chan = randi(length(IEDcontam_chanind),1);
    plot(data_ECoGfiltref_trials.trial{i_trial}(IEDcontam_chanind(i_chan),:)','r')
    hold on;
    plot(Test2.trial{i_trial}(i_chan,:)','b')
    hold on;
    plot(Test3.trial{i_trial}(i_chan,:)','k')
    
%11.3) Fuse original and cleaned/NaNed data
%Interpolation will be done later on
for i_trial = 1:length(data_ECoGfiltref_trials.trial)
    for i_chan = 1:length(IEDcontam_chanind)
        data_ECoGfiltref_trials.trial{i_trial}(IEDcontam_chanind(i_chan),:) ...
            = data_cleanedIEDchannels.trial{i_trial}(i_chan,:);
    end
end

%11.4) Visual Check
cfg             = [];
cfg.ylim        = 'maxabs';
cfg.viewmode    = 'vertical';
% cfg.channel = data_ECoGfiltref_trials.label(IEDcontam_chanind);
cfg.channel     = setdiff(data_ECoGfiltref_trials.cfg.info_elec.selected.index4EDF, ...
    subs_PreProcSettings.(sub).rejectedChan_index); %selected channels - chans prev determined as bad
ft_databrowser(cfg,data_ECoGfiltref_trials);
%Enter choice of rejected trials into: NASTD_ECoG_SubPreprocSettings under rejectedTrials

%% 12) Save summary file to indiv. preproc path
DataClean_AllTrials = data_ECoGfiltref_trials;

savefile = [preproc_dir sub '_DataClean_AllTrials.mat'];
save(savefile, 'DataClean_AllTrials','-v7.3');

