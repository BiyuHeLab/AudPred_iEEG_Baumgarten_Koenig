%% 6) Compute amplitude envelope (block level) for different frequency bands
%For each block, filter signal to get Alpha-, Beta-, and Gamma-band activity,
%then compute Hilbert-transform to get the amplitude envelope.
%Save amplitude value normalized by average amplitude across block and
%log-amplitude normalized by average log-amplitude across block.

frequency_bands = [8,12; 15,30; 30, 70; 70, 150]; %Hz [HP,LP]
frequency_band_labels = {'Alpha', 'Beta', 'Gamma', 'HighGamma'};

for i_freqbands = 1:length(frequency_bands)
    
    data_ECoGfiltref_blocks = ...
        NASTD_ECoG_Preproc_AmpEnvelope...
        (sub, frequency_bands(i_freqbands,1), frequency_bands(i_freqbands,2), ...
        frequency_band_labels{i_freqbands}, ...
        data_ECoGfiltref_blocks, ...
        paths_NASTD_ECoG, plot_poststepFigs, save_poststepFigs);
    
end

%% 8) Apply low-pass filter to base signal

filterOrder = 3;
fsample = data_ECoGraw_blocks.fsample; %Sample freq in Hz
LowPassFreq = 35;

for i_block = 1:length(data_ECoGfiltref_blocks.trial)
    
    [lpb,lpa] = butter(filterOrder,(LowPassFreq*2)/fsample,'low');
    data_ECoGfiltref_blocks.trial{i_block} = filtfilt(lpb,lpa,data_ECoGfiltref_blocks.trial{i_block}')';
    
end

%Visual check
figure;
pwelch(data_ECoGfiltref_blocks.trial{i_block}(1,:),[],[],[],512); %unfiltered data
xlim([0 60])

%9.2.1) Workaround for FreqBand-specific subfields
for i_subfields = 1:length(frequency_band_labels)
    
    FieldLabel = [frequency_band_labels{i_subfields} '_NormLogAmp'];
    
    proxystruct_blocks = data_ECoGfiltref_blocks;
    proxystruct_blocks.trial = data_ECoGfiltref_blocks.(FieldLabel);
    proxystruct_trials = ft_redefinetrial(cfg, proxystruct_blocks);
    
    data_ECoGfiltref_trials.(FieldLabel) = proxystruct_trials.trial;
    
    clear proxystruct_trials proxystruct_blocks
end

